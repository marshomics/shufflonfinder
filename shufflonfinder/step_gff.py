"""Step 5: Generate, merge, and filter GFF files for shufflon annotations."""

import csv
import logging
import os
import re
from collections import defaultdict
from dataclasses import dataclass, field
from io import StringIO

import pandas as pd

from .sample_sheet import Sample
from .step_flanking import FlankingRegion, parse_cds_from_gff
from .utils import ensure_dir

logger = logging.getLogger("shufflonfinder")


# ---------------------------------------------------------------------------
# GFF generation from HMM hits (generic, multi-profile)
# ---------------------------------------------------------------------------

def hmm_hits_to_gff(
    hmm_hits_df: pd.DataFrame,
    samples: list[Sample],
    output_dir: str,
) -> dict[str, str]:
    """Convert HMM hit proteins into per-sample GFF annotation files.

    Each hit becomes an ``hmm_hit`` feature (a distinct GFF type, separate
    from the Prokka CDS) with source ``shufflonfinder`` and the matching
    HMM profile recorded in the attributes.  This lets downstream window
    extraction identify HMM hit genes unambiguously and render them as
    their own annotation track.

    Args:
        hmm_hits_df: Filtered HMM hits DataFrame with columns:
                     target_name, genome, hmm_profile.
        samples: List of all Sample objects (for GFF coordinate lookup).
        output_dir: Directory to write per-sample .gff files.

    Returns:
        Dict mapping sample_id -> path to the HMM hits GFF file.
    """
    ensure_dir(output_dir)
    files = {}

    if hmm_hits_df.empty:
        return files

    sample_map = {s.sample_id: s for s in samples}

    # Collect all profiles per protein so we can write them once per gene
    profiles_per_gene = (
        hmm_hits_df.groupby(["genome", "target_name"])["hmm_profile"]
        .apply(lambda x: ";".join(sorted(set(x))))
    )

    grouped = hmm_hits_df.groupby("genome")
    for sample_id, group in grouped:
        sample = sample_map.get(sample_id)
        if sample is None:
            logger.warning("Sample %s not found, skipping GFF generation", sample_id)
            continue

        cds_map = parse_cds_from_gff(sample.gff_path)

        gff_path = os.path.join(output_dir, f"{sample_id}_hmm_hits.gff")
        written = 0
        seen = set()
        with open(gff_path, "w") as fh:
            for _, row in group.iterrows():
                protein_id = row["target_name"]
                if protein_id in seen:
                    continue
                seen.add(protein_id)

                cds = cds_map.get(protein_id)
                if cds is None:
                    continue

                profiles = profiles_per_gene.get((sample_id, protein_id), row["hmm_profile"])

                attrs = (
                    f"ID=hmm_hit_{cds.locus_tag};"
                    f"locus_tag={cds.locus_tag};"
                    f"hmm_profiles={profiles};"
                    f"Name=HMM hit: {profiles}"
                )
                fields = [
                    cds.contig,
                    "shufflonfinder",
                    "hmm_hit",
                    str(cds.start),
                    str(cds.end),
                    ".",
                    cds.strand,
                    ".",
                    attrs,
                ]
                fh.write("\t".join(fields) + "\n")
                written += 1

        files[sample_id] = gff_path
        logger.debug("Wrote HMM hits GFF for %s: %d features", sample_id, written)

    return files


# ---------------------------------------------------------------------------
# GFF generation from IR coordinates
# ---------------------------------------------------------------------------

def ir_to_gff(ir_df: pd.DataFrame, output_dir: str) -> dict[str, str]:
    """Convert inverted repeat table rows into per-sample GFF annotation files.

    Each IR produces two features (left/forward and right/reverse) with
    source='einverted' and type='inverted_repeat'.

    Args:
        ir_df: DataFrame with columns: IR_Chr, LeftIRStart, LeftIRStop,
               RightIRStart, RightIRStop, sample_id.
        output_dir: Directory to write per-sample .gff files.

    Returns:
        Dict mapping sample_id -> path to the IR GFF file.
    """
    ensure_dir(output_dir)
    files = {}

    if ir_df.empty:
        return files

    grouped = ir_df.groupby("sample_id")
    for sample_id, group in grouped:
        gff_path = os.path.join(output_dir, f"{sample_id}_ir.gff")
        with open(gff_path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            for counter, (_, row) in enumerate(group.iterrows(), start=1):
                seq_num = f"{counter:02d}"
                is_unpaired = row.get("unpaired_site", False) is True
                # Include cluster_id in attributes if available
                cluster_attr = ""
                if "cluster_id" in row.index and pd.notna(row.get("cluster_id")):
                    cluster_attr = f";cluster_id={row['cluster_id']}"
                source = "motif_search" if is_unpaired else "einverted"
                # Left (forward) IR — skip if unpaired R-site
                has_left = (
                    pd.notna(row.get("LeftIRStart"))
                    and pd.notna(row.get("LeftIRStop"))
                    and row.get("LeftIRSequence", "") != ""
                )
                has_right = (
                    pd.notna(row.get("RightIRStart"))
                    and pd.notna(row.get("RightIRStop"))
                    and row.get("RightIRSequence", "") != ""
                )
                if has_left:
                    writer.writerow([
                        row["IR_Chr"], source, "inverted_repeat",
                        int(row["LeftIRStart"]), int(row["LeftIRStop"]),
                        ".", "+", ".",
                        f"ID=inverted_repeat_{seq_num}_FOR;Name=inverted_repeat_{seq_num}_FOR{cluster_attr}",
                    ])
                if has_right:
                    writer.writerow([
                        row["IR_Chr"], source, "inverted_repeat",
                        int(row["RightIRStart"]), int(row["RightIRStop"]),
                        ".", "-", ".",
                        f"ID=inverted_repeat_{seq_num}_REV;Name=inverted_repeat_{seq_num}_REV{cluster_attr}",
                    ])
        files[sample_id] = gff_path
        logger.debug("Wrote IR GFF for %s: %d IR entries", sample_id, len(group))

    return files


# ---------------------------------------------------------------------------
# Merge annotations into Prokka GFF
# ---------------------------------------------------------------------------

def merge_gff_into_prokka(
    prokka_gff: str,
    extra_gff_files: list[str],
    output_gff: str,
) -> str:
    """Insert additional GFF features into a Prokka GFF file.

    Prokka GFFs have a ##FASTA section at the bottom. This function inserts
    the extra annotations between the Prokka feature lines and the FASTA block.

    Args:
        prokka_gff: Path to the original Prokka .gff file.
        extra_gff_files: List of paths to GFF files to merge in.
        output_gff: Where to write the merged GFF.

    Returns:
        Path to the merged GFF file.
    """
    ensure_dir(os.path.dirname(output_gff))

    annotation_lines = []
    fasta_lines = []
    in_fasta = False

    with open(prokka_gff) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                in_fasta = True
                fasta_lines.append(line)
            elif in_fasta:
                fasta_lines.append(line)
            else:
                annotation_lines.append(line)

    # Read extra features
    extra_lines = []
    for gff_path in extra_gff_files:
        if not os.path.isfile(gff_path):
            continue
        with open(gff_path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    extra_lines.append(line)

    with open(output_gff, "w") as fh:
        fh.writelines(annotation_lines)
        if extra_lines:
            # Ensure newline separation
            if annotation_lines and not annotation_lines[-1].endswith("\n"):
                fh.write("\n")
            fh.writelines(extra_lines)
        fh.writelines(fasta_lines)

    logger.debug("Merged GFF -> %s (%d extra features)", output_gff, len(extra_lines))
    return output_gff


# ---------------------------------------------------------------------------
# Window-based GFF filtering (IR + CDS co-location)
# ---------------------------------------------------------------------------

def _parse_gff_attributes(attr_string: str) -> dict[str, str]:
    """Parse a GFF3 attributes column into a dict.

    Handles semicolon-separated key=value pairs.  Returns a dict with
    lowercased keys for case-insensitive lookup.
    """
    attrs = {}
    for part in attr_string.split(";"):
        part = part.strip()
        if "=" in part:
            key, _, val = part.partition("=")
            attrs[key] = val
            attrs[key.lower()] = val   # also store lowercase for easy lookup
    return attrs


class Feature:
    """Lightweight representation of a GFF feature for spatial operations."""

    __slots__ = ("seqid", "source", "type", "start", "end", "strand", "attributes")

    def __init__(self, seqid, source, type, start, end, strand, attributes):
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.strand = strand
        self.attributes = attributes

    def overlaps(self, other):
        return self.seqid == other.seqid and self.start <= other.end and self.end >= other.start

    def distance_to(self, other):
        if self.seqid != other.seqid:
            return float("inf")
        if self.overlaps(other):
            return 0
        return min(abs(self.start - other.end), abs(self.end - other.start))

    @property
    def parsed_attributes(self) -> dict[str, str]:
        return _parse_gff_attributes(self.attributes)

    @property
    def locus_tag(self) -> str:
        attrs = self.parsed_attributes
        return attrs.get("ID", attrs.get("locus_tag", ""))

    @property
    def product(self) -> str:
        attrs = self.parsed_attributes
        return attrs.get("product", attrs.get("Product", ""))

    @property
    def name(self) -> str:
        attrs = self.parsed_attributes
        return attrs.get("Name", attrs.get("name", ""))


@dataclass
class ShufflonWindow:
    """Structured representation of one extracted shufflon window."""
    sample_id: str
    window_id: str
    contig: str
    window_start: int   # 0-based
    window_end: int     # 1-based
    n_ir_pairs: int
    ir_features: list[Feature] = field(default_factory=list)
    cds_features: list[Feature] = field(default_factory=list)
    hmm_hit_features: list[Feature] = field(default_factory=list)
    gff_path: str = ""


def parse_gff_with_fasta(gff_path: str):
    """Parse a GFF file (with embedded FASTA) into features and sequences.

    Returns:
        Tuple of (ir_by_contig, cds_by_contig, hmm_by_contig, sequences)
        where sequences is a dict of contig_id -> sequence_str.
    """
    ir_by_contig = defaultdict(list)
    cds_by_contig = defaultdict(list)
    hmm_by_contig = defaultdict(list)
    fasta_lines = []
    in_fasta = False

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if in_fasta:
                fasta_lines.append(line)
                continue
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feat = Feature(
                seqid=fields[0],
                source=fields[1],
                type=fields[2],
                start=int(fields[3]) - 1,  # convert to 0-based
                end=int(fields[4]),
                strand=fields[6],
                attributes=fields[8],
            )

            if feat.type == "inverted_repeat":
                ir_by_contig[feat.seqid].append(feat)
            elif feat.type == "hmm_hit":
                hmm_by_contig[feat.seqid].append(feat)
            elif feat.type == "CDS":
                cds_by_contig[feat.seqid].append(feat)

    # Parse FASTA
    sequences = {}
    current_id = None
    current_seq = []
    for line in fasta_lines:
        if line.startswith(">"):
            if current_id:
                sequences[current_id] = "".join(current_seq)
            current_id = line[1:].strip().split()[0]
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_id:
        sequences[current_id] = "".join(current_seq)

    return ir_by_contig, cds_by_contig, hmm_by_contig, sequences


def group_features_by_window(features: list[Feature], window_size: int = 3000):
    """Group features that are within window_size bp of each other."""
    if not features:
        return []

    sorted_feats = sorted(features, key=lambda f: f.start)
    groups = [[sorted_feats[0]]]

    for feat in sorted_feats[1:]:
        in_window = any(
            feat.distance_to(g) <= window_size for g in groups[-1]
        )
        if in_window:
            groups[-1].append(feat)
        else:
            groups.append([feat])

    return groups


def find_intersecting_cds(ir_group, cds_features):
    """Find CDS features overlapping the span of an IR group."""
    group_start = min(ir.start for ir in ir_group)
    group_end = max(ir.end for ir in ir_group)
    seqid = ir_group[0].seqid
    window = Feature(seqid, ".", "window", group_start, group_end, "+", "")
    return [cds for cds in cds_features if window.overlaps(cds)]


def _find_invertible_cassettes(
    ir_group: list[Feature],
) -> list[tuple[int, int]]:
    """Identify invertible cassettes from IR pairs in a shufflon cluster.

    Each einverted-detected IR pair defines one invertible cassette: the
    DNA spanning from the left arm through to the right arm (inclusive).
    The cassette boundaries include the repeat arms themselves because
    the coding sequences within shufflon cassettes typically overlap the
    sfx recognition sites.

    Overlapping cassettes (from einverted detecting different pairwise
    combinations of sfx sites) are merged.

    Returns a list of (start, end) tuples in the same coordinate space
    as the input features (0-based start, 1-based end).
    """
    # Group IR features into left/right pairs by their ID prefix.
    # Supports both old (_F/_R) and new (_FOR/_REV) naming conventions.
    pairs = {}
    for ir in ir_group:
        attrs = _parse_gff_attributes(ir.attributes)
        ir_id = attrs.get("ID", "")
        if ir_id.endswith("_FOR"):
            pairs.setdefault(ir_id[:-4], {})["F"] = ir
        elif ir_id.endswith("_REV"):
            pairs.setdefault(ir_id[:-4], {})["R"] = ir
        elif ir_id.endswith("_F"):
            pairs.setdefault(ir_id[:-2], {})["F"] = ir
        elif ir_id.endswith("_R"):
            pairs.setdefault(ir_id[:-2], {})["R"] = ir
        else:
            pairs.setdefault(ir_id, {})["?"] = ir

    # Each complete pair (both F and R) defines a cassette span.
    # Unpaired sites (from motif refinement) are skipped here — they are
    # recognition sites that do not define an invertible segment on their own.
    raw_cassettes = []
    for pair_id, sides in pairs.items():
        f = sides.get("F")
        r = sides.get("R")
        if f and r:
            span_start = min(f.start, r.start)
            span_end = max(f.end, r.end)
            raw_cassettes.append((span_start, span_end))

    if not raw_cassettes:
        return []

    # Sort by start position but do NOT merge overlapping cassettes.
    # Adjacent shufflon cassettes share a boundary repeat (~19 bp overlap),
    # and merging would collapse distinct invertible units into one.
    raw_cassettes.sort()
    return raw_cassettes


def _longest_orf_per_strand(
    seq: str,
    seg_start: int,
    seg_end: int,
    min_aa: int = 20,
) -> list[tuple[int, int, str]]:
    """Find the longest ORF on each strand within a DNA segment.

    Shufflon cassettes can carry one or two coding sequences (one per
    strand).  These are partial genes (truncated C-terminal segments)
    that may lack a canonical start or stop codon at the cassette
    boundaries.  For each strand, this function finds the longest run
    of codons without an internal stop across all three reading frames.

    Returns a list of up to 2 (start, end, strand) tuples in genome
    coordinates.
    """
    segment_seq = seq[seg_start:seg_end]
    if len(segment_seq) < min_aa * 3:
        return []

    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    results = []

    for strand, search_seq in [("+", segment_seq),
                                ("-", segment_seq[::-1].translate(complement))]:
        best = None  # (length_aa, genome_start, genome_end)
        for frame in range(3):
            run_start = frame
            i = frame
            while i + 2 < len(search_seq):
                codon = search_seq[i:i + 3].upper()
                if codon in ("TAA", "TAG", "TGA"):
                    run_len_aa = (i - run_start) // 3
                    if run_len_aa >= min_aa:
                        if strand == "+":
                            gs = seg_start + run_start
                            ge = seg_start + i
                        else:
                            gs = seg_start + len(search_seq) - i
                            ge = seg_start + len(search_seq) - run_start
                        if best is None or run_len_aa > best[0]:
                            best = (run_len_aa, gs, ge)
                    run_start = i + 3
                i += 3

            # Handle run extending to end of segment (partial gene)
            run_len_aa = (i - run_start) // 3
            if run_len_aa >= min_aa:
                if strand == "+":
                    gs = seg_start + run_start
                    ge = seg_start + i
                else:
                    gs = seg_start + len(search_seq) - i
                    ge = seg_start + len(search_seq) - run_start
                if best is None or run_len_aa > best[0]:
                    best = (run_len_aa, gs, ge)

        if best is not None:
            results.append((best[1], best[2], strand))

    return results


def _find_nearby_hmm_hits(
    ir_group: list[Feature],
    hmm_features: list[Feature],
    window_size: int,
) -> list[Feature]:
    """Find hmm_hit features near an IR cluster.

    An HMM hit is included if it overlaps the IR span or is within
    window_size bp of any IR feature in the group.
    """
    group_start = min(ir.start for ir in ir_group)
    group_end = max(ir.end for ir in ir_group)
    seqid = ir_group[0].seqid
    expanded = Feature(
        seqid, ".", "window",
        max(0, group_start - window_size),
        group_end + window_size,
        "+", "",
    )
    return [h for h in hmm_features if expanded.overlaps(h)]


def _group_irs_by_cluster(
    ir_features: list[Feature],
    window_size: int,
) -> list[list[Feature]]:
    """Group IR features by their cluster_id attribute.

    If cluster_id attributes are present, IRs sharing the same cluster_id
    form one group.  This preserves the cluster assignments made during
    shufflon candidate filtering, so that adjacent-but-distinct shufflon
    clusters produce separate windows.

    Falls back to proximity-based grouping (``group_features_by_window``)
    when no cluster_id attributes are found.
    """
    cluster_map: dict[str, list[Feature]] = defaultdict(list)
    unclustered: list[Feature] = []

    for ir in ir_features:
        attrs = _parse_gff_attributes(ir.attributes)
        cid = attrs.get("cluster_id", "")
        if cid:
            cluster_map[cid].append(ir)
        else:
            unclustered.append(ir)

    if not cluster_map:
        # No cluster_id on any feature — fall back to proximity grouping
        return group_features_by_window(ir_features, window_size)

    groups = [sorted(feats, key=lambda f: f.start) for feats in cluster_map.values()]

    # Handle any unclustered IRs via proximity grouping
    if unclustered:
        groups.extend(group_features_by_window(unclustered, window_size))

    return groups


# ── Window classification helpers ────────────────────────────────────────────

_IR_NAME_RE = re.compile(r"inverted_repeat_(\d+)_(FOR|REV)", re.IGNORECASE)


def _ir_pairs_from_features(ir_features: list[Feature]) -> dict[int, dict]:
    """Group IR features into pairs by their pair number.

    Returns a dict mapping pair_number -> {"for": Feature, "rev": Feature}.
    Only pairs with both a FOR and REV arm are included.
    """
    pairs: dict[int, dict] = {}
    for ir in ir_features:
        name = ir.name or ""
        m = _IR_NAME_RE.search(name)
        if not m:
            continue
        pair_num = int(m.group(1))
        direction = m.group(2).upper()
        pairs.setdefault(pair_num, {})
        pairs[pair_num]["for" if direction == "FOR" else "rev"] = ir

    # Keep only complete pairs
    return {k: v for k, v in pairs.items() if "for" in v and "rev" in v}


def _is_shufflon_like(ir_features: list[Feature], max_gap: int = 10) -> bool:
    """Check whether IR features have shufflon-like organisation.

    In a shufflon, sfx recognition sites occur at the boundaries between
    invertible cassettes.  At each boundary, the end of one IR arm and
    the start of the next are immediately adjacent (within ``max_gap``
    bp).  This function scans ALL IR arms in the group (including
    unpaired sites from motif refinement) for these boundary clusters
    and requires at least two distinct boundaries, which is the minimum
    for a three-cassette shufflon.

    Requires at least 3 complete pairs.
    """
    pairs = _ir_pairs_from_features(ir_features)
    if len(pairs) < 3:
        return False

    # Sort all IR arms (including unpaired) by start position.
    # Using every arm rather than just complete-pair spans means
    # unpaired recognition sites discovered by motif search still
    # contribute to boundary detection.
    arms = sorted(ir_features, key=lambda f: f.start)

    # Scan for boundary clusters: positions where the end of one arm
    # and the start of the next are within max_gap.  Adjacent arms
    # that are part of the same boundary cluster are counted once.
    boundaries = 0
    i = 0
    while i < len(arms) - 1:
        gap = arms[i + 1].start - arms[i].end
        if abs(gap) <= max_gap:
            boundaries += 1
            # Skip past all arms belonging to this boundary cluster
            while (
                i < len(arms) - 1
                and abs(arms[i + 1].start - arms[i].end) <= max_gap
            ):
                i += 1
        i += 1

    return boundaries >= 2


def _ir_near_hmm_hit(
    ir_features: list[Feature],
    hmm_features: list[Feature],
    max_distance: int = 2000,
) -> bool:
    """Check if at least one IR arm is within max_distance of an HMM hit."""
    if not hmm_features:
        return False
    for ir in ir_features:
        for hmm in hmm_features:
            if ir.distance_to(hmm) <= max_distance:
                return True
    return False


def extract_shufflon_windows(
    merged_gff: str,
    shufflon_output_dir: str,
    inverton_output_dir: str,
    sample_id: str = "",
    window_size: int = 3000,
    min_ir_pairs: int = 3,
    shufflon_max_gap: int = 10,
    inverton_hmm_distance: int = 2000,
) -> tuple[list[ShufflonWindow], list[ShufflonWindow]]:
    """Extract and classify GFF+FASTA windows from merged annotations.

    Each cluster of nearby IRs is classified as shufflon-like or
    inverton-like:

      * Shufflon-like — at least *min_ir_pairs* complete pairs **and**
        every consecutive pair (sorted by position) is adjacent within
        *shufflon_max_gap* bp (the R64 / TP114 organisation where sfx
        sites sit at cassette boundaries in quick succession).

      * Inverton-like — at least one complete IR pair within
        *inverton_hmm_distance* bp of an HMM-hit CDS.

    Windows matching neither category are skipped.

    Args:
        merged_gff: Merged GFF (Prokka + HMM + IR, with embedded FASTA).
        shufflon_output_dir: Directory for shufflon-like window GFFs.
        inverton_output_dir: Directory for inverton-like window GFFs.
        sample_id: Sample identifier for ShufflonWindow objects.
        window_size: Distance (bp) for clustering IRs / finding nearby HMM hits.
        min_ir_pairs: Minimum IR pairs for shufflon-like classification.
        shufflon_max_gap: Maximum gap (bp) between consecutive pair spans
            for shufflon-like adjacency.
        inverton_hmm_distance: Maximum distance (bp) from any IR arm to an
            HMM hit for inverton-like classification.

    Returns:
        ``(shufflon_windows, inverton_windows)`` — each a list of
        :class:`ShufflonWindow`.
    """
    ensure_dir(shufflon_output_dir)
    ensure_dir(inverton_output_dir)
    input_prefix = os.path.splitext(os.path.basename(merged_gff))[0]

    ir_by_contig, cds_by_contig, hmm_by_contig, sequences = parse_gff_with_fasta(merged_gff)

    if not sequences:
        logger.warning("No FASTA sequences found in %s", merged_gff)
        return [], []

    shufflon_windows: list[ShufflonWindow] = []
    inverton_windows: list[ShufflonWindow] = []
    window_counter = defaultdict(int)

    for contig_id, ir_features in ir_by_contig.items():
        ir_groups = _group_irs_by_cluster(ir_features, window_size)
        contig_cds = cds_by_contig.get(contig_id, [])
        contig_hmm = hmm_by_contig.get(contig_id, [])
        seq = sequences.get(contig_id, "")
        seq_len = len(seq)

        for ir_group in ir_groups:
            n_ir_pairs = len(ir_group) // 2
            if n_ir_pairs < 1:
                continue

            # ── Classify ────────────────────────────────────────────
            nearby_hmm = _find_nearby_hmm_hits(ir_group, contig_hmm, window_size)

            is_shufflon = (
                n_ir_pairs >= min_ir_pairs
                and _is_shufflon_like(ir_group, max_gap=shufflon_max_gap)
            )

            if is_shufflon:
                category = "shufflon_like"
                output_dir = shufflon_output_dir
            elif _ir_near_hmm_hit(ir_group, contig_hmm, max_distance=inverton_hmm_distance):
                category = "inverton_like"
                output_dir = inverton_output_dir
            else:
                logger.debug(
                    "IR group on %s (%d pairs) matches neither shufflon nor inverton criteria, skipping",
                    contig_id, n_ir_pairs,
                )
                continue

            # ── Build window ────────────────────────────────────────
            intersecting_cds = find_intersecting_cds(ir_group, contig_cds)

            if not intersecting_cds and not nearby_hmm:
                continue

            window_counter[contig_id] += 1
            wnum = window_counter[contig_id]

            # Window bounds encompassing all feature types
            all_starts = [ir.start for ir in ir_group]
            all_ends = [ir.end for ir in ir_group]
            for c in intersecting_cds:
                all_starts.append(c.start)
                all_ends.append(c.end)
            for h in nearby_hmm:
                all_starts.append(h.start)
                all_ends.append(h.end)

            win_start = max(0, min(all_starts))
            win_end = min(seq_len, max(all_ends))

            if win_start >= win_end:
                continue

            # Re-collect CDS within the (possibly expanded) window
            win_span = Feature(contig_id, ".", "window", win_start, win_end, "+", "")
            intersecting_cds = [c for c in contig_cds if win_span.overlaps(c)]

            seq_id = f"{input_prefix}_{contig_id}_{wnum}"
            window_id = f"{sample_id}_{contig_id}_w{wnum}" if sample_id else seq_id
            out_path = os.path.join(
                output_dir,
                f"{input_prefix}_contig_{contig_id}_window_{wnum}.gff",
            )

            cassettes = _find_invertible_cassettes(ir_group)

            # ── Write per-window GFF ────────────────────────────────
            with open(out_path, "w") as fh:
                fh.write("##gff-version 3\n")

                for ir in ir_group:
                    ir_source = ir.source if ir.source else "einverted"
                    fh.write(
                        f"{seq_id}\t{ir_source}\tinverted_repeat\t"
                        f"{ir.start - win_start + 1}\t{ir.end - win_start}\t.\t"
                        f"{ir.strand}\t.\t{ir.attributes}\n"
                    )

                for hmm in nearby_hmm:
                    fh.write(
                        f"{seq_id}\tshufflonfinder\thmm_hit\t"
                        f"{hmm.start - win_start + 1}\t{hmm.end - win_start}\t.\t"
                        f"{hmm.strand}\t.\t{hmm.attributes}\n"
                    )

                for cds in intersecting_cds:
                    fh.write(
                        f"{seq_id}\t{cds.source}\tCDS\t"
                        f"{cds.start - win_start + 1}\t{cds.end - win_start}\t.\t"
                        f"{cds.strand}\t.\t{cds.attributes}\n"
                    )

                for si, (seg_s, seg_e) in enumerate(cassettes, start=1):
                    local_s = seg_s - win_start + 1
                    local_e = seg_e - win_start
                    seg_len = seg_e - seg_s
                    fh.write(
                        f"{seq_id}\tshufflonfinder\tinvertible_segment\t"
                        f"{local_s}\t{local_e}\t.\t.\t.\t"
                        f"ID=invertible_segment_{si:03d};"
                        f"Name=Invertible segment {si} ({seg_len} bp)\n"
                    )

                fh.write("##FASTA\n")
                fh.write(f">{seq_id}\n{seq[win_start:win_end]}\n")

            win = ShufflonWindow(
                sample_id=sample_id,
                window_id=window_id,
                contig=contig_id,
                window_start=win_start,
                window_end=win_end,
                n_ir_pairs=n_ir_pairs,
                ir_features=list(ir_group),
                cds_features=list(intersecting_cds),
                hmm_hit_features=list(nearby_hmm),
                gff_path=out_path,
            )

            if category == "shufflon_like":
                shufflon_windows.append(win)
            else:
                inverton_windows.append(win)

    logger.info(
        "Extracted %d shufflon-like + %d inverton-like windows from %s",
        len(shufflon_windows), len(inverton_windows), merged_gff,
    )
    return shufflon_windows, inverton_windows


def shufflon_windows_to_tsv(
    windows: list[ShufflonWindow],
    output_path: str,
) -> pd.DataFrame:
    """Write a summary table of all shufflon windows with IR, CDS, and HMM hit details.

    Each row represents one CDS feature within a window.  Window-level info
    (IR coordinates, HMM hit gene) is repeated on every row.  CDS features
    whose locus_tag matches an HMM hit in the same window are flagged with
    ``is_hmm_hit=True`` and the matching ``hmm_profiles``.

    Args:
        windows: List of ShufflonWindow objects.
        output_path: Where to write the TSV.

    Returns:
        The summary DataFrame.
    """
    rows = []
    for win in windows:
        # Build compact IR coordinate string
        ir_pairs = defaultdict(dict)
        for ir in win.ir_features:
            attrs = _parse_gff_attributes(ir.attributes)
            ir_id = attrs.get("ID", "")
            if ir_id.endswith("_FOR"):
                ir_pairs[ir_id[:-4]]["F"] = ir
            elif ir_id.endswith("_REV"):
                ir_pairs[ir_id[:-4]]["R"] = ir
            elif ir_id.endswith("_F"):
                ir_pairs[ir_id[:-2]]["F"] = ir
            elif ir_id.endswith("_R"):
                ir_pairs[ir_id[:-2]]["R"] = ir
            else:
                ir_pairs[ir_id]["?"] = ir

        ir_coord_parts = []
        for pair_id, sides in ir_pairs.items():
            f = sides.get("F")
            r = sides.get("R")
            if f and r:
                ir_coord_parts.append(
                    f"{f.start + 1}-{f.end}..{r.start + 1}-{r.end}"
                )
            elif f:
                ir_coord_parts.append(f"{f.start + 1}-{f.end}")
            elif r:
                ir_coord_parts.append(f"{r.start + 1}-{r.end}")
        ir_coords_str = "; ".join(ir_coord_parts) if ir_coord_parts else ""

        # Build a lookup: locus_tag -> hmm_profiles for HMM hit genes in
        # this window, so we can flag matching CDS rows
        hmm_lookup = {}
        for hmm in win.hmm_hit_features:
            attrs = hmm.parsed_attributes
            tag = attrs.get("locus_tag", "")
            profiles = attrs.get("hmm_profiles", "")
            if tag:
                hmm_lookup[tag] = profiles

        # Shared window-level columns
        base = {
            "sample_id": win.sample_id,
            "window_id": win.window_id,
            "contig": win.contig,
            "window_start": win.window_start + 1,  # 1-based
            "window_end": win.window_end,
            "window_length_bp": win.window_end - win.window_start,
            "n_ir_pairs": win.n_ir_pairs,
            "ir_coords": ir_coords_str,
        }

        if not win.cds_features:
            row = dict(base)
            row.update({
                "locus_tag": "",
                "cds_start": "",
                "cds_end": "",
                "strand": "",
                "product": "",
                "cds_source": "",
                "is_hmm_hit": False,
                "hmm_profiles": "",
                "gff_path": win.gff_path,
            })
            rows.append(row)
        else:
            for cds in win.cds_features:
                tag = cds.locus_tag
                hit_profiles = hmm_lookup.get(tag, "")
                row = dict(base)
                row.update({
                    "locus_tag": tag,
                    "cds_start": cds.start + 1,
                    "cds_end": cds.end,
                    "strand": cds.strand,
                    "product": cds.product,
                    "cds_source": cds.source,
                    "is_hmm_hit": bool(hit_profiles),
                    "hmm_profiles": hit_profiles,
                    "gff_path": win.gff_path,
                })
                rows.append(row)

    df = pd.DataFrame(rows)
    if not df.empty:
        df.to_csv(output_path, sep="\t", index=False)
    logger.info("Wrote shufflon window summary: %d rows -> %s", len(df), output_path)
    return df
