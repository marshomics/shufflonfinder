"""Step 4: Detect inverted repeats via EMBOSS einverted on flanking regions.

Runs einverted at two sensitivity thresholds (51 and 75), merges the
results with coordinate-based deduplication (preserving distinct pairs in
dense clusters), extracts arm sequences, computes percent identity, and
applies shufflon-specific filters: minimum cluster density and at least
one IR arm overlapping a CDS.
"""

import logging
import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import pandas as pd

from .sample_sheet import Sample
from .step_flanking import FlankingRegion, parse_fasta_file
from .utils import ensure_dir, run_cmd

logger = logging.getLogger("shufflonfinder")

# Tolerance (bp) for considering two IR pairs as the same detection
_COORD_DEDUP_TOLERANCE = 3


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class InvertedRepeat:
    """One pair of inverted repeat arms detected by einverted."""
    chrom: str
    left_start: int   # 0-based
    left_stop: int     # 1-based
    right_start: int   # 0-based
    right_stop: int    # 1-based
    left_seq: str = ""
    right_seq: str = ""
    middle_seq: str = ""
    percent_identity: float = 0.0


# ---------------------------------------------------------------------------
# Running einverted
# ---------------------------------------------------------------------------

def _run_einverted(
    fasta_path: str,
    out_dir: str,
    threshold: int,
    match: int,
    mismatch: int,
) -> tuple[str, str]:
    """Run EMBOSS einverted with given parameters.

    Returns (outfile_path, outseq_path).
    """
    ensure_dir(out_dir)
    outfile = os.path.join(out_dir, f"einverted.{threshold}.outfile")
    outseq = os.path.join(out_dir, f"einverted.{threshold}.outseq")

    cmd = [
        "einverted",
        "-maxrepeat", "750",
        "-gap", "100",
        "-threshold", str(threshold),
        "-match", str(match),
        "-mismatch", str(mismatch),
        "-outfile", outfile,
        "-outseq", outseq,
        "-sequence", fasta_path,
    ]
    run_cmd(cmd, description=f"einverted (threshold={threshold})")
    return outfile, outseq


def run_einverted_dual(
    fasta_path: str,
    out_dir: str,
    cpus: int = 4,
) -> tuple[str, str]:
    """Run einverted at two thresholds (51 and 75) in parallel.

    Returns (outfile_51, outfile_75).
    """
    intermediate = ensure_dir(os.path.join(out_dir, "intermediate"))

    def _run51():
        return _run_einverted(fasta_path, intermediate, 51, match=5, mismatch=-9)

    def _run75():
        return _run_einverted(fasta_path, intermediate, 75, match=5, mismatch=-15)

    with ThreadPoolExecutor(max_workers=min(2, cpus)) as pool:
        f51 = pool.submit(_run51)
        f75 = pool.submit(_run75)
        outfile_51 = f51.result()[0]
        outfile_75 = f75.result()[0]

    return outfile_51, outfile_75


# ---------------------------------------------------------------------------
# Parsing einverted output
# ---------------------------------------------------------------------------

def _parse_einverted_outfile(path: str) -> list[tuple[str, int, int, int, int]]:
    """Parse an einverted .outfile into coordinate tuples.

    einverted output is grouped in 5-line blocks:
        line 0: blank / separator
        line 1: "<contig>: <info>"
        line 2: "<left_start>  <seq>  <left_stop>"
        line 3: alignment characters
        line 4: "<right_stop>  <seq>  <right_start>"

    Returns list of (chrom, left_start, left_stop, right_start, right_stop)
    with 0-based start, 1-based stop.
    """
    if not os.path.isfile(path):
        return []

    irs = []
    with open(path) as fh:
        lines = []
        for line in fh:
            lines.append(line)
            if len(lines) >= 5:
                chrom = lines[1].split(":")[0].strip()
                left_parts = lines[2].split()
                right_parts = lines[4].split()

                left_start = int(left_parts[0]) - 1   # to 0-based
                left_stop = int(left_parts[-1])

                right_start = int(right_parts[-1]) - 1  # to 0-based
                right_stop = int(right_parts[0])

                irs.append((chrom, left_start, left_stop, right_start, right_stop))
                lines = []

    return irs


def _is_coordinate_duplicate(
    a: tuple[str, int, int, int, int],
    b: tuple[str, int, int, int, int],
    tolerance: int = _COORD_DEDUP_TOLERANCE,
) -> bool:
    """Check if two IR tuples represent the same detection.

    Two IRs are duplicates if they're on the same contig and all four
    coordinate values are within `tolerance` bp of each other.  This is
    much more conservative than the old overlap-based check, which
    incorrectly discarded distinct pairs in dense clusters.
    """
    if a[0] != b[0]:
        return False
    return (
        abs(a[1] - b[1]) <= tolerance
        and abs(a[2] - b[2]) <= tolerance
        and abs(a[3] - b[3]) <= tolerance
        and abs(a[4] - b[4]) <= tolerance
    )


def merge_einverted_results(
    outfile_51: str,
    outfile_75: str,
) -> list[tuple[str, int, int, int, int]]:
    """Merge results from two einverted runs with coordinate-based dedup.

    Pools all IRs from both threshold runs.  When two IRs from different
    runs have nearly identical coordinates (both arms within 3 bp), the
    duplicate is removed.  Unlike the old overlap-based merge, this
    preserves distinct IR pairs that happen to share a dense genomic
    region (exactly the pattern seen in shufflons).
    """
    irs_75 = _parse_einverted_outfile(outfile_75)
    irs_51 = _parse_einverted_outfile(outfile_51)

    # Start with all threshold-75 results
    merged = list(irs_75)

    # Add threshold-51 results that aren't coordinate-duplicates of a t75 hit
    for ir51 in irs_51:
        if not any(_is_coordinate_duplicate(ir51, existing) for existing in merged):
            merged.append(ir51)

    return merged


# ---------------------------------------------------------------------------
# Sequence extraction and identity computation
# ---------------------------------------------------------------------------

def _compute_percent_identity(seq_a: str, seq_b: str) -> float:
    """Compute percent identity between two IR arm sequences.

    Compares seq_a to the reverse complement of seq_b (since IR arms are
    complementary-reversed).
    """
    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    rev_comp_b = seq_b[::-1].translate(complement)

    compare_len = min(len(seq_a), len(rev_comp_b))
    if compare_len == 0:
        return 0.0

    matches = sum(
        1 for i in range(compare_len)
        if seq_a[i].upper() == rev_comp_b[i].upper()
    )
    return (matches / compare_len) * 100.0


def annotate_ir_sequences(
    raw_irs: list[tuple[str, int, int, int, int]],
    sequences: dict[str, str],
    min_middle_bp: int = 30,
) -> list[InvertedRepeat]:
    """Extract arm/middle sequences for each IR and compute percent identity.

    Also applies the minimum inter-arm distance filter (default 30 bp).

    Args:
        raw_irs: Coordinate tuples from merged einverted output.
        sequences: Dict of chrom -> DNA sequence.
        min_middle_bp: Minimum length of the sequence between the two arms.

    Returns:
        List of fully annotated InvertedRepeat objects.
    """
    results = []
    for chrom, ls, le, rs, re_ in raw_irs:
        seq = sequences.get(chrom)
        if seq is None:
            logger.debug("Contig '%s' not in FASTA, skipping IR", chrom)
            continue

        left_seq = seq[ls:le]
        right_seq = seq[rs:re_]
        middle_seq = seq[le:rs]

        if len(middle_seq) < min_middle_bp:
            continue

        pct_id = _compute_percent_identity(left_seq, right_seq)

        results.append(InvertedRepeat(
            chrom=chrom,
            left_start=ls,
            left_stop=le,
            right_start=rs,
            right_stop=re_,
            left_seq=left_seq,
            right_seq=right_seq,
            middle_seq=middle_seq,
            percent_identity=pct_id,
        ))

    return results


def irs_to_dataframe(irs: list[InvertedRepeat]) -> pd.DataFrame:
    """Convert InvertedRepeat objects to a DataFrame."""
    rows = []
    for ir in irs:
        rows.append({
            "IR_Chr": ir.chrom,
            "LeftIRStart": ir.left_start,
            "LeftIRStop": ir.left_stop,
            "RightIRStart": ir.right_start,
            "RightIRStop": ir.right_stop,
            "LeftIRSequence": ir.left_seq,
            "InvertibleSequence": ir.middle_seq,
            "RightIRSequence": ir.right_seq,
            "PercentIdentity": round(ir.percent_identity, 2),
        })
    return pd.DataFrame(rows)


def export_irs_tsv(ir_df: pd.DataFrame, out_dir: str) -> str:
    """Write IRs.tsv to the data/ subdirectory."""
    data_dir = ensure_dir(os.path.join(out_dir, "data"))
    path = os.path.join(data_dir, "IRs.tsv")
    ir_df.to_csv(path, sep="\t", index=False)
    logger.info("Wrote %d IRs to %s", len(ir_df), path)
    return path


# ---------------------------------------------------------------------------
# Top-level detection entry point
# ---------------------------------------------------------------------------

def detect_inverted_repeats(
    sample: Sample,
    flanking_fasta: str,
    outdir: str,
    cpus: int = 4,
) -> str:
    """Detect inverted repeats in flanking DNA via einverted.

    Runs einverted at two thresholds, merges with coordinate-based
    deduplication, extracts sequences, computes percent identity, and
    writes IRs.tsv.

    Args:
        sample: The Sample being processed.
        flanking_fasta: Path to the multi-record FASTA of flanking regions.
        outdir: Base output directory.
        cpus: Number of threads.

    Returns:
        Path to the output directory for this sample (contains data/IRs.tsv).
    """
    if not flanking_fasta or not os.path.isfile(flanking_fasta):
        logger.warning("No flanking FASTA for %s, skipping IR detection", sample.sample_id)
        return ""

    with open(flanking_fasta) as fh:
        if not fh.read(1):
            logger.warning("Empty flanking FASTA for %s, skipping IR detection", sample.sample_id)
            return ""

    ensure_dir(outdir)

    outfile_51, outfile_75 = run_einverted_dual(flanking_fasta, outdir, cpus=cpus)

    raw_irs = merge_einverted_results(outfile_51, outfile_75)
    logger.info(
        "einverted found %d IRs for %s (after coordinate dedup)",
        len(raw_irs), sample.sample_id,
    )

    sequences = parse_fasta_file(flanking_fasta)
    annotated = annotate_ir_sequences(raw_irs, sequences)
    logger.info(
        "%d IRs passed min-middle-distance filter for %s",
        len(annotated), sample.sample_id,
    )

    ir_df = irs_to_dataframe(annotated)
    export_irs_tsv(ir_df, outdir)

    logger.info("IR detection complete for %s", sample.sample_id)
    return outdir


# ---------------------------------------------------------------------------
# Loading, remapping, combining
# ---------------------------------------------------------------------------

def load_ir_table(ir_dir: str) -> pd.DataFrame:
    """Load the inverted repeats table from an output directory."""
    if not ir_dir:
        return pd.DataFrame()

    ir_path = os.path.join(ir_dir, "data", "IRs.tsv")
    if not os.path.isfile(ir_path):
        logger.warning("No IRs.tsv found at %s", ir_path)
        return pd.DataFrame()

    df = pd.read_csv(ir_path, sep="\t", dtype={"IR_Chr": str})
    for col in ("LeftIRStart", "LeftIRStop", "RightIRStart", "RightIRStop"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    logger.info("Loaded %d inverted repeats from %s", len(df), ir_path)
    return df


def remap_ir_to_genome_coords(
    ir_df: pd.DataFrame,
    flanking_regions: list[FlankingRegion],
) -> pd.DataFrame:
    """Convert IR coordinates from flanking-region-local to genome-absolute."""
    if ir_df.empty:
        return ir_df

    region_map = {r.hit_id: r for r in flanking_regions}

    rows = []
    for _, ir_row in ir_df.iterrows():
        ir_chr = str(ir_row["IR_Chr"])
        region = region_map.get(ir_chr)
        if region is None:
            logger.debug("IR_Chr '%s' doesn't match any flanking region, skipping", ir_chr)
            continue

        offset = region.flank_start - 1

        remapped = ir_row.copy()
        remapped["IR_Chr"] = region.contig
        remapped["LeftIRStart"] = ir_row["LeftIRStart"] + offset
        remapped["LeftIRStop"] = ir_row["LeftIRStop"] + offset
        remapped["RightIRStart"] = ir_row["RightIRStart"] + offset
        remapped["RightIRStop"] = ir_row["RightIRStop"] + offset

        remapped["sample_id"] = region.sample_id
        remapped["hmm_profiles"] = ";".join(region.hmm_profiles)
        remapped["locus_tag"] = region.locus_tag
        remapped["hit_id"] = region.hit_id

        rows.append(remapped)

    if not rows:
        logger.warning("No IRs could be remapped to genome coordinates")
        return pd.DataFrame()

    result = pd.DataFrame(rows)
    logger.info("Remapped %d IRs to genome coordinates", len(result))
    return result


def filter_ir_table(
    ir_df: pd.DataFrame,
    min_arm_length: int = 0,
    max_arm_length: int = 0,
    min_identity: float = 0.0,
) -> pd.DataFrame:
    """Filter inverted repeats by arm length and percent identity."""
    if ir_df.empty:
        return ir_df

    before = len(ir_df)

    if min_arm_length > 0:
        left_len = (ir_df["LeftIRStop"] - ir_df["LeftIRStart"]).abs() + 1
        right_len = (ir_df["RightIRStop"] - ir_df["RightIRStart"]).abs() + 1
        mask = (left_len >= min_arm_length) & (right_len >= min_arm_length)
        ir_df = ir_df.loc[mask].copy()
        logger.info(
            "Arm-length filter (>= %d bp): %d -> %d IRs",
            min_arm_length, before, len(ir_df),
        )

    if max_arm_length > 0 and not ir_df.empty:
        before_max = len(ir_df)
        left_len = (ir_df["LeftIRStop"] - ir_df["LeftIRStart"]).abs() + 1
        right_len = (ir_df["RightIRStop"] - ir_df["RightIRStart"]).abs() + 1
        mask = (left_len <= max_arm_length) & (right_len <= max_arm_length)
        ir_df = ir_df.loc[mask].copy()
        logger.info(
            "Arm-length filter (<= %d bp): %d -> %d IRs",
            max_arm_length, before_max, len(ir_df),
        )

    if min_identity > 0 and not ir_df.empty:
        identity_col = None
        for candidate in ("PercentIdentity", "Percent", "Match", "Score", "Identity"):
            if candidate in ir_df.columns:
                identity_col = candidate
                break

        if identity_col is None:
            logger.warning(
                "No percent-identity column found in IR table. "
                "Skipping identity filter."
            )
        else:
            before_id = len(ir_df)
            vals = pd.to_numeric(ir_df[identity_col], errors="coerce")
            ir_df = ir_df.loc[vals >= min_identity].copy()
            logger.info(
                "Identity filter (>= %.1f%%, column '%s'): %d -> %d IRs",
                min_identity, identity_col, before_id, len(ir_df),
            )

    removed = before - len(ir_df)
    if removed:
        logger.info("IR filtering removed %d of %d records", removed, before)

    return ir_df


def filter_irs_in_hmm_hits(
    ir_df: pd.DataFrame,
    flanking_regions: list[FlankingRegion],
) -> pd.DataFrame:
    """Remove IRs where either arm falls inside an HMM-hit CDS.

    sfx recognition sites sit between invertible gene cassettes, not
    inside the recombinase gene itself.  einverted sometimes reports
    palindromic structure within the recombinase coding sequence as
    an inverted repeat; these are false positives for shufflon
    annotation purposes.

    An IR is removed if *either* arm (left or right) is fully
    contained within the CDS coordinates of any HMM hit on the same
    sample and contig.

    Args:
        ir_df: IR DataFrame with genome-absolute coordinates (must
            contain ``sample_id``, ``IR_Chr``, ``LeftIRStart``,
            ``LeftIRStop``, ``RightIRStart``, ``RightIRStop``).
        flanking_regions: List of FlankingRegion objects from step 3,
            each carrying ``sample_id``, ``contig``, ``cds_start``,
            ``cds_end`` for the HMM-hit CDS.

    Returns:
        Filtered DataFrame with offending IRs removed.
    """
    if ir_df.empty or not flanking_regions:
        return ir_df

    # Build lookup: (sample_id, contig) -> list of (cds_start, cds_end)
    hmm_cds: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for region in flanking_regions:
        key = (region.sample_id, region.contig)
        hmm_cds.setdefault(key, []).append((region.cds_start, region.cds_end))

    def _arm_inside_any_cds(row) -> bool:
        key = (row["sample_id"], row["IR_Chr"])
        cds_list = hmm_cds.get(key)
        if not cds_list:
            return False
        ls, le = int(row["LeftIRStart"]), int(row["LeftIRStop"])
        rs, re = int(row["RightIRStart"]), int(row["RightIRStop"])
        for cds_s, cds_e in cds_list:
            left_inside = (ls >= cds_s and le <= cds_e)
            right_inside = (rs >= cds_s and re <= cds_e)
            if left_inside or right_inside:
                return True
        return False

    mask = ir_df.apply(_arm_inside_any_cds, axis=1)
    n_removed = mask.sum()
    if n_removed:
        logger.info(
            "Removed %d IR(s) with arm(s) inside HMM-hit CDS regions",
            n_removed,
        )
    return ir_df.loc[~mask].reset_index(drop=True)


def _deduplicate_ir_by_coords(ir_df: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate IRs that share identical genome-absolute coordinates.

    The same physical IR pair can appear multiple times when it falls
    within the flanking regions of several HMM hits.  After remapping
    to genome-absolute coordinates these rows are identical except for
    the ``locus_tag`` / ``hit_id`` / ``hmm_profiles`` columns.  We
    keep one representative row per unique coordinate set, merging the
    HMM profile and locus tag information so nothing is lost.
    """
    if ir_df.empty:
        return ir_df

    coord_cols = [
        "sample_id", "IR_Chr",
        "LeftIRStart", "LeftIRStop",
        "RightIRStart", "RightIRStop",
    ]
    # Verify all coordinate columns are present
    missing = [c for c in coord_cols if c not in ir_df.columns]
    if missing:
        logger.warning(
            "Cannot deduplicate IRs: missing columns %s", missing,
        )
        return ir_df

    before = len(ir_df)

    def _merge_group(group: pd.DataFrame) -> pd.Series:
        """Collapse duplicate rows into one, merging metadata."""
        row = group.iloc[0].copy()
        if "hmm_profiles" in group.columns:
            all_profiles = set()
            for val in group["hmm_profiles"].dropna():
                for p in str(val).split(";"):
                    p = p.strip()
                    if p:
                        all_profiles.add(p)
            row["hmm_profiles"] = ";".join(sorted(all_profiles))
        if "locus_tag" in group.columns:
            tags = sorted(set(group["locus_tag"].dropna().astype(str)))
            row["locus_tag"] = ";".join(tags)
        if "hit_id" in group.columns:
            ids = sorted(set(group["hit_id"].dropna().astype(str)))
            row["hit_id"] = ";".join(ids)
        return row

    deduped = (
        ir_df.groupby(coord_cols, sort=False)
        .apply(_merge_group, include_groups=False)
        .reset_index(drop=True)
    )

    removed = before - len(deduped)
    if removed:
        logger.info(
            "Coordinate deduplication: %d -> %d IRs "
            "(%d duplicates from overlapping flanking regions removed)",
            before, len(deduped), removed,
        )
    return deduped


def combine_ir_tables(
    ir_results: list[tuple[str, str, list[FlankingRegion]]],
    output_path: str,
) -> pd.DataFrame:
    """Combine and remap IR tables from multiple einverted runs."""
    frames = []
    for ir_dir, sid, flanking_regions in ir_results:
        raw_ir = load_ir_table(ir_dir)
        if raw_ir.empty:
            continue
        remapped = remap_ir_to_genome_coords(raw_ir, flanking_regions)
        if not remapped.empty:
            frames.append(remapped)

    if not frames:
        logger.warning("No inverted repeats found across any samples.")
        combined = pd.DataFrame()
    else:
        combined = pd.concat(frames, ignore_index=True)
        combined = _deduplicate_ir_by_coords(combined)

    combined.to_csv(output_path, sep="\t", index=False)
    logger.info("Combined IR table: %d records -> %s", len(combined), output_path)
    return combined


# ---------------------------------------------------------------------------
# Shufflon candidate filtering: density + CDS overlap
# ---------------------------------------------------------------------------

def _ir_span(row: pd.Series) -> tuple[int, int]:
    """Return the full span (min_coord, max_coord) of an IR pair."""
    coords = [row["LeftIRStart"], row["LeftIRStop"],
              row["RightIRStart"], row["RightIRStop"]]
    return int(min(coords)), int(max(coords))


def _cluster_ir_rows(
    ir_df: pd.DataFrame,
    window_size: int,
) -> list[list[int]]:
    """Cluster IR rows by proximity on the same contig.

    Returns a list of clusters, where each cluster is a list of row indices.
    Two IR pairs are in the same cluster if any of their arms is within
    window_size bp of any arm in the cluster.
    """
    if ir_df.empty:
        return []

    # Sort by contig then by leftmost coordinate
    ir_df = ir_df.copy()
    ir_df["_span_min"] = ir_df.apply(lambda r: _ir_span(r)[0], axis=1)
    sorted_idx = ir_df.sort_values(["IR_Chr", "_span_min"]).index.tolist()

    clusters = []
    current_cluster = [sorted_idx[0]]
    current_contig = ir_df.loc[sorted_idx[0], "IR_Chr"]
    current_max = _ir_span(ir_df.loc[sorted_idx[0]])[1]

    for idx in sorted_idx[1:]:
        row = ir_df.loc[idx]
        span_min, span_max = _ir_span(row)
        contig = row["IR_Chr"]

        if contig == current_contig and span_min <= current_max + window_size:
            current_cluster.append(idx)
            current_max = max(current_max, span_max)
        else:
            clusters.append(current_cluster)
            current_cluster = [idx]
            current_contig = contig
            current_max = span_max

    clusters.append(current_cluster)
    return clusters



def _cluster_density(ir_df: pd.DataFrame, cluster_idx: list[int]) -> float:
    """Compute IR pair density (pairs per kilobase) for a cluster."""
    if len(cluster_idx) < 2:
        return float("inf")  # single-pair clusters are infinitely dense

    spans = [_ir_span(ir_df.loc[idx]) for idx in cluster_idx]
    cluster_start = min(s[0] for s in spans)
    cluster_end = max(s[1] for s in spans)
    span_kb = max((cluster_end - cluster_start) / 1000.0, 0.001)
    return len(cluster_idx) / span_kb


def _remove_nested_pairs(
    ir_df: pd.DataFrame,
    cluster_idx: list[int],
) -> list[int]:
    """Remove IR pairs whose span is entirely contained within another pair.

    An IR pair's span runs from its leftmost arm coordinate to its
    rightmost arm coordinate.  If pair B's span falls entirely within
    pair A's span, B is considered nested and is removed.  When two
    pairs have identical spans, neither is treated as nested (they are
    kept as potential distinct sfx combinations at the same locus).

    Args:
        ir_df: Full IR DataFrame (indexed so .loc works).
        cluster_idx: Row indices belonging to one cluster.

    Returns:
        Filtered list of row indices with nested pairs removed.
    """
    if len(cluster_idx) <= 1:
        return cluster_idx

    spans = {idx: _ir_span(ir_df.loc[idx]) for idx in cluster_idx}
    nested = set()

    for idx_b in cluster_idx:
        b_min, b_max = spans[idx_b]
        for idx_a in cluster_idx:
            if idx_a == idx_b:
                continue
            a_min, a_max = spans[idx_a]
            # B is strictly nested inside A (not identical spans)
            if a_min <= b_min and a_max >= b_max and (a_min, a_max) != (b_min, b_max):
                nested.add(idx_b)
                break

    if nested:
        logger.debug(
            "Removed %d nested IR pair(s) from cluster of %d",
            len(nested), len(cluster_idx),
        )
    return [idx for idx in cluster_idx if idx not in nested]


def filter_shufflon_candidates(
    ir_df: pd.DataFrame,
    samples: list[Sample],
    cluster_distance: int = 1000,
    min_ir_pairs: int = 3,
    min_ir_density: float = 1.0,
    window_size: int = 3000,
) -> pd.DataFrame:
    """Filter IRs to keep only dense, shufflon-like clusters.

    For each cluster the pipeline:
      1. Removes nested IR pairs (pairs whose span falls entirely
         within another pair's span).
      2. Requires at least ``min_ir_pairs`` remaining pairs.
      3. Requires a density of at least ``min_ir_density`` pairs/kb.

    Every IR in the input already comes from the flanking region of a
    recombinase HMM hit (within ±flank_bp of the gene), so proximity to a
    recombinase is guaranteed by construction.  The count and density
    requirements are sufficient to distinguish real shufflon recognition
    site clusters from sparse transposon or IS-element IRs.

    Note: CDS overlap is intentionally NOT required.  Shufflon cassettes
    contain partial genes (truncated C-terminal segments) that Prodigal
    cannot predict, so a CDS overlap test would reject genuine shufflons.
    The pipeline's own ORF finder (step_gff) discovers these partial
    CDS after the cluster has been accepted.

    Args:
        ir_df: Remapped IR DataFrame (genome-absolute coordinates).
        samples: List of all Sample objects.
        cluster_distance: Maximum gap (bp) between IR pairs for chaining
            into one cluster.  Default 1000 (tighter than the window
            extraction distance, to isolate the dense shufflon core).
        min_ir_pairs: Minimum number of IR pairs per cluster.
        min_ir_density: Minimum IR pairs per kilobase of cluster span.
        window_size: Passed through for compatibility but NOT used for
            clustering; use ``cluster_distance`` instead.

    Returns:
        Filtered DataFrame containing only IRs from qualifying clusters,
        plus a new ``cluster_id`` column.
    """
    if ir_df.empty:
        return ir_df

    before = len(ir_df)
    n_nested_total = 0
    keep_indices = []
    cluster_labels = {}

    for sample_id, sample_irs in ir_df.groupby("sample_id"):
        clusters = _cluster_ir_rows(sample_irs, cluster_distance)

        for ci, cluster_idx in enumerate(clusters):
            # --- remove nested pairs first ---
            before_nesting = len(cluster_idx)
            cluster_idx = _remove_nested_pairs(ir_df, cluster_idx)
            n_nested_total += before_nesting - len(cluster_idx)

            # --- count check (applied after nesting removal) ---
            if len(cluster_idx) < min_ir_pairs:
                continue

            # --- density check ---
            density = _cluster_density(ir_df, cluster_idx)
            if density < min_ir_density:
                spans = [_ir_span(ir_df.loc[idx]) for idx in cluster_idx]
                span_bp = max(s[1] for s in spans) - min(s[0] for s in spans)
                logger.debug(
                    "Cluster with %d pairs over %d bp (%.2f pairs/kb) "
                    "rejected by density filter (min %.2f)",
                    len(cluster_idx), span_bp, density, min_ir_density,
                )
                continue

            cid = f"{sample_id}_cluster_{ci + 1}"
            for idx in cluster_idx:
                keep_indices.append(idx)
                cluster_labels[idx] = cid

    if n_nested_total:
        logger.info("Removed %d nested IR pair(s) across all clusters", n_nested_total)

    if not keep_indices:
        logger.warning("No IR clusters met the shufflon candidate criteria")
        return pd.DataFrame(columns=list(ir_df.columns) + ["cluster_id"])

    result = ir_df.loc[keep_indices].copy()
    result["cluster_id"] = [cluster_labels[i] for i in keep_indices]

    removed = before - len(result)
    n_clusters = result["cluster_id"].nunique()
    logger.info(
        "Shufflon candidate filter: %d -> %d IRs in %d cluster(s) "
        "(removed %d IRs: %d nested + %d from sparse/small regions)",
        before, len(result), n_clusters, removed,
        n_nested_total, removed - n_nested_total,
    )
    return result


# ---------------------------------------------------------------------------
# Motif-based refinement: detect sfx recognition sites missed by einverted
# ---------------------------------------------------------------------------

def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq[::-1].translate(complement)


def _derive_core_motif(arm_sequences: list[str], min_core: int = 8) -> str:
    """Derive the longest shared core substring from a set of sfx arm sequences.

    All sfx-type recognition sites in a shufflon share a conserved core
    (e.g., ``GCCAATCCGG`` in R64).  This function finds that core by
    identifying the longest substring common to all input arm sequences.

    Only forward-orientation arms should be passed in (the function
    reverse-complements R-arms internally if needed).

    Args:
        arm_sequences: List of arm DNA sequences (all same orientation).
        min_core: Minimum core length to accept.

    Returns:
        The longest common substring, or "" if none ≥ min_core.
    """
    if not arm_sequences:
        return ""

    # Use the shortest sequence as reference
    ref = min(arm_sequences, key=len)
    best = ""

    for length in range(len(ref), min_core - 1, -1):
        for start in range(len(ref) - length + 1):
            candidate = ref[start:start + length].upper()
            if all(candidate in s.upper() for s in arm_sequences):
                return candidate

    return best


def refine_sfx_sites(
    ir_df: pd.DataFrame,
    sequences: dict[str, str],
    search_margin: int = 500,
) -> pd.DataFrame:
    """Detect sfx recognition sites missed by einverted within shufflon clusters.

    einverted finds IR pairs — sequences that are reverse complements of
    each other.  In shufflons, the recognition sites (sfx sites) all share
    a conserved core motif, but cross-type pairings may have too many
    mismatches for einverted to detect.  This function:

    1. Derives a consensus core motif from the arms already detected.
    2. Searches the cluster region for all occurrences of the core
       (forward and reverse strand).
    3. Adds any site not already covered by a detected arm as a new
       IR row, paired with its best-matching existing arm.

    Args:
        ir_df: Shufflon candidate IR DataFrame (with ``cluster_id``).
        sequences: Dict of contig_id -> full genome sequence.
        search_margin: Extra bp to search beyond the cluster span on
            each side.

    Returns:
        Augmented DataFrame with additional IR rows for newly detected
        sites.
    """
    if ir_df.empty or "cluster_id" not in ir_df.columns:
        return ir_df

    new_rows = []

    for cluster_id, cluster in ir_df.groupby("cluster_id"):
        contig = cluster["IR_Chr"].iloc[0]
        seq = sequences.get(contig, "")
        if not seq:
            continue

        # Collect all arm sequences, normalised to forward orientation
        # Forward arms: LeftIRSequence; Reverse arms: RC of RightIRSequence
        fwd_arms = []
        for _, row in cluster.iterrows():
            fwd_arms.append(str(row["LeftIRSequence"]).upper())
            fwd_arms.append(_reverse_complement(str(row["RightIRSequence"])).upper())

        core = _derive_core_motif(fwd_arms)
        if not core or len(core) < 8:
            logger.debug(
                "Cluster %s: could not derive core motif (best=%r)",
                cluster_id, core,
            )
            continue

        core_rc = _reverse_complement(core)

        # Determine search region
        all_coords = []
        for _, row in cluster.iterrows():
            all_coords.extend([
                int(row["LeftIRStart"]), int(row["LeftIRStop"]),
                int(row["RightIRStart"]), int(row["RightIRStop"]),
            ])
        region_start = max(0, min(all_coords) - search_margin)
        region_end = min(len(seq), max(all_coords) + search_margin)
        region_seq = seq[region_start:region_end].upper()

        # Identify "extended" IR pairs whose arms bridge inner boundaries
        # (spanning from one sfx site into an adjacent one).  These are
        # einverted artifacts from matching long stretches that include
        # parts of two neighbouring sites.  Detection: any pair whose
        # BOTH arms exceed 1.5× the median arm length is likely extended.
        # We require 3+ pairs to have a meaningful median.
        import statistics
        all_arm_lens = []
        for _, row in cluster.iterrows():
            all_arm_lens.append(int(row["LeftIRStop"]) - int(row["LeftIRStart"]))
            all_arm_lens.append(int(row["RightIRStop"]) - int(row["RightIRStart"]))
        median_arm = statistics.median(all_arm_lens) if all_arm_lens else 0

        extended_indices = set()
        if len(cluster) >= 3 and median_arm > 0:
            for idx, row in cluster.iterrows():
                l_len = int(row["LeftIRStop"]) - int(row["LeftIRStart"])
                r_len = int(row["RightIRStop"]) - int(row["RightIRStart"])
                if l_len > median_arm * 1.4 and r_len > median_arm * 1.4:
                    extended_indices.add(idx)

        if extended_indices:
            logger.info(
                "Cluster %s: removing %d extended IR pair(s) that bridge "
                "inner boundaries (arm ratio > 1.4× median %dbp)",
                cluster_id, len(extended_indices), median_arm,
            )

        # Build existing_spans from non-extended pairs only.
        # This allows the motif search to find individual sites within
        # regions previously covered by extended arms.
        existing_spans = set()
        for idx, row in cluster.iterrows():
            if idx in extended_indices:
                continue
            existing_spans.add((int(row["LeftIRStart"]), int(row["LeftIRStop"])))
            existing_spans.add((int(row["RightIRStart"]), int(row["RightIRStop"])))

        def _overlaps_existing(start: int, end: int, tolerance: int = 5) -> bool:
            for es, ee in existing_spans:
                if abs(start - es) <= tolerance and abs(end - ee) <= tolerance:
                    return True
                # Check if the new site is contained within an existing arm
                if start >= es - tolerance and end <= ee + tolerance:
                    return True
                # Check reciprocal overlap: if >50% of either span overlaps
                overlap_start = max(start, es)
                overlap_end = min(end, ee)
                if overlap_end > overlap_start:
                    overlap_len = overlap_end - overlap_start
                    new_len = max(1, end - start)
                    exist_len = max(1, ee - es)
                    if overlap_len / new_len > 0.5 or overlap_len / exist_len > 0.5:
                        return True
            return False

        # Compute typical individual sfx site length from detected arms.
        # einverted arms may span two abutting sites, so use the MINIMUM
        # arm length from NON-EXTENDED arms as the best estimate for
        # single-site width.
        arm_lengths = []
        for idx, row in cluster.iterrows():
            if idx in extended_indices:
                continue
            arm_lengths.append(int(row["LeftIRStop"]) - int(row["LeftIRStart"]))
            arm_lengths.append(int(row["RightIRStop"]) - int(row["RightIRStart"]))
        if not arm_lengths:
            # All pairs are extended; fall back to core length + 6
            # (typical sfx site is core + 6bp extension)
            arm_lengths = [len(core) + 6]
        site_len = min(arm_lengths)

        # Collect existing R-arm RCs and F-arm sequences for partner matching.
        # This lets us determine the correct offset of the core within a site
        # by trying multiple extensions and picking the one with best identity.
        r_arm_rcs = []  # (rc_seq, row) for each existing pair
        f_arm_seqs = []
        for _, row in cluster.iterrows():
            r_arm_rcs.append((_reverse_complement(str(row["RightIRSequence"])), row))
            f_arm_seqs.append((str(row["LeftIRSequence"]), row))

        def _best_site_boundaries(
            core_genome_pos: int, strand: str,
        ) -> tuple[int, int, str, float] | None:
            """Try multiple offsets of the core within a site-length window.

            Returns (site_start, site_end, site_seq, best_identity) or None.
            """
            best = None
            # The core can appear at any offset 0..site_len-len(core)
            for offset in range(max(1, site_len - len(core) + 1)):
                s = core_genome_pos - offset
                e = s + site_len
                s = max(0, s)
                e = min(len(seq), e)
                if _overlaps_existing(s, e):
                    continue
                candidate = seq[s:e]

                # Compare against partner arms
                if strand == "+":
                    targets = [rc for rc, _ in r_arm_rcs]
                else:
                    targets = [_reverse_complement(fseq) for fseq, _ in f_arm_seqs]

                for target in targets:
                    minlen = min(len(candidate), len(target))
                    if minlen == 0:
                        continue
                    matches = sum(
                        1 for i in range(minlen)
                        if candidate[i].upper() == target[i].upper()
                    )
                    identity = (matches / minlen) * 100.0
                    if best is None or identity > best[3]:
                        best = (s, e, candidate, identity)

            return best

        # Search for core motif occurrences (forward strand)
        found_sites = []  # (genome_start, genome_end, strand, site_seq)
        pos = 0
        while True:
            idx = region_seq.find(core, pos)
            if idx == -1:
                break
            genome_pos = region_start + idx
            result = _best_site_boundaries(genome_pos, "+")
            if result is not None and result[3] >= 70.0:
                found_sites.append((result[0], result[1], "+", result[2]))
            pos = idx + 1

        # Search reverse strand
        pos = 0
        while True:
            idx = region_seq.find(core_rc, pos)
            if idx == -1:
                break
            genome_pos = region_start + idx
            result = _best_site_boundaries(genome_pos, "-")
            if result is not None and result[3] >= 70.0:
                found_sites.append((result[0], result[1], "-", result[2]))
            pos = idx + 1

        if not found_sites:
            continue

        # For each new site, compute identity against existing arms to
        # confirm it is a genuine sfx recognition site.  The site is stored
        # as an *unpaired* recognition site (not duplicating an existing
        # arm) because cross-type sfx pairings do not map 1:1 onto a
        # single partner arm.
        template_row = cluster.iloc[0].copy()

        for site_start, site_end, strand, site_seq in found_sites:
            best_identity = 0.0

            for _, row in cluster.iterrows():
                if strand == "+":
                    partner_seq = str(row["RightIRSequence"])
                else:
                    partner_seq = str(row["LeftIRSequence"])

                rc_partner = _reverse_complement(partner_seq)
                minlen = min(len(site_seq), len(rc_partner))
                if minlen == 0:
                    continue
                matches = sum(
                    1 for i in range(minlen)
                    if site_seq[i].upper() == rc_partner[i].upper()
                )
                identity = (matches / minlen) * 100.0
                if identity > best_identity:
                    best_identity = identity

            if best_identity < 70.0:
                continue

            # Store as unpaired recognition site.
            # F-sites go in Left columns, R-sites in Right columns;
            # the opposite side is left empty.
            new_row = template_row.copy()
            if strand == "+":
                new_row["LeftIRStart"] = site_start
                new_row["LeftIRStop"] = site_end
                new_row["LeftIRSequence"] = site_seq
                new_row["RightIRStart"] = pd.NA
                new_row["RightIRStop"] = pd.NA
                new_row["RightIRSequence"] = ""
                new_row["InvertibleSequence"] = ""
            else:
                new_row["LeftIRStart"] = pd.NA
                new_row["LeftIRStop"] = pd.NA
                new_row["LeftIRSequence"] = ""
                new_row["RightIRStart"] = site_start
                new_row["RightIRStop"] = site_end
                new_row["RightIRSequence"] = site_seq
                new_row["InvertibleSequence"] = ""

            new_row["PercentIdentity"] = round(best_identity, 2)
            new_row["unpaired_site"] = True
            new_rows.append(new_row)

            logger.info(
                "Cluster %s: motif refinement found unpaired sfx site at "
                "%s:%d-%d (%s) with %.1f%% identity to nearest arm",
                cluster_id, contig, site_start + 1, site_end, strand,
                best_identity,
            )

    # Remove extended pairs that bridge inner boundaries.
    # Re-detect using the same median-ratio logic applied per-cluster above.
    import statistics
    all_extended = set()
    for cluster_id, cluster in ir_df.groupby("cluster_id"):
        if len(cluster) < 3:
            continue
        arm_lens = []
        for _, row in cluster.iterrows():
            arm_lens.append(int(row["LeftIRStop"]) - int(row["LeftIRStart"]))
            arm_lens.append(int(row["RightIRStop"]) - int(row["RightIRStart"]))
        med_a = statistics.median(arm_lens) if arm_lens else 0
        if med_a <= 0:
            continue
        for idx, row in cluster.iterrows():
            l_len = int(row["LeftIRStop"]) - int(row["LeftIRStart"])
            r_len = int(row["RightIRStop"]) - int(row["RightIRStart"])
            if l_len > med_a * 1.4 and r_len > med_a * 1.4:
                all_extended.add(idx)

    # Drop extended pairs from the original DataFrame
    if all_extended:
        ir_df = ir_df.drop(index=all_extended).reset_index(drop=True)
        logger.info(
            "Motif refinement: removed %d extended IR pair(s)",
            len(all_extended),
        )

    if not new_rows and not all_extended:
        logger.info("Motif refinement: no changes")
        return ir_df

    if new_rows:
        new_df = pd.DataFrame(new_rows)
        result = pd.concat([ir_df, new_df], ignore_index=True)
    else:
        result = ir_df

    logger.info(
        "Motif refinement: %d extended pairs removed, %d sites added "
        "(%d IR entries total)",
        len(all_extended), len(new_rows), len(result),
    )
    return result
