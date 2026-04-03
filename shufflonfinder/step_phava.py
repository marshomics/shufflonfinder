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
from .step_flanking import FlankingRegion, parse_fasta_file, parse_cds_from_gff
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

    sample_dir = ensure_dir(os.path.join(outdir, sample.sample_id))

    outfile_51, outfile_75 = run_einverted_dual(flanking_fasta, sample_dir, cpus=cpus)

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
    export_irs_tsv(ir_df, sample_dir)

    logger.info("IR detection complete for %s", sample.sample_id)
    return sample_dir


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


def _arm_overlaps_cds(
    arm_start: int,
    arm_stop: int,
    cds_list: list,
) -> bool:
    """Check if an IR arm overlaps any CDS feature.

    Uses 0-based coordinates throughout.  An arm overlaps a CDS if
    [arm_start, arm_stop) intersects [cds.start-1, cds.end).
    """
    for cds in cds_list:
        cds_s = cds.start - 1  # CdsCoords are 1-based
        cds_e = cds.end
        if arm_stop > cds_s and arm_start < cds_e:
            return True
    return False


def filter_shufflon_candidates(
    ir_df: pd.DataFrame,
    samples: list[Sample],
    window_size: int = 3000,
    min_ir_pairs: int = 2,
) -> pd.DataFrame:
    """Filter IRs to keep only dense, shufflon-like clusters.

    A valid shufflon candidate must have:
      1. At least `min_ir_pairs` IR pairs within a `window_size` bp cluster.
      2. At least one IR arm in the cluster overlapping a CDS.

    IRs that don't belong to any qualifying cluster are removed.

    Args:
        ir_df: Remapped IR DataFrame (genome-absolute coordinates).
        samples: List of all Sample objects (for GFF CDS lookup).
        window_size: Maximum gap (bp) for clustering nearby IR pairs.
        min_ir_pairs: Minimum number of IR pairs per cluster.

    Returns:
        Filtered DataFrame containing only IRs from qualifying clusters,
        plus a new ``cluster_id`` column.
    """
    if ir_df.empty:
        return ir_df

    # Build a CDS lookup per (sample, contig)
    sample_map = {s.sample_id: s for s in samples}
    cds_cache: dict[str, dict[str, list]] = {}  # sample_id -> contig -> [CdsCoords]

    def _get_cds(sample_id: str, contig: str) -> list:
        if sample_id not in cds_cache:
            sample = sample_map.get(sample_id)
            if sample and sample.gff_path:
                by_contig: dict[str, list] = {}
                for cds in parse_cds_from_gff(sample.gff_path).values():
                    by_contig.setdefault(cds.contig, []).append(cds)
                cds_cache[sample_id] = by_contig
            else:
                cds_cache[sample_id] = {}
        return cds_cache[sample_id].get(contig, [])

    before = len(ir_df)
    keep_indices = []
    cluster_labels = {}

    # Process each sample independently
    for sample_id, sample_irs in ir_df.groupby("sample_id"):
        clusters = _cluster_ir_rows(sample_irs, window_size)

        for ci, cluster_idx in enumerate(clusters):
            # Density check
            if len(cluster_idx) < min_ir_pairs:
                continue

            # CDS overlap check: at least one arm in the cluster must
            # overlap a CDS on the same contig
            has_cds_overlap = False
            for idx in cluster_idx:
                row = ir_df.loc[idx]
                contig = str(row["IR_Chr"])
                cds_list = _get_cds(sample_id, contig)
                if (_arm_overlaps_cds(int(row["LeftIRStart"]), int(row["LeftIRStop"]), cds_list)
                        or _arm_overlaps_cds(int(row["RightIRStart"]), int(row["RightIRStop"]), cds_list)):
                    has_cds_overlap = True
                    break

            if not has_cds_overlap:
                continue

            cid = f"{sample_id}_cluster_{ci + 1}"
            for idx in cluster_idx:
                keep_indices.append(idx)
                cluster_labels[idx] = cid

    if not keep_indices:
        logger.warning("No IR clusters met the shufflon candidate criteria")
        return pd.DataFrame(columns=list(ir_df.columns) + ["cluster_id"])

    result = ir_df.loc[keep_indices].copy()
    result["cluster_id"] = [cluster_labels[i] for i in keep_indices]

    removed = before - len(result)
    n_clusters = result["cluster_id"].nunique()
    logger.info(
        "Shufflon candidate filter: %d -> %d IRs in %d cluster(s) "
        "(removed %d IRs from sparse or non-coding regions)",
        before, len(result), n_clusters, removed,
    )
    return result
