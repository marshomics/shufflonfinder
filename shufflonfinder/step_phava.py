"""Step 4: Detect inverted repeats via EMBOSS einverted on flanking regions.

Replaces the external PHAVA dependency with direct einverted calls.  Runs
einverted at two sensitivity thresholds (51 and 75), merges the results
with overlap deduplication, filters by minimum inter-arm distance, extracts
arm sequences, and computes percent identity between arms.
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


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class InvertedRepeat:
    """One pair of inverted repeat arms detected by einverted."""
    chrom: str
    left_start: int   # 0-based
    left_stop: int     # 1-based (exclusive-style, as PHAVA used it)
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
    with 0-based start, 1-based stop (matching PHAVA's convention).
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

                # left arm: first token is start, last token is stop
                left_start = int(left_parts[0]) - 1   # to 0-based
                left_stop = int(left_parts[-1])

                # right arm: last token is start (smaller coord), first is stop
                right_start = int(right_parts[-1]) - 1  # to 0-based
                right_stop = int(right_parts[0])

                irs.append((chrom, left_start, left_stop, right_start, right_stop))
                lines = []

    return irs


def _coords_overlap(a: tuple, b: tuple) -> bool:
    """Check if two IR coordinate tuples overlap on the same contig.

    Each tuple is (chrom, left_start, left_stop, right_start, right_stop).
    Overlap is tested across the full span [left_start, right_stop].
    """
    if a[0] != b[0]:
        return False
    a_min, a_max = a[1], a[4]
    b_min, b_max = b[1], b[4]
    return a_max >= b_min and b_max >= a_min


def merge_einverted_results(
    outfile_51: str,
    outfile_75: str,
) -> list[tuple[str, int, int, int, int]]:
    """Merge results from two einverted runs with overlap deduplication.

    All threshold-75 results are kept.  Threshold-51 results are added only
    if they don't overlap any threshold-75 result on the same contig.
    """
    irs_75 = _parse_einverted_outfile(outfile_75)
    irs_51 = _parse_einverted_outfile(outfile_51)

    merged = list(irs_75)

    for ir51 in irs_51:
        if not any(_coords_overlap(ir51, ir75) for ir75 in irs_75):
            merged.append(ir51)

    return merged


# ---------------------------------------------------------------------------
# Sequence extraction and identity computation
# ---------------------------------------------------------------------------

def _compute_percent_identity(seq_a: str, seq_b: str) -> float:
    """Compute percent identity between two sequences of equal length.

    For IR arms the sequences are complementary-reversed, so we compare
    seq_a to the reverse complement of seq_b.  If the lengths differ
    (shouldn't happen with einverted output, but just in case), we
    compare up to the shorter length.
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

    Also applies the minimum inter-arm distance filter (default 30 bp,
    matching PHAVA's hardcoded filter).

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

        # Minimum inter-arm distance filter
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
    """Convert InvertedRepeat objects to a DataFrame matching the old IRs.tsv
    format, plus a PercentIdentity column."""
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
    """Write IRs.tsv to the data/ subdirectory (same layout PHAVA produced)."""
    data_dir = ensure_dir(os.path.join(out_dir, "data"))
    path = os.path.join(data_dir, "IRs.tsv")
    ir_df.to_csv(path, sep="\t", index=False)
    logger.info("Wrote %d IRs to %s", len(ir_df), path)
    return path


# ---------------------------------------------------------------------------
# Top-level entry point (replaces run_phava_on_flanking)
# ---------------------------------------------------------------------------

def detect_inverted_repeats(
    sample: Sample,
    flanking_fasta: str,
    outdir: str,
    cpus: int = 4,
) -> str:
    """Detect inverted repeats in flanking DNA via einverted.

    Replaces the old run_phava_on_flanking(). Runs einverted at two
    thresholds, merges with deduplication, extracts sequences, computes
    percent identity, and writes IRs.tsv.

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

    # Run einverted at two thresholds
    outfile_51, outfile_75 = run_einverted_dual(flanking_fasta, sample_dir, cpus=cpus)

    # Merge with overlap deduplication
    raw_irs = merge_einverted_results(outfile_51, outfile_75)
    logger.info(
        "einverted found %d IRs for %s (after overlap dedup)",
        len(raw_irs), sample.sample_id,
    )

    # Load flanking sequences for arm extraction
    sequences = parse_fasta_file(flanking_fasta)

    # Annotate sequences and compute identity
    annotated = annotate_ir_sequences(raw_irs, sequences)
    logger.info(
        "%d IRs passed min-middle-distance filter for %s",
        len(annotated), sample.sample_id,
    )

    # Export
    ir_df = irs_to_dataframe(annotated)
    export_irs_tsv(ir_df, sample_dir)

    logger.info("IR detection complete for %s", sample.sample_id)
    return sample_dir


# ---------------------------------------------------------------------------
# Loading, remapping, filtering, combining (unchanged interface)
# ---------------------------------------------------------------------------

def load_ir_table(ir_dir: str) -> pd.DataFrame:
    """Load the inverted repeats table from an output directory.

    Looks for data/IRs.tsv within the directory.

    Returns:
        DataFrame with columns: IR_Chr, LeftIRStart, LeftIRStop,
        RightIRStart, RightIRStop, PercentIdentity, etc.
    """
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
    """Convert IR coordinates from flanking-region-local to genome-absolute.

    einverted reports IR positions relative to the flanking FASTA records.
    This function maps them back to the original contig coordinate space.

    Args:
        ir_df: IR table (coordinates relative to flanking records).
        flanking_regions: The FlankingRegion objects used to build the FASTA.

    Returns:
        IR DataFrame with coordinates translated to genome-absolute positions,
        plus additional columns: sample_id, hmm_profiles, locus_tag, contig.
    """
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

        # Offset: flanking region starts at flank_start (1-based)
        offset = region.flank_start - 1  # convert to 0-based for arithmetic

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
    """Filter inverted repeats by arm length and percent identity.

    Arm length is computed from the coordinate columns.
    Percent identity uses the PercentIdentity column (computed from arm
    sequences during detection).  Falls back to searching for legacy column
    names if PercentIdentity isn't present.

    Args:
        ir_df: IR DataFrame (pre- or post-remapping).
        min_arm_length: Minimum arm length in bp.  0 disables the filter.
        min_identity: Minimum percent identity (0-100).  0 disables the filter.

    Returns:
        Filtered DataFrame.
    """
    if ir_df.empty:
        return ir_df

    before = len(ir_df)

    # ---- arm length filter ----
    if min_arm_length > 0:
        left_len = (ir_df["LeftIRStop"] - ir_df["LeftIRStart"]).abs() + 1
        right_len = (ir_df["RightIRStop"] - ir_df["RightIRStart"]).abs() + 1
        mask = (left_len >= min_arm_length) & (right_len >= min_arm_length)
        ir_df = ir_df.loc[mask].copy()
        logger.info(
            "Arm-length filter (>= %d bp): %d -> %d IRs",
            min_arm_length, before, len(ir_df),
        )

    # ---- percent identity filter ----
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
    """Combine and remap IR tables from multiple einverted runs.

    Args:
        ir_results: List of (ir_dir, sample_id, flanking_regions) tuples.
        output_path: Where to write the combined TSV.

    Returns:
        Combined DataFrame with genome-absolute coordinates.
    """
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
