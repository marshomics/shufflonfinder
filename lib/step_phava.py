"""Step 4: Detect inverted repeats with PHAVA on flanking regions around HMM hits."""

import logging
import os

import pandas as pd

from .sample_sheet import Sample
from .step_flanking import FlankingRegion
from .utils import ensure_dir, run_cmd

logger = logging.getLogger("shufflon-pipeline")


def run_phava_on_flanking(
    sample: Sample,
    flanking_fasta: str,
    outdir: str,
    cpus: int = 4,
) -> str:
    """Run PHAVA locate + summarize on the flanking-region FASTA for a sample.

    Instead of scanning the entire genome, this targets only the DNA flanking
    each HMM hit protein.

    Args:
        sample: The Sample being processed.
        flanking_fasta: Path to the multi-record FASTA of flanking regions.
        outdir: Base output directory. PHAVA writes to outdir/<sample_id>/.
        cpus: Number of threads.

    Returns:
        Path to the PHAVA output directory for this sample.
    """
    if not flanking_fasta or not os.path.isfile(flanking_fasta):
        logger.warning("No flanking FASTA for %s, skipping PHAVA", sample.sample_id)
        return ""

    # Check the FASTA isn't empty
    with open(flanking_fasta) as fh:
        if not fh.read(1):
            logger.warning("Empty flanking FASTA for %s, skipping PHAVA", sample.sample_id)
            return ""

    sample_dir = ensure_dir(os.path.join(outdir, sample.sample_id))

    cmd_locate = [
        "phava", "locate",
        "--fasta", flanking_fasta,
        "--dir", sample_dir,
        "--cpus", str(cpus),
    ]
    run_cmd(cmd_locate, description=f"PHAVA locate (flanking): {sample.sample_id}")

    cmd_summarize = [
        "phava", "summarize",
        "--dir", sample_dir,
        "--cpus", str(cpus),
    ]
    run_cmd(cmd_summarize, description=f"PHAVA summarize (flanking): {sample.sample_id}")

    logger.info("PHAVA complete for %s (flanking regions)", sample.sample_id)
    return sample_dir


def load_ir_table(phava_dir: str) -> pd.DataFrame:
    """Load the inverted repeats table from a PHAVA output directory.

    Looks for data/IRs.tsv within the PHAVA output directory.

    Returns:
        DataFrame with columns: IR_Chr, LeftIRStart, LeftIRStop,
        RightIRStart, RightIRStop, etc.
    """
    if not phava_dir:
        return pd.DataFrame()

    ir_path = os.path.join(phava_dir, "data", "IRs.tsv")
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

    PHAVA reports IR positions relative to the flanking FASTA records.
    This function maps them back to the original contig coordinate space.

    Args:
        ir_df: IR table from PHAVA (coordinates relative to flanking records).
        flanking_regions: The FlankingRegion objects used to build the FASTA.

    Returns:
        IR DataFrame with coordinates translated to genome-absolute positions,
        plus additional columns: sample_id, hmm_profile, locus_tag, contig.
    """
    if ir_df.empty:
        return ir_df

    # Build a lookup from hit_id -> FlankingRegion
    region_map = {r.hit_id: r for r in flanking_regions}

    rows = []
    for _, ir_row in ir_df.iterrows():
        ir_chr = str(ir_row["IR_Chr"])

        # The IR_Chr should match a hit_id from the flanking FASTA headers
        region = region_map.get(ir_chr)
        if region is None:
            logger.debug("IR_Chr '%s' doesn't match any flanking region, skipping", ir_chr)
            continue

        # Offset: flanking region starts at flank_start (1-based)
        offset = region.flank_start - 1  # convert to 0-based for arithmetic

        remapped = ir_row.copy()
        remapped["IR_Chr"] = region.contig  # Map back to original contig
        remapped["LeftIRStart"] = ir_row["LeftIRStart"] + offset
        remapped["LeftIRStop"] = ir_row["LeftIRStop"] + offset
        remapped["RightIRStart"] = ir_row["RightIRStart"] + offset
        remapped["RightIRStop"] = ir_row["RightIRStop"] + offset

        # Add metadata
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


def combine_ir_tables(
    phava_results: list[tuple[str, str, list[FlankingRegion]]],
    output_path: str,
) -> pd.DataFrame:
    """Combine and remap IR tables from multiple PHAVA runs.

    Args:
        phava_results: List of (phava_dir, sample_id, flanking_regions) tuples.
        output_path: Where to write the combined TSV.

    Returns:
        Combined DataFrame with genome-absolute coordinates.
    """
    frames = []
    for pdir, sid, flanking_regions in phava_results:
        raw_ir = load_ir_table(pdir)
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
