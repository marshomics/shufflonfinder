"""Step 9: KOfamscan annotation of Prokka proteins.

Runs KOfamscan (exec_annotation) on each sample's .faa file to assign
KEGG Orthology (KO) accessions to every protein where possible.  Supports
auto-detection of pre-computed KOfamscan output in either *detail* or
*mapper* format so the run can be skipped when results already exist.
"""

import csv
import logging
import os
import re

import pandas as pd

from .sample_sheet import Sample
from .utils import ensure_dir, run_cmd

logger = logging.getLogger("shufflonfinder")

# ── Column names for the normalised KO table ──────────────────────────
KO_COLUMNS = ["gene_id", "ko_accession", "threshold", "score", "e_value"]


# -----------------------------------------------------------------------
# Format detection
# -----------------------------------------------------------------------

def _detect_kofamscan_format(path: str) -> str:
    """Auto-detect whether a KOfamscan output file is *detail* or *mapper*.

    Detail format starts with a ``#`` header row and has significance marks
    (``*`` or nothing) in the first column.  Mapper format has no header and
    consists of simple ``gene_id\\tKO`` lines.

    Returns:
        ``"detail"`` or ``"mapper"``.
    """
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            # Detail format always starts with a header comment
            if stripped.startswith("#"):
                return "detail"
            # Mapper lines are two tab-separated fields (gene \t KO)
            parts = stripped.split("\t")
            if len(parts) == 2 and parts[1].startswith("K"):
                return "mapper"
            # Detail data lines begin with * or space
            if stripped.startswith("*") or re.match(r"^\s+\S", stripped):
                return "detail"
            break
    # Fall back to detail; the parser will handle errors gracefully
    return "detail"


# -----------------------------------------------------------------------
# Parsers
# -----------------------------------------------------------------------

def _parse_detail_format(path: str) -> pd.DataFrame:
    """Parse KOfamscan *detail* (``--format detail``) output.

    Detail format columns (fixed-width, space-separated):
        significance  gene_name  KO  thrshld  score  E-value  KO_definition

    Only rows marked with ``*`` (above the profile-specific threshold) are
    kept by default.  Rows without a ``*`` are included but flagged with
    ``threshold="below"``.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            # Fixed-width: significance is column 0 (1 char), then whitespace
            sig = line[0].strip()
            parts = line[1:].split()
            if len(parts) < 5:
                continue
            gene_id = parts[0]
            ko = parts[1]
            if not ko.startswith("K"):
                continue
            threshold_status = "above" if sig == "*" else "below"
            score = _safe_float(parts[3])
            e_value = _safe_float(parts[4])
            rows.append({
                "gene_id": gene_id,
                "ko_accession": ko,
                "threshold": threshold_status,
                "score": score,
                "e_value": e_value,
            })
    return pd.DataFrame(rows, columns=KO_COLUMNS)


def _parse_mapper_format(path: str) -> pd.DataFrame:
    """Parse KOfamscan *mapper* output.

    Mapper format: one line per gene–KO assignment, tab-separated.  Only
    above-threshold assignments are emitted by KOfamscan in mapper mode.
    """
    rows = []
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for parts in reader:
            if len(parts) < 2 or parts[0].startswith("#"):
                continue
            gene_id = parts[0].strip()
            ko = parts[1].strip()
            if not ko.startswith("K"):
                continue
            rows.append({
                "gene_id": gene_id,
                "ko_accession": ko,
                "threshold": "above",
                "score": float("nan"),
                "e_value": float("nan"),
            })
    return pd.DataFrame(rows, columns=KO_COLUMNS)


def parse_kofamscan_output(path: str) -> pd.DataFrame:
    """Parse a KOfamscan output file, auto-detecting the format.

    Returns:
        DataFrame with columns: gene_id, ko_accession, threshold, score, e_value.
    """
    fmt = _detect_kofamscan_format(path)
    logger.debug("KOfamscan format detected as '%s' for %s", fmt, path)
    if fmt == "mapper":
        return _parse_mapper_format(path)
    return _parse_detail_format(path)


# -----------------------------------------------------------------------
# Running exec_annotation
# -----------------------------------------------------------------------

def run_kofamscan(
    sample: Sample,
    outdir: str,
    ko_profiles_dir: str,
    cpus: int = 4,
    tmp_dir: str | None = None,
) -> str:
    """Run KOfamscan (exec_annotation) on a sample's protein FASTA.

    Args:
        sample: Sample with a valid ``faa_path``.
        outdir: Directory for the KOfamscan output.
        ko_profiles_dir: Path to the KOfamscan profile + ko_list directory.
            Must contain ``profiles/`` and ``ko_list`` as produced by the
            KOfamscan setup script.
        cpus: Number of parallel workers for hmmsearch.
        tmp_dir: Temporary directory for exec_annotation intermediates.

    Returns:
        Path to the detail-format output file.
    """
    ensure_dir(outdir)
    output_path = os.path.join(outdir, f"{sample.sample_id}_kofamscan.txt")

    profile_dir = os.path.join(ko_profiles_dir, "profiles")
    ko_list = os.path.join(ko_profiles_dir, "ko_list")

    cmd = [
        "exec_annotation",
        "-f", "detail",
        "--profile", profile_dir,
        "--ko-list", ko_list,
        "--cpu", str(cpus),
        "-o", output_path,
    ]
    if tmp_dir:
        ensure_dir(tmp_dir)
        cmd.extend(["--tmp-dir", tmp_dir])

    cmd.append(sample.faa_path)

    run_cmd(cmd, description=f"KOfamscan ({sample.sample_id})")
    return output_path


# -----------------------------------------------------------------------
# Per-sample KO table
# -----------------------------------------------------------------------

def load_ko_table(
    sample: Sample,
    kofamscan_path: str,
    above_threshold_only: bool = True,
) -> pd.DataFrame:
    """Load and normalise a KO table for one sample.

    Args:
        sample: The Sample (used for the sample_id column).
        kofamscan_path: Path to the KOfamscan output file.
        above_threshold_only: If True, keep only above-threshold assignments.

    Returns:
        DataFrame with columns: sample_id, gene_id, ko_accession, threshold,
        score, e_value.
    """
    if not kofamscan_path or not os.path.isfile(kofamscan_path):
        logger.info("No KOfamscan output for %s", sample.sample_id)
        return pd.DataFrame(columns=["sample_id"] + KO_COLUMNS)

    df = parse_kofamscan_output(kofamscan_path)
    if above_threshold_only:
        df = df[df["threshold"] == "above"].copy()

    df.insert(0, "sample_id", sample.sample_id)

    # Deduplicate: keep the best-scoring KO per gene if multiple pass
    if not df.empty and "score" in df.columns:
        df = df.sort_values("score", ascending=False).drop_duplicates(
            subset=["sample_id", "gene_id"], keep="first",
        )

    logger.info(
        "KOfamscan: %d KO assignments for %s (%d above threshold)",
        len(df), sample.sample_id, len(df[df["threshold"] == "above"]),
    )
    return df


# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

def _safe_float(val: str) -> float:
    try:
        return float(val)
    except (ValueError, TypeError):
        return float("nan")
