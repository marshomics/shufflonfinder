"""Step 9: KOfamscan annotation of IR-overlapping CDS.

Extracts protein sequences for CDS that overlap at least one inverted
repeat within shufflon-like or inverton-like windows, runs KOfamscan
(exec_annotation) on that subset, and produces a combined output table
with Prokka annotation, KO accession, and window category.
"""

import csv
import logging
import os
import re

import pandas as pd

from .step_gff import Feature, ShufflonWindow
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
            if stripped.startswith("#"):
                return "detail"
            parts = stripped.split("\t")
            if len(parts) == 2 and parts[1].startswith("K"):
                return "mapper"
            if stripped.startswith("*") or re.match(r"^\s+\S", stripped):
                return "detail"
            break
    return "detail"


# -----------------------------------------------------------------------
# Parsers
# -----------------------------------------------------------------------

def _parse_detail_format(path: str) -> pd.DataFrame:
    """Parse KOfamscan *detail* (``--format detail``) output.

    Detail format columns (fixed-width, space-separated):
        significance  gene_name  KO  thrshld  score  E-value  KO_definition

    Rows marked with ``*`` are above the profile-specific threshold.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
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
    """Parse KOfamscan *mapper* output."""
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
# Identify IR-overlapping CDS within windows
# -----------------------------------------------------------------------

def _cds_overlaps_any_ir(cds: Feature, ir_features: list[Feature]) -> bool:
    """Check if a CDS feature overlaps at least one inverted repeat."""
    for ir in ir_features:
        if cds.overlaps(ir):
            return True
    return False


def identify_ir_cds_in_windows(
    windows: list[ShufflonWindow],
    category: str,
) -> list[dict]:
    """Find CDS with at least one IR overlap in the given windows.

    Args:
        windows: List of ShufflonWindow objects.
        category: ``"shufflon_like"`` or ``"inverton_like"``.

    Returns:
        List of dicts with keys: sample_id, window_id, contig, locus_tag,
        cds_start, cds_end, strand, product, category.
    """
    rows = []
    for win in windows:
        for cds in win.cds_features:
            if not _cds_overlaps_any_ir(cds, win.ir_features):
                continue
            rows.append({
                "sample_id": win.sample_id,
                "window_id": win.window_id,
                "contig": win.contig,
                "locus_tag": cds.locus_tag,
                "cds_start": cds.start + 1,  # 1-based
                "cds_end": cds.end,
                "strand": cds.strand,
                "product": cds.product,
                "category": category,
            })
    return rows


# -----------------------------------------------------------------------
# Extract protein FASTA for IR-overlapping CDS
# -----------------------------------------------------------------------

def _parse_faa(faa_path: str) -> dict[str, str]:
    """Parse a protein FASTA into a dict of id -> sequence."""
    proteins: dict[str, str] = {}
    current_id = ""
    current_seq: list[str] = []
    with open(faa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id:
                    proteins[current_id] = "".join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
    if current_id:
        proteins[current_id] = "".join(current_seq)
    return proteins


def extract_ir_cds_fasta(
    ir_cds_rows: list[dict],
    faa_path: str,
    output_fasta: str,
) -> set[str]:
    """Write a FASTA containing only the IR-overlapping proteins.

    Args:
        ir_cds_rows: List of dicts from identify_ir_cds_in_windows
            (must have 'locus_tag' key).
        faa_path: Path to the full Prokka .faa file.
        output_fasta: Where to write the subset FASTA.

    Returns:
        Set of locus_tags that were found and written.
    """
    locus_tags = {row["locus_tag"] for row in ir_cds_rows}
    if not locus_tags:
        return set()

    proteins = _parse_faa(faa_path)
    written = set()

    with open(output_fasta, "w") as fh:
        for tag in sorted(locus_tags):
            seq = proteins.get(tag, "")
            if seq:
                fh.write(f">{tag}\n{seq}\n")
                written.add(tag)

    if locus_tags - written:
        logger.debug(
            "%d IR-CDS locus tags not found in %s: %s",
            len(locus_tags - written), faa_path,
            ", ".join(sorted(locus_tags - written)[:5]),
        )
    return written


# -----------------------------------------------------------------------
# Run KOfamscan on a protein FASTA
# -----------------------------------------------------------------------

def run_kofamscan(
    fasta_path: str,
    output_path: str,
    ko_profiles_dir: str,
    cpus: int = 4,
    tmp_dir: str | None = None,
    label: str = "",
) -> str:
    """Run KOfamscan (exec_annotation) on a protein FASTA.

    Args:
        fasta_path: Path to input protein FASTA.
        output_path: Where to write the KOfamscan output.
        ko_profiles_dir: Path to the KOfamscan profile + ko_list directory.
        cpus: Number of parallel workers for hmmsearch.
        tmp_dir: Temporary directory for exec_annotation intermediates.
        label: Label for log messages.

    Returns:
        Path to the output file.
    """
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

    cmd.append(fasta_path)

    run_cmd(cmd, description=f"KOfamscan ({label})" if label else "KOfamscan")
    return output_path


# -----------------------------------------------------------------------
# Build the combined IR-CDS KO table
# -----------------------------------------------------------------------

def build_ir_cds_ko_table(
    ir_cds_rows: list[dict],
    kofamscan_path: str | None,
) -> pd.DataFrame:
    """Combine IR-CDS info with KOfamscan results into one table.

    Args:
        ir_cds_rows: Rows from identify_ir_cds_in_windows (both categories).
        kofamscan_path: Path to the KOfamscan output file (may be None).

    Returns:
        DataFrame with columns: sample_id, window_id, contig, locus_tag,
        cds_start, cds_end, strand, product, ko_accession, category.
    """
    columns = [
        "sample_id", "window_id", "contig", "locus_tag",
        "cds_start", "cds_end", "strand", "product",
        "ko_accession", "category",
    ]

    if not ir_cds_rows:
        return pd.DataFrame(columns=columns)

    df = pd.DataFrame(ir_cds_rows)
    # Deduplicate by locus_tag + category (same CDS can appear in multiple windows)
    df = df.drop_duplicates(subset=["locus_tag", "category"])

    # Parse KOfamscan output and join
    ko_lookup: dict[str, str] = {}
    if kofamscan_path and os.path.isfile(kofamscan_path):
        ko_df = parse_kofamscan_output(kofamscan_path)
        # Keep only above-threshold
        ko_df = ko_df[ko_df["threshold"] == "above"]
        # Deduplicate: best score per gene
        if not ko_df.empty:
            ko_df = ko_df.sort_values("score", ascending=False).drop_duplicates(
                subset=["gene_id"], keep="first",
            )
            ko_lookup = dict(zip(ko_df["gene_id"], ko_df["ko_accession"]))

    df["ko_accession"] = df["locus_tag"].map(ko_lookup).fillna("")
    df = df[columns]

    n_with_ko = (df["ko_accession"] != "").sum()
    logger.info(
        "IR-CDS KO table: %d CDS total, %d with KO accession",
        len(df), n_with_ko,
    )
    return df


# -----------------------------------------------------------------------
# Orchestration (per-sample)
# -----------------------------------------------------------------------

def run_ir_cds_kofamscan(
    sample_id: str,
    faa_path: str,
    shufflon_windows: list[ShufflonWindow],
    inverton_windows: list[ShufflonWindow],
    ko_profiles_dir: str,
    outdir: str,
    cpus: int = 4,
) -> pd.DataFrame:
    """Full IR-CDS KOfamscan workflow for one sample.

    1. Identify CDS overlapping IRs in both window types.
    2. Extract those proteins from the .faa.
    3. Run KOfamscan on the extracted subset.
    4. Build and write the combined table.

    Args:
        sample_id: Sample identifier.
        faa_path: Path to the full Prokka .faa.
        shufflon_windows: Shufflon-like windows for this sample.
        inverton_windows: Inverton-like windows for this sample.
        ko_profiles_dir: KOfamscan profile directory.
        outdir: Output directory (within 07_shufflon_windows).
        cpus: Thread count for exec_annotation.

    Returns:
        The combined IR-CDS KO DataFrame.
    """
    ensure_dir(outdir)

    # Step 1: identify IR-overlapping CDS
    shuf_rows = identify_ir_cds_in_windows(shufflon_windows, "shufflon_like")
    inv_rows = identify_ir_cds_in_windows(inverton_windows, "inverton_like")
    all_rows = shuf_rows + inv_rows

    if not all_rows:
        logger.info("No IR-overlapping CDS found for %s", sample_id)
        table = build_ir_cds_ko_table([], None)
        table_path = os.path.join(outdir, f"{sample_id}_ir_cds_ko.tsv")
        table.to_csv(table_path, sep="\t", index=False)
        return table

    # Step 2: extract proteins
    fasta_path = os.path.join(outdir, f"{sample_id}_ir_cds.faa")
    written = extract_ir_cds_fasta(all_rows, faa_path, fasta_path)
    logger.info(
        "Extracted %d IR-CDS proteins for %s (%d unique locus tags)",
        len(written), sample_id,
        len({r["locus_tag"] for r in all_rows}),
    )

    # Step 3: run KOfamscan
    kofamscan_path = None
    if written:
        kofamscan_path = os.path.join(outdir, f"{sample_id}_ir_cds_kofamscan.txt")
        run_kofamscan(
            fasta_path, kofamscan_path, ko_profiles_dir,
            cpus=cpus, label=sample_id,
        )

    # Step 4: build combined table
    table = build_ir_cds_ko_table(all_rows, kofamscan_path)
    table_path = os.path.join(outdir, f"{sample_id}_ir_cds_ko.tsv")
    table.to_csv(table_path, sep="\t", index=False)
    logger.info("Wrote IR-CDS KO table: %s", table_path)

    return table


# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

def _safe_float(val: str) -> float:
    try:
        return float(val)
    except (ValueError, TypeError):
        return float("nan")
