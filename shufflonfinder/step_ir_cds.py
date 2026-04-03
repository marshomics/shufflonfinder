"""Step 4: Integrate inverted repeat coordinates with CDS gene annotations."""

import logging
import os

import pandas as pd

from .sample_sheet import Sample
from .utils import ensure_dir

logger = logging.getLogger("shufflonfinder")


def annotate_genes_with_ir(
    gene_table: pd.DataFrame,
    ir_table: pd.DataFrame,
) -> pd.DataFrame:
    """Flag genes that contain inverted repeat regions within their boundaries.

    For each gene, checks whether either the left or right IR region is
    completely contained within the gene's coordinate span (min/max of
    start/end, handling both strands).

    Args:
        gene_table: DataFrame with columns: contig, start, end, locus_tag, ...
        ir_table: DataFrame with columns: IR_Chr, LeftIRStart, LeftIRStop,
                  RightIRStart, RightIRStop.

    Returns:
        gene_table with an added 'ir' column ('ir' if the gene contains an
        inverted repeat, '' otherwise).
    """
    if ir_table.empty or gene_table.empty:
        gene_table["ir"] = ""
        return gene_table

    gene_table = gene_table.copy()
    gene_table["start"] = pd.to_numeric(gene_table["start"], errors="coerce")
    gene_table["end"] = pd.to_numeric(gene_table["end"], errors="coerce")
    gene_table["gene_lower"] = gene_table[["start", "end"]].min(axis=1)
    gene_table["gene_upper"] = gene_table[["start", "end"]].max(axis=1)

    merged = pd.merge(
        gene_table.reset_index(),
        ir_table,
        how="left",
        left_on="contig",
        right_on="IR_Chr",
    )

    cond_left = (merged["LeftIRStart"] >= merged["gene_lower"]) & (
        merged["LeftIRStop"] <= merged["gene_upper"]
    )
    cond_right = (merged["RightIRStart"] >= merged["gene_lower"]) & (
        merged["RightIRStop"] <= merged["gene_upper"]
    )
    merged["is_ir"] = cond_left | cond_right

    ir_flags = (
        merged.groupby("index")["is_ir"]
        .any()
        .reindex(gene_table.index, fill_value=False)
    )

    gene_table["ir"] = ir_flags.map(lambda x: "ir" if x else "")
    gene_table.drop(columns=["gene_lower", "gene_upper"], inplace=True)
    return gene_table


def run_ir_cds_integration(
    sample: Sample,
    ir_table: pd.DataFrame,
    gene_table_path: str,
    outdir: str,
) -> str:
    """Run IR-CDS integration for one sample.

    Args:
        sample: The Sample being processed.
        ir_table: IR DataFrame (already filtered to this sample if needed).
        gene_table_path: Path to tab-separated gene annotation table.
        outdir: Where to write the annotated table.

    Returns:
        Path to the output annotated table.
    """
    ensure_dir(outdir)

    gene_df = pd.read_csv(gene_table_path, sep="\t", dtype={"contig": str})
    annotated = annotate_genes_with_ir(gene_df, ir_table)

    out_path = os.path.join(outdir, f"{sample.sample_id}_ir_annotated.tsv")
    annotated.to_csv(out_path, sep="\t", index=False)
    logger.info("IR-CDS integration complete for %s -> %s", sample.sample_id, out_path)
    return out_path
