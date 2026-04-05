"""Step 10: KEGG enrichment analysis for IR-containing CDS.

Identifies CDS that overlap inverted repeats within shufflon-like and
inverton-like windows, maps them to KEGG Orthology (KO) accessions via
KOfamscan output, and performs hypergeometric enrichment tests for KEGG
pathways and modules among inverton-associated CDS.  Produces per-sample
output tables and horizontal bar plots of enriched terms.
"""

import logging
import math
import os
from collections import defaultdict
from io import StringIO
from urllib.request import urlopen
from urllib.error import URLError

import pandas as pd

from .step_gff import Feature, ShufflonWindow
from .utils import ensure_dir

logger = logging.getLogger("shufflonfinder")


# -----------------------------------------------------------------------
# KO-to-pathway / KO-to-module mapping via KEGG REST API
# -----------------------------------------------------------------------

_KO_PATHWAY_CACHE: dict[str, set[str]] | None = None
_KO_MODULE_CACHE: dict[str, set[str]] | None = None
_PATHWAY_NAME_CACHE: dict[str, str] = {}
_MODULE_NAME_CACHE: dict[str, str] = {}


def _fetch_kegg_link(target: str, source: str) -> dict[str, set[str]]:
    """Fetch KEGG LINK data (source -> set of targets).

    Uses the KEGG REST API: https://rest.kegg.jp/link/{target}/{source}
    """
    url = f"https://rest.kegg.jp/link/{target}/{source}"
    mapping: dict[str, set[str]] = defaultdict(set)
    try:
        logger.info("Fetching KEGG link data: %s -> %s", source, target)
        with urlopen(url, timeout=60) as resp:
            for line in resp.read().decode("utf-8").splitlines():
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    src_id = parts[0].split(":")[-1]  # e.g. ko:K00001 -> K00001
                    tgt_id = parts[1].split(":")[-1]  # e.g. path:ko00010 -> ko00010
                    mapping[src_id].add(tgt_id)
    except (URLError, OSError) as exc:
        logger.warning("Could not fetch KEGG link (%s -> %s): %s", source, target, exc)
    return dict(mapping)


def _fetch_kegg_list(database: str) -> dict[str, str]:
    """Fetch KEGG LIST data (id -> name).

    Uses the KEGG REST API: https://rest.kegg.jp/list/{database}
    """
    url = f"https://rest.kegg.jp/list/{database}"
    names: dict[str, str] = {}
    try:
        logger.info("Fetching KEGG list: %s", database)
        with urlopen(url, timeout=60) as resp:
            for line in resp.read().decode("utf-8").splitlines():
                parts = line.strip().split("\t", 1)
                if len(parts) == 2:
                    entry_id = parts[0].split(":")[-1]
                    names[entry_id] = parts[1]
    except (URLError, OSError) as exc:
        logger.warning("Could not fetch KEGG list (%s): %s", database, exc)
    return names


def get_ko_pathway_map() -> dict[str, set[str]]:
    """Return KO -> set of pathway IDs mapping (cached)."""
    global _KO_PATHWAY_CACHE
    if _KO_PATHWAY_CACHE is None:
        _KO_PATHWAY_CACHE = _fetch_kegg_link("pathway", "ko")
    return _KO_PATHWAY_CACHE


def get_ko_module_map() -> dict[str, set[str]]:
    """Return KO -> set of module IDs mapping (cached)."""
    global _KO_MODULE_CACHE
    if _KO_MODULE_CACHE is None:
        _KO_MODULE_CACHE = _fetch_kegg_link("module", "ko")
    return _KO_MODULE_CACHE


def get_pathway_names() -> dict[str, str]:
    """Return pathway ID -> name mapping (cached)."""
    global _PATHWAY_NAME_CACHE
    if not _PATHWAY_NAME_CACHE:
        _PATHWAY_NAME_CACHE = _fetch_kegg_list("pathway")
    return _PATHWAY_NAME_CACHE


def get_module_names() -> dict[str, str]:
    """Return module ID -> name mapping (cached)."""
    global _MODULE_NAME_CACHE
    if not _MODULE_NAME_CACHE:
        _MODULE_NAME_CACHE = _fetch_kegg_list("module")
    return _MODULE_NAME_CACHE


# -----------------------------------------------------------------------
# Identify CDS overlapping IRs within windows
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
    ko_df: pd.DataFrame,
) -> pd.DataFrame:
    """Find CDS with at least one IR overlap in the given windows and map to KO.

    Args:
        windows: List of ShufflonWindow objects.
        category: ``"shufflon_like"`` or ``"inverton_like"``.
        ko_df: KO table for this sample (columns: sample_id, gene_id, ko_accession).

    Returns:
        DataFrame with columns: sample_id, window_id, contig, locus_tag,
        cds_start, cds_end, strand, product, ko_accession, category.
    """
    if not windows:
        return pd.DataFrame(columns=[
            "sample_id", "window_id", "contig", "locus_tag",
            "cds_start", "cds_end", "strand", "product",
            "ko_accession", "category",
        ])

    # Build gene_id -> ko lookup
    ko_lookup: dict[str, str] = {}
    if not ko_df.empty:
        for _, row in ko_df.iterrows():
            ko_lookup[row["gene_id"]] = row["ko_accession"]

    rows = []
    for win in windows:
        for cds in win.cds_features:
            if not _cds_overlaps_any_ir(cds, win.ir_features):
                continue
            tag = cds.locus_tag
            rows.append({
                "sample_id": win.sample_id,
                "window_id": win.window_id,
                "contig": win.contig,
                "locus_tag": tag,
                "cds_start": cds.start + 1,  # 1-based
                "cds_end": cds.end,
                "strand": cds.strand,
                "product": cds.product,
                "ko_accession": ko_lookup.get(tag, ""),
                "category": category,
            })

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------
# Hypergeometric enrichment test
# -----------------------------------------------------------------------

def _hypergeom_pvalue(k: int, M: int, n: int, N: int) -> float:
    """Compute hypergeometric survival function P(X >= k).

    Args:
        k: Number of successes in sample (observed).
        M: Population size.
        n: Number of successes in population.
        N: Sample size (number drawn).

    Uses scipy if available, otherwise a pure-Python fallback.
    """
    try:
        from scipy.stats import hypergeom
        # P(X >= k) = 1 - P(X <= k-1) = hypergeom.sf(k-1, M, n, N)
        return float(hypergeom.sf(k - 1, M, n, N))
    except ImportError:
        return _hypergeom_pvalue_fallback(k, M, n, N)


def _hypergeom_pvalue_fallback(k: int, M: int, n: int, N: int) -> float:
    """Pure-Python hypergeometric survival function using log-space math."""
    def _log_comb(a, b):
        if b < 0 or b > a:
            return float("-inf")
        if b == 0 or b == a:
            return 0.0
        b = min(b, a - b)
        result = 0.0
        for i in range(b):
            result += math.log(a - i) - math.log(i + 1)
        return result

    log_total = _log_comb(M, N)
    pval = 0.0
    for x in range(k, min(n, N) + 1):
        log_p = _log_comb(n, x) + _log_comb(M - n, N - x) - log_total
        pval += math.exp(log_p)
    return min(pval, 1.0)


def _bh_correction(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    adjusted = [0.0] * n
    cummin = 1.0
    for rank_idx in range(n - 1, -1, -1):
        orig_idx, pval = indexed[rank_idx]
        rank = rank_idx + 1
        adj = pval * n / rank
        cummin = min(cummin, adj)
        adjusted[orig_idx] = min(cummin, 1.0)
    return adjusted


def kegg_enrichment_test(
    ir_ko_accessions: set[str],
    background_ko_accessions: set[str],
    term_type: str = "pathway",
) -> pd.DataFrame:
    """Run hypergeometric enrichment for KEGG pathways or modules.

    Args:
        ir_ko_accessions: KO accessions of IR-overlapping CDS.
        background_ko_accessions: KO accessions of all CDS in the genome.
        term_type: ``"pathway"`` or ``"module"``.

    Returns:
        DataFrame with columns: term_id, term_name, k, M, n, N,
        fold_change, pvalue, padj, neg_log10_pval.
    """
    if not ir_ko_accessions or not background_ko_accessions:
        return pd.DataFrame(columns=[
            "term_id", "term_name", "k", "M", "n", "N",
            "fold_change", "pvalue", "padj", "neg_log10_pval",
        ])

    if term_type == "pathway":
        ko_term_map = get_ko_pathway_map()
        term_names = get_pathway_names()
    else:
        ko_term_map = get_ko_module_map()
        term_names = get_module_names()

    # Population: all background KOs
    M = len(background_ko_accessions)
    # Sample: IR-overlapping KOs
    N = len(ir_ko_accessions)

    if M == 0 or N == 0:
        return pd.DataFrame(columns=[
            "term_id", "term_name", "k", "M", "n", "N",
            "fold_change", "pvalue", "padj", "neg_log10_pval",
        ])

    # Count KOs per term in background and in IR set
    term_bg_count: dict[str, int] = defaultdict(int)
    term_ir_count: dict[str, int] = defaultdict(int)

    for ko in background_ko_accessions:
        for term in ko_term_map.get(ko, set()):
            term_bg_count[term] += 1

    for ko in ir_ko_accessions:
        for term in ko_term_map.get(ko, set()):
            term_ir_count[term] += 1

    rows = []
    for term_id, k in term_ir_count.items():
        if k < 1:
            continue
        n = term_bg_count.get(term_id, 0)
        if n == 0:
            continue

        pval = _hypergeom_pvalue(k, M, n, N)
        expected = (n / M) * N if M > 0 else 0
        fold_change = k / expected if expected > 0 else float("inf")

        rows.append({
            "term_id": term_id,
            "term_name": term_names.get(term_id, term_id),
            "k": k,
            "M": M,
            "n": n,
            "N": N,
            "fold_change": fold_change,
            "pvalue": pval,
        })

    if not rows:
        return pd.DataFrame(columns=[
            "term_id", "term_name", "k", "M", "n", "N",
            "fold_change", "pvalue", "padj", "neg_log10_pval",
        ])

    df = pd.DataFrame(rows)
    df["padj"] = _bh_correction(df["pvalue"].tolist())
    df["neg_log10_pval"] = df["pvalue"].apply(
        lambda p: -math.log10(p) if p > 0 else 300.0
    )
    df = df.sort_values("pvalue").reset_index(drop=True)
    return df


# -----------------------------------------------------------------------
# Bar plots
# -----------------------------------------------------------------------

def plot_enrichment_bars(
    enrichment_df: pd.DataFrame,
    output_path: str,
    title: str = "KEGG enrichment",
    max_terms: int = 25,
    pval_cutoff: float = 0.05,
) -> str | None:
    """Produce a horizontal bar plot of enriched KEGG terms.

    Bars show ``-log10(p-value) × fold-change`` for each term.  Only terms
    with ``pvalue < pval_cutoff`` are plotted.  Bar colour encodes fold
    change.

    Args:
        enrichment_df: Output of :func:`kegg_enrichment_test`.
        output_path: Path for the output PNG file.
        title: Plot title.
        max_terms: Maximum number of terms to show.
        pval_cutoff: Only plot terms below this raw p-value.

    Returns:
        Path to the plot file, or ``None`` if nothing to plot.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
    except ImportError:
        logger.warning("matplotlib not available; skipping enrichment plot")
        return None

    sig = enrichment_df[enrichment_df["pvalue"] < pval_cutoff].copy()
    if sig.empty:
        logger.info("No terms below p-value cutoff %.2g for %s", pval_cutoff, title)
        return None

    sig = sig.head(max_terms).copy()
    sig["score"] = sig["neg_log10_pval"] * sig["fold_change"]
    sig = sig.sort_values("score", ascending=True)  # bottom-to-top

    # Truncate long term names
    sig["label"] = sig["term_name"].apply(
        lambda s: s[:60] + "..." if len(s) > 63 else s
    )

    fig_height = max(3, 0.35 * len(sig) + 1.5)
    fig, ax = plt.subplots(figsize=(8, fig_height))

    norm = Normalize(vmin=sig["fold_change"].min(), vmax=sig["fold_change"].max())
    cmap = plt.cm.YlOrRd
    colours = [cmap(norm(fc)) for fc in sig["fold_change"]]

    bars = ax.barh(range(len(sig)), sig["score"], color=colours, edgecolor="none")
    ax.set_yticks(range(len(sig)))
    ax.set_yticklabels(sig["label"], fontsize=8)
    ax.set_xlabel("-log₁₀(p-value) × fold change", fontsize=10)
    ax.set_title(title, fontsize=11, pad=10)

    # Colour bar for fold change
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, aspect=30)
    cbar.set_label("Fold change", fontsize=9)

    # Annotate bars with k/n counts
    for i, (_, row) in enumerate(sig.iterrows()):
        ax.text(
            row["score"] + sig["score"].max() * 0.01, i,
            f'{row["k"]}/{row["n"]}',
            va="center", fontsize=7, color="grey",
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    ensure_dir(os.path.dirname(output_path))
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Also save SVG
    svg_path = output_path.rsplit(".", 1)[0] + ".svg"
    fig2, ax2 = plt.subplots(figsize=(8, fig_height))
    bars2 = ax2.barh(range(len(sig)), sig["score"], color=colours, edgecolor="none")
    ax2.set_yticks(range(len(sig)))
    ax2.set_yticklabels(sig["label"], fontsize=8)
    ax2.set_xlabel("-log₁₀(p-value) × fold change", fontsize=10)
    ax2.set_title(title, fontsize=11, pad=10)
    sm2 = ScalarMappable(cmap=cmap, norm=norm)
    sm2.set_array([])
    cbar2 = fig2.colorbar(sm2, ax=ax2, pad=0.02, aspect=30)
    cbar2.set_label("Fold change", fontsize=9)
    for i, (_, row) in enumerate(sig.iterrows()):
        ax2.text(
            row["score"] + sig["score"].max() * 0.01, i,
            f'{row["k"]}/{row["n"]}',
            va="center", fontsize=7, color="grey",
        )
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    plt.tight_layout()
    fig2.savefig(svg_path, bbox_inches="tight")
    plt.close(fig2)

    logger.info("Enrichment plot: %s", output_path)
    return output_path


# -----------------------------------------------------------------------
# Main orchestration (per-sample)
# -----------------------------------------------------------------------

def run_kegg_enrichment_for_sample(
    sample_id: str,
    shufflon_windows: list[ShufflonWindow],
    inverton_windows: list[ShufflonWindow],
    ko_df: pd.DataFrame,
    background_ko_accessions: set[str],
    outdir: str,
) -> tuple[pd.DataFrame, pd.DataFrame | None, pd.DataFrame | None]:
    """Run the full KEGG enrichment workflow for one sample.

    1. Identify IR-overlapping CDS in both window types.
    2. Map CDS to KO accessions.
    3. For inverton-like windows, run hypergeometric tests against the
       genome-wide background.
    4. Produce enrichment bar plots.

    Args:
        sample_id: Sample identifier.
        shufflon_windows: Shufflon-like windows for this sample.
        inverton_windows: Inverton-like windows for this sample.
        ko_df: Full KO table for this sample.
        background_ko_accessions: Set of all KO accessions in this genome.
        outdir: Output directory for tables and plots.

    Returns:
        Tuple of (ir_cds_ko_table, pathway_enrichment_df, module_enrichment_df).
        Enrichment DataFrames are None if no enrichment was possible.
    """
    ensure_dir(outdir)
    plot_dir = ensure_dir(os.path.join(outdir, "plots"))

    # Step 1: identify IR-overlapping CDS and their KOs
    shuf_ir_cds = identify_ir_cds_in_windows(shufflon_windows, "shufflon_like", ko_df)
    inv_ir_cds = identify_ir_cds_in_windows(inverton_windows, "inverton_like", ko_df)
    ir_cds_table = pd.concat([shuf_ir_cds, inv_ir_cds], ignore_index=True)

    # Deduplicate by locus_tag (same CDS can appear in multiple windows)
    if not ir_cds_table.empty:
        ir_cds_table = ir_cds_table.drop_duplicates(subset=["locus_tag", "category"])

    # Write the combined table
    table_path = os.path.join(outdir, f"{sample_id}_ir_cds_ko.tsv")
    ir_cds_table.to_csv(table_path, sep="\t", index=False)
    logger.info(
        "IR-CDS KO table for %s: %d CDS (%d with KO)",
        sample_id, len(ir_cds_table),
        ir_cds_table["ko_accession"].astype(bool).sum(),
    )

    # Step 2: enrichment for inverton-like windows only
    pathway_df = None
    module_df = None

    inv_kos = set(inv_ir_cds["ko_accession"].dropna()) - {""}
    if inv_kos and background_ko_accessions:
        # Pathway enrichment
        pathway_df = kegg_enrichment_test(inv_kos, background_ko_accessions, "pathway")
        if not pathway_df.empty:
            pw_path = os.path.join(outdir, f"{sample_id}_inverton_pathway_enrichment.tsv")
            pathway_df.to_csv(pw_path, sep="\t", index=False)
            plot_enrichment_bars(
                pathway_df,
                os.path.join(plot_dir, f"{sample_id}_inverton_pathway_enrichment.png"),
                title=f"{sample_id} — Inverton KEGG pathway enrichment",
            )

        # Module enrichment
        module_df = kegg_enrichment_test(inv_kos, background_ko_accessions, "module")
        if not module_df.empty:
            mod_path = os.path.join(outdir, f"{sample_id}_inverton_module_enrichment.tsv")
            module_df.to_csv(mod_path, sep="\t", index=False)
            plot_enrichment_bars(
                module_df,
                os.path.join(plot_dir, f"{sample_id}_inverton_module_enrichment.png"),
                title=f"{sample_id} — Inverton KEGG module enrichment",
            )
    else:
        logger.info(
            "Skipping KEGG enrichment for %s: %d inverton KOs, %d background KOs",
            sample_id, len(inv_kos), len(background_ko_accessions),
        )

    return ir_cds_table, pathway_df, module_df
