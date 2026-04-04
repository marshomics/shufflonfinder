"""Step 8: Generate gene-organisation plots for shufflon windows.

Uses dna_features_viewer to produce two-track PNG and SVG figures:
  - Top track: directional arrows for CDS, inverted repeats, and recombinases.
  - Bottom track: undirectional boxes for invertible segments.

Colour scheme is consistent across all windows.
"""

import logging
import os

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for headless rendering
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from dna_features_viewer import GraphicFeature, GraphicRecord

from .utils import ensure_dir

logger = logging.getLogger("shufflonfinder")

# ── Colour scheme ────────────────────────────────────────────────────────────
COLOUR_RECOMBINASE = "#E63946"        # red
COLOUR_INVERTED_REPEAT = "#457B9D"    # steel blue
COLOUR_INVERTIBLE_SEGMENT = "#F4A261" # sandy orange
COLOUR_CDS_WITH_IR = "#2A9D8F"       # teal
COLOUR_OTHER_CDS = "#CCCCCC"         # light grey

CATEGORY_COLOURS = {
    "recombinase": COLOUR_RECOMBINASE,
    "inverted_repeat": COLOUR_INVERTED_REPEAT,
    "invertible_segment": COLOUR_INVERTIBLE_SEGMENT,
    "cds_with_ir": COLOUR_CDS_WITH_IR,
    "other_cds": COLOUR_OTHER_CDS,
}

LEGEND_LABELS = {
    "recombinase": "Recombinase",
    "inverted_repeat": "Inverted repeat",
    "invertible_segment": "Invertible segment",
    "cds_with_ir": "CDS containing IR",
    "other_cds": "CDS",
}


# ── GFF parsing ──────────────────────────────────────────────────────────────

def _parse_gff_attributes(attr_string: str) -> dict[str, str]:
    attrs = {}
    for part in attr_string.split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs


def _parse_window_gff(gff_path: str):
    """Parse a shufflon-window GFF with embedded FASTA.

    Returns (features, sequences) where features is a list of dicts and
    sequences maps seqid -> nucleotide string.
    """
    features = []
    sequences = {}
    in_fasta = False
    current_seq_id = None
    current_seq_parts = []

    with open(gff_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "##FASTA":
                in_fasta = True
                continue
            if in_fasta:
                if line.startswith(">"):
                    if current_seq_id is not None:
                        sequences[current_seq_id] = "".join(current_seq_parts)
                    current_seq_id = line[1:].split()[0]
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line.strip())
                continue
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            attrs = _parse_gff_attributes(parts[8])
            features.append({
                "seqid": parts[0],
                "source": parts[1],
                "type": parts[2],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "score": parts[5],
                "strand": parts[6],
                "phase": parts[7],
                "attributes": parts[8],
                "attrs": attrs,
            })

    if current_seq_id is not None:
        sequences[current_seq_id] = "".join(current_seq_parts)

    return features, sequences


# ── Feature classification ───────────────────────────────────────────────────

def _classify_features(features: list[dict]) -> list[dict]:
    """Assign a category and display label to each feature."""
    ir_spans = [(f["start"], f["end"]) for f in features if f["type"] == "inverted_repeat"]
    hmm_spans = {(f["start"], f["end"]) for f in features if f["type"] == "hmm_hit"}

    classified = []
    for f in features:
        ftype = f["type"]
        source = f["source"]
        attrs = f["attrs"]

        if ftype == "hmm_hit":
            f["category"] = "recombinase"
            gene = attrs.get("gene", "")
            locus = attrs.get("locus_tag", "")
            f["label"] = f"Recombinase ({gene})" if gene else (
                f"Recombinase ({locus})" if locus else "Recombinase"
            )
            classified.append(f)

        elif ftype == "inverted_repeat":
            f["category"] = "inverted_repeat"
            f["label"] = attrs.get("Name", attrs.get("ID", "IR"))
            classified.append(f)

        elif ftype == "invertible_segment":
            f["category"] = "invertible_segment"
            f["label"] = attrs.get("Name", attrs.get("ID", "invertible_segment"))
            classified.append(f)

        elif ftype == "CDS":
            if source == "shufflonfinder":
                continue
            if (f["start"], f["end"]) in hmm_spans:
                continue

            overlaps_ir = any(
                f["start"] <= ir_e and f["end"] >= ir_s
                for ir_s, ir_e in ir_spans
            )
            f["category"] = "cds_with_ir" if overlaps_ir else "other_cds"

            gene = attrs.get("gene", "")
            product = attrs.get("product", "hypothetical protein")
            locus = attrs.get("locus_tag", "")
            if gene:
                f["label"] = gene
            elif product and product != "hypothetical protein":
                f["label"] = product
            elif locus:
                f["label"] = locus
            else:
                f["label"] = "CDS"
            classified.append(f)

    return classified


# ── Plot generation ──────────────────────────────────────────────────────────

def _strand_int(strand_str: str) -> int:
    if strand_str == "+":
        return +1
    if strand_str == "-":
        return -1
    return 0


def _build_plot(classified: list[dict], seq_length: int, title: str) -> plt.Figure:
    """Build a two-track matplotlib figure.

    Top track:  CDS, inverted repeats, recombinase (directional arrows).
    Bottom track: invertible segments (undirectional boxes, no arrows).
    """
    top_features = []
    bottom_features = []

    for f in classified:
        cat = f["category"]
        colour = CATEGORY_COLOURS[cat]
        # GFF is 1-based inclusive; dna_features_viewer is 0-based
        start = f["start"] - 1
        end = f["end"]

        if cat == "invertible_segment":
            bottom_features.append(
                GraphicFeature(
                    start=start,
                    end=end,
                    strand=0,           # undirectional
                    color=colour,
                    label=f["label"],
                    linewidth=0.5,
                    thickness=12,
                )
            )
        else:
            strand = _strand_int(f["strand"])
            top_features.append(
                GraphicFeature(
                    start=start,
                    end=end,
                    strand=strand,
                    color=colour,
                    label=f["label"],
                    linewidth=0.5,
                    thickness=14,
                )
            )

    has_bottom = len(bottom_features) > 0

    # Scale width to sequence length so short windows aren't stretched
    fig_width = max(12, min(20, seq_length / 250))

    if has_bottom:
        fig, (ax_top, ax_bot) = plt.subplots(
            2, 1,
            figsize=(fig_width, 4.5),
            gridspec_kw={"height_ratios": [3, 1.2]},
            sharex=True,
        )
    else:
        fig, ax_top = plt.subplots(1, 1, figsize=(fig_width, 2.5))
        ax_bot = None

    # Top track
    record_top = GraphicRecord(sequence_length=seq_length, features=top_features)
    record_top.plot(
        ax=ax_top,
        with_ruler=not has_bottom,
        draw_line=True,
        annotate_inline=True,
        strand_in_label_threshold=4,
    )
    ax_top.set_title(title, fontsize=10, loc="left", fontweight="bold")

    # Bottom track
    if has_bottom:
        record_bot = GraphicRecord(sequence_length=seq_length, features=bottom_features)
        record_bot.plot(
            ax=ax_bot,
            with_ruler=True,
            draw_line=True,
            annotate_inline=True,
            strand_in_label_threshold=4,
        )

    # Legend
    categories_present = {f["category"] for f in classified}
    # Fixed ordering for consistency
    legend_order = [
        "recombinase", "inverted_repeat", "invertible_segment",
        "cds_with_ir", "other_cds",
    ]
    handles = []
    for cat in legend_order:
        if cat in categories_present:
            handles.append(
                mpatches.Patch(
                    facecolor=CATEGORY_COLOURS[cat],
                    edgecolor="black",
                    linewidth=0.5,
                    label=LEGEND_LABELS[cat],
                )
            )

    legend_ax = ax_bot if has_bottom else ax_top
    legend_ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.4),
        ncol=len(handles),
        fontsize=8,
        frameon=False,
    )

    plt.tight_layout()
    return fig


# ── Main entry points ────────────────────────────────────────────────────────

def generate_shufflon_plot(gff_path: str, output_dir: str) -> list[str]:
    """Generate PNG and SVG plots for a single shufflon-window GFF.

    Args:
        gff_path: Path to a window GFF file (with embedded FASTA).
        output_dir: Directory to write plot files.

    Returns:
        List of paths to generated plot files, or empty list on failure.
    """
    ensure_dir(output_dir)
    basename = os.path.splitext(os.path.basename(gff_path))[0]

    features, sequences = _parse_window_gff(gff_path)
    if not features or not sequences:
        logger.warning("No features or sequences in %s, skipping plot", gff_path)
        return []

    classified = _classify_features(features)
    if not classified:
        logger.warning("No classifiable features in %s, skipping plot", gff_path)
        return []

    # Use first (typically only) sequence's length
    seq_length = len(next(iter(sequences.values())))

    # Build a short title from the basename
    title = basename.replace("_contig_", " — ").replace("_window_", " window ")

    fig = _build_plot(classified, seq_length, title)

    outputs = []
    for ext in ("png", "svg"):
        path = os.path.join(output_dir, f"{basename}.{ext}")
        fig.savefig(path, dpi=200, bbox_inches="tight")
        outputs.append(path)
        logger.info("Plot saved: %s", path)

    plt.close(fig)
    return outputs


def generate_shufflon_plots(window_dir: str, plots_dir: str) -> list[str]:
    """Generate plots for all window GFFs in a directory tree.

    Walks window_dir for .gff files and writes PNG/SVG plots into a
    mirrored subdirectory structure under plots_dir.

    Args:
        window_dir: Top-level 07_shufflon_windows directory containing GFFs.
        plots_dir: Top-level output directory for plot files.

    Returns:
        List of paths to all generated plot files.
    """
    all_outputs = []
    for root, _dirs, files in os.walk(window_dir):
        for fname in sorted(files):
            if not fname.endswith(".gff"):
                continue
            gff_path = os.path.join(root, fname)
            rel = os.path.relpath(root, window_dir)
            out_subdir = os.path.join(plots_dir, rel)
            outputs = generate_shufflon_plot(gff_path, out_subdir)
            all_outputs.extend(outputs)
    return all_outputs
