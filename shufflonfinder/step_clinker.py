"""Step 8: Generate gene-organisation plots for shufflon windows.

Uses dna_features_viewer to produce three-track PNG and SVG figures:
  - Top track: inverted repeat arrows (above the annotation line).
  - Middle track: CDS and recombinase arrows (main annotation line).
  - Bottom track: undirectional invertible segment boxes (below).

Features are colour-coded; identification is via legend only (no inline
labels). Colour scheme is consistent across all windows.
"""

import logging
import os

import matplotlib
matplotlib.use("Agg")
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

# Fixed legend ordering so every plot is consistent
LEGEND_ORDER = [
    "recombinase", "inverted_repeat", "invertible_segment",
    "cds_with_ir", "other_cds",
]


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
    """Assign a category to each feature."""
    ir_spans = [(f["start"], f["end"]) for f in features if f["type"] == "inverted_repeat"]
    hmm_spans = {(f["start"], f["end"]) for f in features if f["type"] == "hmm_hit"}

    classified = []
    for f in features:
        ftype = f["type"]
        source = f["source"]

        if ftype == "hmm_hit":
            f["category"] = "recombinase"
            classified.append(f)

        elif ftype == "inverted_repeat":
            f["category"] = "inverted_repeat"
            classified.append(f)

        elif ftype == "invertible_segment":
            f["category"] = "invertible_segment"
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
    """Build a three-track matplotlib figure.

    Top track:    inverted repeats (directional arrows, no line).
    Middle track: CDS + recombinase (directional arrows, with backbone line).
    Bottom track: invertible segments (undirectional boxes, no line).

    No inline labels — colour legend only.
    """
    ir_features = []
    cds_features = []
    seg_features = []

    for f in classified:
        cat = f["category"]
        colour = CATEGORY_COLOURS[cat]
        start = f["start"] - 1  # GFF 1-based -> 0-based
        end = f["end"]

        if cat == "inverted_repeat":
            ir_features.append(
                GraphicFeature(
                    start=start, end=end,
                    strand=_strand_int(f["strand"]),
                    color=colour, thickness=8, linewidth=0.5,
                )
            )
        elif cat == "invertible_segment":
            seg_features.append(
                GraphicFeature(
                    start=start, end=end,
                    strand=0, color=colour,
                    thickness=10, linewidth=0.5,
                )
            )
        else:
            cds_features.append(
                GraphicFeature(
                    start=start, end=end,
                    strand=_strand_int(f["strand"]),
                    color=colour, thickness=14, linewidth=0.5,
                )
            )

    has_ir = len(ir_features) > 0
    has_seg = len(seg_features) > 0

    # Build subplot layout depending on which tracks are populated
    ratios = []
    track_keys = []
    if has_ir:
        ratios.append(0.6)
        track_keys.append("ir")
    ratios.append(1.2)
    track_keys.append("cds")
    if has_seg:
        ratios.append(0.6)
        track_keys.append("seg")

    fig_width = max(12, min(20, seq_length / 250))
    fig_height = 1.0 + 0.9 * len(ratios)

    fig, axes = plt.subplots(
        len(ratios), 1,
        figsize=(fig_width, fig_height),
        gridspec_kw={"height_ratios": ratios, "hspace": 0.05},
        sharex=True,
        squeeze=False,
    )
    axes = [ax for ax in axes.flat]

    ax_map = dict(zip(track_keys, axes))

    # IR track (above CDS)
    if "ir" in ax_map:
        rec = GraphicRecord(sequence_length=seq_length, features=ir_features)
        rec.plot(ax=ax_map["ir"], with_ruler=False, draw_line=False, annotate_inline=False)

    # CDS track (main annotation line)
    rec = GraphicRecord(sequence_length=seq_length, features=cds_features)
    rec.plot(ax=ax_map["cds"], with_ruler=not has_seg, draw_line=True, annotate_inline=False)
    ax_map["cds"].set_title(title, fontsize=10, loc="left", fontweight="bold")

    # Segment track (below CDS)
    if "seg" in ax_map:
        rec = GraphicRecord(sequence_length=seq_length, features=seg_features)
        rec.plot(ax=ax_map["seg"], with_ruler=True, draw_line=False, annotate_inline=False)

    # Legend
    categories_present = {f["category"] for f in classified}
    handles = [
        mpatches.Patch(
            facecolor=CATEGORY_COLOURS[cat],
            edgecolor="black", linewidth=0.5,
            label=LEGEND_LABELS[cat],
        )
        for cat in LEGEND_ORDER if cat in categories_present
    ]

    bottom_ax = axes[-1]
    bottom_ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.5),
        ncol=len(handles),
        fontsize=8,
        frameon=False,
    )

    fig.subplots_adjust(hspace=0.05)
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

    seq_length = len(next(iter(sequences.values())))
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

    Walks window_dir for .gff files and writes PNG/SVG plots into
    per-sample subdirectories under plots_dir, mirroring the
    subdirectory structure from window_dir.

    Args:
        window_dir: Top-level GFF directory (e.g. 07_shufflon_windows/gffs/).
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
