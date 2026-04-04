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
COLOUR_INVERTIBLE_SEGMENT = "#F4A261" # sandy orange
COLOUR_CDS_WITH_IR = "#2A9D8F"       # teal
COLOUR_OTHER_CDS = "#CCCCCC"         # light grey

# Distinct palette for IR pairs — each pair gets its own colour.
# Wraps if there are more pairs than palette entries.
IR_PAIR_PALETTE = [
    "#457B9D",  # steel blue
    "#8338EC",  # violet
    "#06D6A0",  # mint
    "#FF6B6B",  # coral
    "#FFD166",  # gold
    "#118AB2",  # cerulean
    "#EF476F",  # pink
    "#073B4C",  # dark teal
    "#B5179E",  # magenta
    "#7209B7",  # purple
]

CATEGORY_COLOURS = {
    "recombinase": COLOUR_RECOMBINASE,
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

import re

_IR_PAIR_RE = re.compile(r"inverted_repeat_(\d+)")


def _extract_ir_pair_number(feature: dict) -> int:
    """Extract the pair number from an IR feature's ID attribute.

    Parses names like ``inverted_repeat_01_FOR`` / ``inverted_repeat_03_REV``
    and returns the integer pair index (1, 3, …).  Returns 0 when the name
    doesn't match the expected pattern.
    """
    name = feature.get("attrs", {}).get("ID", "") or feature.get("attrs", {}).get("Name", "")
    m = _IR_PAIR_RE.search(name)
    return int(m.group(1)) if m else 0


def _ir_pair_colour(pair_number: int) -> str:
    """Return a deterministic colour for the given IR pair number."""
    idx = (pair_number - 1) % len(IR_PAIR_PALETTE) if pair_number > 0 else 0
    return IR_PAIR_PALETTE[idx]


def _strand_int(strand_str: str) -> int:
    if strand_str == "+":
        return +1
    if strand_str == "-":
        return -1
    return 0


def _draw_ir_arrows(ax, ir_data: list[dict], seq_length: int, y_center: float = 0.55):
    """Draw IR features as small directional arrows directly on *ax*.

    All IRs sit on the same horizontal line at *y_center* (in data
    coordinates of the CDS axes).  Each entry in *ir_data* is a dict
    with keys: start, end, strand, colour.

    Arrows point right for ``+``, left for ``-``, and are drawn as
    simple filled polygons (a rectangle body with a triangular tip).
    """
    arrow_h = 0.14  # half-height of the arrow body
    head_len_frac = 0.3  # fraction of feature length used for arrowhead

    for ir in ir_data:
        s, e = ir["start"], ir["end"]
        colour = ir["colour"]
        strand = ir["strand"]
        length = e - s
        head_len = min(length * head_len_frac, seq_length * 0.008)

        if strand == +1:
            # Body rectangle then right-pointing head
            body_end = e - head_len
            verts = [
                (s, y_center - arrow_h),
                (body_end, y_center - arrow_h),
                (body_end, y_center - arrow_h * 1.6),
                (e, y_center),
                (body_end, y_center + arrow_h * 1.6),
                (body_end, y_center + arrow_h),
                (s, y_center + arrow_h),
                (s, y_center - arrow_h),
            ]
        elif strand == -1:
            body_start = s + head_len
            verts = [
                (e, y_center - arrow_h),
                (e, y_center + arrow_h),
                (body_start, y_center + arrow_h),
                (body_start, y_center + arrow_h * 1.6),
                (s, y_center),
                (body_start, y_center - arrow_h * 1.6),
                (body_start, y_center - arrow_h),
                (e, y_center - arrow_h),
            ]
        else:
            # Undirectional rectangle
            verts = [
                (s, y_center - arrow_h),
                (e, y_center - arrow_h),
                (e, y_center + arrow_h),
                (s, y_center + arrow_h),
                (s, y_center - arrow_h),
            ]

        from matplotlib.patches import Polygon
        patch = Polygon(verts, closed=True, facecolor=colour,
                        edgecolor="black", linewidth=0.5, zorder=5)
        ax.add_patch(patch)


def _build_plot(classified: list[dict], seq_length: int, title: str) -> plt.Figure:
    """Build a multi-track matplotlib figure.

    IRs are drawn as small arrow patches directly above the CDS
    backbone line (no separate subplot), so they sit tight against
    the annotation track.  Invertible segments go in a subplot below.

    No inline labels — colour legend only.
    """
    ir_data = []       # dicts with start/end/strand/colour for manual drawing
    cds_features = []  # GraphicFeature objects for dna_features_viewer
    seg_features = []
    ir_pair_numbers = set()

    for f in classified:
        cat = f["category"]
        start = f["start"] - 1  # GFF 1-based -> 0-based
        end = f["end"]

        if cat == "inverted_repeat":
            pair_num = _extract_ir_pair_number(f)
            ir_pair_numbers.add(pair_num)
            ir_data.append({
                "start": start, "end": end,
                "strand": _strand_int(f["strand"]),
                "colour": _ir_pair_colour(pair_num),
            })
        elif cat == "invertible_segment":
            seg_features.append(
                GraphicFeature(
                    start=start, end=end,
                    strand=0, color=CATEGORY_COLOURS[cat],
                    thickness=10, linewidth=0.5,
                )
            )
        else:
            cds_features.append(
                GraphicFeature(
                    start=start, end=end,
                    strand=_strand_int(f["strand"]),
                    color=CATEGORY_COLOURS[cat],
                    thickness=14, linewidth=0.5,
                )
            )

    has_ir = len(ir_data) > 0
    has_seg = len(seg_features) > 0

    # Subplot layout: CDS (+ IRs drawn on same axes), optionally segments below
    n_panels = 1 + int(has_seg)
    ratios = [1.0]
    if has_seg:
        ratios.append(0.4)

    fig_width = max(12, min(20, seq_length / 250))
    fig_height = 1.6 + (0.8 if has_seg else 0)

    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(fig_width, fig_height),
        gridspec_kw={"height_ratios": ratios, "hspace": 0.02},
        sharex=True,
        squeeze=False,
    )
    axes = [ax for ax in axes.flat]
    ax_cds = axes[0]

    # CDS track (main annotation line)
    rec = GraphicRecord(sequence_length=seq_length, features=cds_features)
    rec.plot(ax=ax_cds, with_ruler=not has_seg, draw_line=True, annotate_inline=False)
    ax_cds.set_title(title, fontsize=10, loc="left", fontweight="bold")

    # Draw IR arrows directly above the CDS backbone on the same axes
    if has_ir:
        _draw_ir_arrows(ax_cds, ir_data, seq_length, y_center=0.55)
        # Expand y-limits upward to make room for the IR row
        ymin, ymax = ax_cds.get_ylim()
        ax_cds.set_ylim(ymin, max(ymax, 0.85))

    # Segment track (below CDS)
    if has_seg:
        ax_seg = axes[1]
        rec = GraphicRecord(sequence_length=seq_length, features=seg_features)
        rec.plot(ax=ax_seg, with_ruler=True, draw_line=False, annotate_inline=False)

    # Legend — CDS categories + per-pair IR colours
    handles = []
    for cat in LEGEND_ORDER:
        if cat == "inverted_repeat":
            if not has_ir:
                continue
            for pair_num in sorted(ir_pair_numbers):
                label = f"IR pair {pair_num:02d}" if pair_num > 0 else "IR (unknown pair)"
                handles.append(
                    mpatches.Patch(
                        facecolor=_ir_pair_colour(pair_num),
                        edgecolor="black", linewidth=0.5,
                        label=label,
                    )
                )
        else:
            present = any(f["category"] == cat for f in classified)
            if present:
                handles.append(
                    mpatches.Patch(
                        facecolor=CATEGORY_COLOURS[cat],
                        edgecolor="black", linewidth=0.5,
                        label=LEGEND_LABELS[cat],
                    )
                )

    bottom_ax = axes[-1]
    bottom_ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.4),
        ncol=min(len(handles), 6),
        fontsize=8,
        frameon=False,
    )

    fig.subplots_adjust(hspace=0.02)
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
    """Generate plots for all window GFFs in a directory.

    Scans window_dir for .gff files and writes PNG/SVG plots directly
    into plots_dir (flat, no subdirectories).

    Args:
        window_dir: Directory containing window GFF files.
        plots_dir: Output directory for plot files.

    Returns:
        List of paths to all generated plot files.
    """
    all_outputs = []
    if not os.path.isdir(window_dir):
        return all_outputs
    for fname in sorted(os.listdir(window_dir)):
        if not fname.endswith(".gff"):
            continue
        gff_path = os.path.join(window_dir, fname)
        outputs = generate_shufflon_plot(gff_path, plots_dir)
        all_outputs.extend(outputs)
    return all_outputs
