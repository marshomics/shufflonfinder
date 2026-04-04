"""Step 8: Generate Clinker gene-cluster comparison plots for shufflon windows."""

import csv
import logging
import os
import subprocess
import tempfile
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .utils import ensure_dir

logger = logging.getLogger("shufflonfinder")

# ── Colour scheme ────────────────────────────────────────────────────────────
# Consistent across all plots.
COLOUR_RECOMBINASE = "#E63946"    # red
COLOUR_INVERTED_REPEAT = "#457B9D"  # steel blue
COLOUR_INVERTIBLE_SEGMENT = "#F4A261"  # sandy orange
COLOUR_CDS_WITH_IR = "#2A9D8F"   # teal
COLOUR_OTHER_CDS = "#CCCCCC"     # light grey

CATEGORY_LABELS = {
    "recombinase": "Recombinase",
    "inverted_repeat": "Inverted repeat",
    "invertible_segment": "Invertible segment",
    "cds_with_ir": "CDS containing inverted repeat",
    "other_cds": "CDS",
}

CATEGORY_COLOURS = {
    "recombinase": COLOUR_RECOMBINASE,
    "inverted_repeat": COLOUR_INVERTED_REPEAT,
    "invertible_segment": COLOUR_INVERTIBLE_SEGMENT,
    "cds_with_ir": COLOUR_CDS_WITH_IR,
    "other_cds": COLOUR_OTHER_CDS,
}


# ── GFF parsing helpers ─────────────────────────────────────────────────────

def _parse_gff_attributes(attr_string: str) -> dict[str, str]:
    """Parse GFF9 attribute column into a dict."""
    attrs = {}
    for part in attr_string.split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs


def _parse_window_gff(gff_path: str):
    """Parse a shufflon-window GFF with embedded FASTA.

    Returns:
        features: list of dicts with keys: seqid, source, type, start, end,
                  score, strand, phase, attributes (raw string), attrs (dict)
        sequences: dict mapping seqid -> sequence string
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

    # Flush last FASTA sequence
    if current_seq_id is not None:
        sequences[current_seq_id] = "".join(current_seq_parts)

    return features, sequences


# ── Feature classification ───────────────────────────────────────────────────

def _classify_features(features: list[dict]) -> list[dict]:
    """Assign a category and display label to each feature.

    Categories: recombinase, inverted_repeat, invertible_segment,
                cds_with_ir, other_cds.
    """
    # Collect IR coordinate spans for overlap detection
    ir_spans = []
    for f in features:
        if f["type"] == "inverted_repeat":
            ir_spans.append((f["start"], f["end"]))

    # Collect hmm_hit spans so we can skip duplicate CDS at same coords
    hmm_spans = set()
    for f in features:
        if f["type"] == "hmm_hit":
            hmm_spans.add((f["start"], f["end"]))

    classified = []
    for f in features:
        ftype = f["type"]
        source = f["source"]
        attrs = f["attrs"]

        if ftype == "hmm_hit":
            f["category"] = "recombinase"
            locus = attrs.get("locus_tag", "")
            gene = attrs.get("gene", "")
            if gene:
                f["label"] = f"Recombinase ({gene})"
            elif locus:
                f["label"] = f"Recombinase ({locus})"
            else:
                f["label"] = "Recombinase"
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
            # Skip shufflon cassette CDS from shufflonfinder (legacy outputs)
            if source == "shufflonfinder":
                continue
            # Skip CDS that duplicates an hmm_hit at the same coordinates
            if (f["start"], f["end"]) in hmm_spans:
                continue

            # Check if this CDS overlaps any inverted repeat
            cds_start, cds_end = f["start"], f["end"]
            overlaps_ir = False
            for ir_s, ir_e in ir_spans:
                if cds_start <= ir_e and cds_end >= ir_s:
                    overlaps_ir = True
                    break

            if overlaps_ir:
                f["category"] = "cds_with_ir"
            else:
                f["category"] = "other_cds"

            product = attrs.get("product", "hypothetical protein")
            gene = attrs.get("gene", "")
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


# ── GenBank conversion ───────────────────────────────────────────────────────

def _gff_to_genbank(features: list[dict], sequences: dict, output_gbk: str):
    """Convert classified GFF features + sequence into a GenBank file.

    Each feature becomes a CDS in the GenBank (clinker only draws CDS
    features), with /label and /product qualifiers for labelling.

    Returns:
        gene_categories: dict mapping gene label -> category string
    """
    gene_categories = {}
    records = []

    for seqid, seq_str in sequences.items():
        record = SeqRecord(
            Seq(seq_str),
            id=seqid,
            name=seqid[:16],  # GenBank name limit
            description=f"Shufflon window {seqid}",
        )
        record.annotations["molecule_type"] = "DNA"

        seq_features = [f for f in features if f["seqid"] == seqid]

        for f in seq_features:
            # GFF is 1-based inclusive; Biopython FeatureLocation is 0-based half-open
            start = f["start"] - 1
            end = f["end"]

            if f["strand"] == "+":
                strand = 1
            elif f["strand"] == "-":
                strand = -1
            else:
                strand = 0

            location = FeatureLocation(start, end, strand=strand)

            label = f["label"]
            category = f["category"]

            # Clinker uses CDS features for drawing
            feature = SeqFeature(
                location=location,
                type="CDS",
                qualifiers={
                    "label": [label],
                    "product": [label],
                    "gene": [label],
                    "locus_tag": [f["attrs"].get("ID", label)],
                },
            )
            record.features.append(feature)
            gene_categories[label] = category

        records.append(record)

    with open(output_gbk, "w") as fh:
        SeqIO.write(records, fh, "genbank")

    return gene_categories


# ── Clinker CSV helpers ──────────────────────────────────────────────────────

def _write_gene_functions_csv(gene_categories: dict, output_csv: str):
    """Write a 2-column CSV: gene_label,function_category."""
    with open(output_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        for gene, cat in gene_categories.items():
            writer.writerow([gene, CATEGORY_LABELS.get(cat, cat)])


def _write_colour_map_csv(gene_categories: dict, output_csv: str):
    """Write a 2-column CSV: gene_label,hex_colour."""
    with open(output_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        for gene, cat in gene_categories.items():
            colour = CATEGORY_COLOURS.get(cat, COLOUR_OTHER_CDS)
            writer.writerow([gene, colour])


# ── Main entry point ─────────────────────────────────────────────────────────

def generate_clinker_plot(gff_path: str, output_dir: str) -> str | None:
    """Generate a Clinker plot for a single shufflon-window GFF.

    Args:
        gff_path: Path to a window GFF file (with embedded FASTA).
        output_dir: Directory to write plot and intermediate files.

    Returns:
        Path to the generated HTML plot, or None on failure.
    """
    ensure_dir(output_dir)
    basename = os.path.splitext(os.path.basename(gff_path))[0]

    # Parse GFF
    features, sequences = _parse_window_gff(gff_path)
    if not features or not sequences:
        logger.warning("No features or sequences in %s, skipping clinker", gff_path)
        return None

    # Classify features
    classified = _classify_features(features)
    if not classified:
        logger.warning("No classifiable features in %s, skipping clinker", gff_path)
        return None

    # Convert to GenBank
    gbk_path = os.path.join(output_dir, f"{basename}.gbk")
    gene_categories = _gff_to_genbank(classified, sequences, gbk_path)

    # Write gene_functions and colour_map CSVs
    func_csv = os.path.join(output_dir, f"{basename}_functions.csv")
    colour_csv = os.path.join(output_dir, f"{basename}_colours.csv")
    _write_gene_functions_csv(gene_categories, func_csv)
    _write_colour_map_csv(gene_categories, colour_csv)

    # Run clinker
    plot_path = os.path.join(output_dir, f"{basename}.html")
    cmd = [
        "clinker",
        gbk_path,
        "-p", plot_path,
        "-gf", func_csv,
        "-cm", colour_csv,
        "-na",  # no alignment (single cluster)
    ]

    logger.info("Running clinker: %s", " ".join(cmd))
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0:
            logger.error("clinker failed for %s: %s", gff_path, result.stderr)
            return None
        logger.info("Clinker plot saved: %s", plot_path)
        return plot_path
    except FileNotFoundError:
        logger.error("clinker not found on PATH. Install with: pip install clinker")
        return None
    except subprocess.TimeoutExpired:
        logger.error("clinker timed out for %s", gff_path)
        return None


def generate_clinker_plots(window_dir: str, clinker_dir: str) -> list[str]:
    """Generate Clinker plots for all window GFFs in a directory tree.

    Walks the window_dir looking for .gff files and writes all clinker
    outputs (GenBank, CSVs, HTML plots) into a mirrored subdirectory
    structure under clinker_dir.

    Args:
        window_dir: Top-level 07_shufflon_windows directory containing GFFs.
        clinker_dir: Top-level output directory for clinker results.

    Returns:
        List of paths to generated HTML plots.
    """
    plots = []
    for root, _dirs, files in os.walk(window_dir):
        for fname in sorted(files):
            if not fname.endswith(".gff"):
                continue
            gff_path = os.path.join(root, fname)
            # Mirror the subdirectory structure from window_dir into clinker_dir
            rel = os.path.relpath(root, window_dir)
            out_subdir = os.path.join(clinker_dir, rel)
            plot = generate_clinker_plot(gff_path, out_subdir)
            if plot:
                plots.append(plot)
    return plots
