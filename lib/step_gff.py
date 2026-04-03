"""Step 5: Generate, merge, and filter GFF files for shufflon annotations."""

import csv
import logging
import os
from collections import defaultdict
from io import StringIO

import pandas as pd

from .sample_sheet import Sample
from .step_flanking import FlankingRegion, parse_cds_from_gff
from .utils import ensure_dir

logger = logging.getLogger("shufflon-pipeline")


# ---------------------------------------------------------------------------
# GFF generation from HMM hits (generic, multi-profile)
# ---------------------------------------------------------------------------

def hmm_hits_to_gff(
    hmm_hits_df: pd.DataFrame,
    samples: list[Sample],
    output_dir: str,
) -> dict[str, str]:
    """Convert HMM hit proteins into per-sample GFF annotation files.

    Each hit becomes a CDS feature with source set to the HMM profile name.
    Coordinates are looked up from the sample's Prokka GFF.

    Args:
        hmm_hits_df: Filtered HMM hits DataFrame with columns:
                     target_name, genome, hmm_profile.
        samples: List of all Sample objects (for GFF coordinate lookup).
        output_dir: Directory to write per-sample .gff files.

    Returns:
        Dict mapping sample_id -> path to the HMM hits GFF file.
    """
    ensure_dir(output_dir)
    files = {}

    if hmm_hits_df.empty:
        return files

    sample_map = {s.sample_id: s for s in samples}

    grouped = hmm_hits_df.groupby("genome")
    for sample_id, group in grouped:
        sample = sample_map.get(sample_id)
        if sample is None:
            logger.warning("Sample %s not found, skipping GFF generation", sample_id)
            continue

        # Parse CDS coordinates from this sample's GFF
        cds_map = parse_cds_from_gff(sample.gff_path)

        gff_path = os.path.join(output_dir, f"{sample_id}_hmm_hits.gff")
        written = 0
        with open(gff_path, "w") as fh:
            for _, row in group.iterrows():
                protein_id = row["target_name"]
                hmm_profile = row["hmm_profile"]

                cds = cds_map.get(protein_id)
                if cds is None:
                    continue

                attrs = (
                    f"ID={cds.locus_tag};"
                    f"Parent={cds.locus_tag}_gene;"
                    f"Name={hmm_profile}"
                )
                fields = [
                    cds.contig,
                    hmm_profile,
                    "CDS",
                    str(cds.start),
                    str(cds.end),
                    ".",
                    cds.strand,
                    ".",
                    attrs,
                ]
                fh.write("\t".join(fields) + "\n")
                written += 1

        files[sample_id] = gff_path
        logger.debug("Wrote HMM hits GFF for %s: %d features", sample_id, written)

    return files


# ---------------------------------------------------------------------------
# GFF generation from IR coordinates
# ---------------------------------------------------------------------------

def ir_to_gff(ir_df: pd.DataFrame, output_dir: str) -> dict[str, str]:
    """Convert inverted repeat table rows into per-sample GFF annotation files.

    Each IR produces two features (left/forward and right/reverse) with
    source='einverted' and type='inverted_repeat'.

    Args:
        ir_df: DataFrame with columns: IR_Chr, LeftIRStart, LeftIRStop,
               RightIRStart, RightIRStop, sample_id.
        output_dir: Directory to write per-sample .gff files.

    Returns:
        Dict mapping sample_id -> path to the IR GFF file.
    """
    ensure_dir(output_dir)
    files = {}

    if ir_df.empty:
        return files

    grouped = ir_df.groupby("sample_id")
    for sample_id, group in grouped:
        gff_path = os.path.join(output_dir, f"{sample_id}_ir.gff")
        with open(gff_path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            for counter, (_, row) in enumerate(group.iterrows(), start=1):
                seq_num = f"{counter:04d}"
                # Left (forward) IR
                writer.writerow([
                    row["IR_Chr"], "einverted", "inverted_repeat",
                    row["LeftIRStart"], row["LeftIRStop"],
                    ".", "+", ".",
                    f"ID=inverted_repeat_{seq_num}_F;Name=inverted_repeat_{seq_num}_F",
                ])
                # Right (reverse) IR
                writer.writerow([
                    row["IR_Chr"], "einverted", "inverted_repeat",
                    row["RightIRStart"], row["RightIRStop"],
                    ".", "-", ".",
                    f"ID=inverted_repeat_{seq_num}_R;Name=inverted_repeat_{seq_num}_R",
                ])
        files[sample_id] = gff_path
        logger.debug("Wrote IR GFF for %s: %d IR pairs", sample_id, len(group))

    return files


# ---------------------------------------------------------------------------
# Merge annotations into Prokka GFF
# ---------------------------------------------------------------------------

def merge_gff_into_prokka(
    prokka_gff: str,
    extra_gff_files: list[str],
    output_gff: str,
) -> str:
    """Insert additional GFF features into a Prokka GFF file.

    Prokka GFFs have a ##FASTA section at the bottom. This function inserts
    the extra annotations between the Prokka feature lines and the FASTA block.

    Args:
        prokka_gff: Path to the original Prokka .gff file.
        extra_gff_files: List of paths to GFF files to merge in.
        output_gff: Where to write the merged GFF.

    Returns:
        Path to the merged GFF file.
    """
    ensure_dir(os.path.dirname(output_gff))

    annotation_lines = []
    fasta_lines = []
    in_fasta = False

    with open(prokka_gff) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                in_fasta = True
                fasta_lines.append(line)
            elif in_fasta:
                fasta_lines.append(line)
            else:
                annotation_lines.append(line)

    # Read extra features
    extra_lines = []
    for gff_path in extra_gff_files:
        if not os.path.isfile(gff_path):
            continue
        with open(gff_path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    extra_lines.append(line)

    with open(output_gff, "w") as fh:
        fh.writelines(annotation_lines)
        if extra_lines:
            # Ensure newline separation
            if annotation_lines and not annotation_lines[-1].endswith("\n"):
                fh.write("\n")
            fh.writelines(extra_lines)
        fh.writelines(fasta_lines)

    logger.debug("Merged GFF -> %s (%d extra features)", output_gff, len(extra_lines))
    return output_gff


# ---------------------------------------------------------------------------
# Window-based GFF filtering (IR + CDS co-location)
# ---------------------------------------------------------------------------

class Feature:
    """Lightweight representation of a GFF feature for spatial operations."""

    __slots__ = ("seqid", "source", "type", "start", "end", "strand", "attributes")

    def __init__(self, seqid, source, type, start, end, strand, attributes):
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.strand = strand
        self.attributes = attributes

    def overlaps(self, other):
        return self.seqid == other.seqid and self.start <= other.end and self.end >= other.start

    def distance_to(self, other):
        if self.seqid != other.seqid:
            return float("inf")
        if self.overlaps(other):
            return 0
        return min(abs(self.start - other.end), abs(self.end - other.start))


def parse_gff_with_fasta(gff_path: str):
    """Parse a GFF file (with embedded FASTA) into features and sequences.

    Returns:
        Tuple of (ir_by_contig, cds_by_contig, sequences)
        where sequences is a dict of contig_id -> (header, sequence_str).
    """
    ir_by_contig = defaultdict(list)
    cds_by_contig = defaultdict(list)
    fasta_lines = []
    in_fasta = False

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if in_fasta:
                fasta_lines.append(line)
                continue
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feat = Feature(
                seqid=fields[0],
                source=fields[1],
                type=fields[2],
                start=int(fields[3]) - 1,  # convert to 0-based
                end=int(fields[4]),
                strand=fields[6],
                attributes=fields[8],
            )

            if feat.type == "inverted_repeat":
                ir_by_contig[feat.seqid].append(feat)
            elif feat.type == "CDS":
                cds_by_contig[feat.seqid].append(feat)

    # Parse FASTA
    sequences = {}
    current_id = None
    current_seq = []
    for line in fasta_lines:
        if line.startswith(">"):
            if current_id:
                sequences[current_id] = "".join(current_seq)
            current_id = line[1:].strip().split()[0]
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_id:
        sequences[current_id] = "".join(current_seq)

    return ir_by_contig, cds_by_contig, sequences


def group_features_by_window(features: list[Feature], window_size: int = 3000):
    """Group features that are within window_size bp of each other."""
    if not features:
        return []

    sorted_feats = sorted(features, key=lambda f: f.start)
    groups = [[sorted_feats[0]]]

    for feat in sorted_feats[1:]:
        in_window = any(
            feat.distance_to(g) <= window_size for g in groups[-1]
        )
        if in_window:
            groups[-1].append(feat)
        else:
            groups.append([feat])

    return groups


def find_intersecting_cds(ir_group, cds_features):
    """Find CDS features overlapping the span of an IR group."""
    group_start = min(ir.start for ir in ir_group)
    group_end = max(ir.end for ir in ir_group)
    seqid = ir_group[0].seqid
    window = Feature(seqid, ".", "window", group_start, group_end, "+", "")
    return [cds for cds in cds_features if window.overlaps(cds)]


def extract_shufflon_windows(
    merged_gff: str,
    output_dir: str,
    window_size: int = 3000,
) -> list[str]:
    """Extract GFF+FASTA windows where inverted repeats co-locate with CDS features.

    This is the final filtering step: for each cluster of nearby IRs, find
    overlapping CDS features and write a self-contained GFF with the windowed
    sequence.

    Args:
        merged_gff: Path to a merged GFF (Prokka + pilV + IR annotations, with FASTA).
        output_dir: Where to write the per-window GFF files.
        window_size: Maximum distance (bp) to cluster IRs together.

    Returns:
        List of paths to the output window GFF files.
    """
    ensure_dir(output_dir)
    input_prefix = os.path.splitext(os.path.basename(merged_gff))[0]

    ir_by_contig, cds_by_contig, sequences = parse_gff_with_fasta(merged_gff)

    if not sequences:
        logger.warning("No FASTA sequences found in %s", merged_gff)
        return []

    output_files = []
    window_counter = defaultdict(int)

    for contig_id, ir_features in ir_by_contig.items():
        ir_groups = group_features_by_window(ir_features, window_size)
        contig_cds = cds_by_contig.get(contig_id, [])
        seq = sequences.get(contig_id, "")
        seq_len = len(seq)

        for ir_group in ir_groups:
            intersecting = find_intersecting_cds(ir_group, contig_cds)
            if not intersecting:
                continue

            window_counter[contig_id] += 1
            wnum = window_counter[contig_id]

            # Calculate window bounds
            ir_start = min(ir.start for ir in ir_group)
            ir_end = max(ir.end for ir in ir_group)
            cds_start = min(c.start for c in intersecting)
            cds_end = max(c.end for c in intersecting)
            win_start = max(0, min(ir_start, cds_start))
            win_end = min(seq_len, max(ir_end, cds_end))

            if win_start >= win_end:
                continue

            seq_id = f"{input_prefix}_{contig_id}_{wnum}"
            out_path = os.path.join(
                output_dir,
                f"{input_prefix}_contig_{contig_id}_window_{wnum}.gff",
            )

            with open(out_path, "w") as fh:
                fh.write("##gff-version 3\n")
                for ir in ir_group:
                    fh.write(
                        f"{seq_id}\t.\tinverted_repeat\t"
                        f"{ir.start - win_start + 1}\t{ir.end - win_start}\t.\t"
                        f"{ir.strand}\t.\t{ir.attributes}\n"
                    )
                for cds in intersecting:
                    fh.write(
                        f"{seq_id}\t.\tCDS\t"
                        f"{cds.start - win_start + 1}\t{cds.end - win_start}\t.\t"
                        f"{cds.strand}\t.\t{cds.attributes}\n"
                    )
                fh.write("##FASTA\n")
                fh.write(f">{seq_id}\n{seq[win_start:win_end]}\n")

            output_files.append(out_path)

    logger.info(
        "Extracted %d shufflon windows from %s", len(output_files), merged_gff
    )
    return output_files
