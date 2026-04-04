"""Step 3: Extract flanking DNA regions around HMM hit proteins for IR detection.

Maps protein HMM hits back to their genomic coordinates via the Prokka GFF,
then extracts up to flank_bp of DNA on each side of the CDS for targeted
inverted repeat scanning.
"""

import logging
import os
from collections import defaultdict
from dataclasses import dataclass

import pandas as pd

from .sample_sheet import Sample
from .utils import ensure_dir

logger = logging.getLogger("shufflonfinder")


@dataclass
class CdsCoords:
    """Genomic coordinates of a CDS feature from a GFF file."""
    contig: str
    start: int   # 1-based, as in GFF
    end: int     # 1-based, inclusive
    strand: str
    locus_tag: str


def parse_cds_from_gff(gff_path: str) -> dict[str, CdsCoords]:
    """Parse CDS features from a Prokka GFF into a locus_tag -> CdsCoords map.

    The locus_tag is extracted from the GFF attributes column (ID= field).

    Returns:
        Dict mapping locus_tag (protein ID) -> CdsCoords.
    """
    cds_map = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break  # Stop at FASTA section
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "CDS":
                continue

            contig = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = fields[8]

            # Extract locus_tag / ID from attributes
            locus_tag = None
            for attr in attrs.split(";"):
                if attr.startswith("ID="):
                    locus_tag = attr[3:]
                    break

            if locus_tag:
                cds_map[locus_tag] = CdsCoords(
                    contig=contig,
                    start=start,
                    end=end,
                    strand=strand,
                    locus_tag=locus_tag,
                )

    logger.debug("Parsed %d CDS features from %s", len(cds_map), gff_path)
    return cds_map


def parse_fasta_from_gff(gff_path: str) -> dict[str, str]:
    """Extract contig sequences from the ##FASTA section of a Prokka GFF.

    Returns:
        Dict mapping contig_id -> sequence string.
    """
    sequences = {}
    in_fasta = False
    current_id = None
    current_seq = []

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if not in_fasta:
                continue
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())

    if current_id:
        sequences[current_id] = "".join(current_seq)

    return sequences


def parse_fasta_file(fasta_path: str) -> dict[str, str]:
    """Parse a standalone FASTA file into a dict of id -> sequence.

    Returns:
        Dict mapping sequence_id -> sequence string.
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())

    if current_id:
        sequences[current_id] = "".join(current_seq)

    return sequences


@dataclass
class FlankingRegion:
    """A flanking DNA region extracted around one unique protein hit.

    When a protein matches multiple HMM profiles, only one flanking region
    is created. The hmm_profiles field records every profile that matched.
    """
    hit_id: str             # Unique identifier: <sample>__<locus_tag>
    sample_id: str
    hmm_profiles: list[str] # All HMM profiles that matched this protein
    locus_tag: str
    contig: str
    cds_start: int          # 1-based
    cds_end: int            # 1-based
    flank_start: int        # 1-based, clipped to contig bounds
    flank_end: int          # 1-based, clipped to contig bounds
    sequence: str           # The flanking DNA sequence


def _deduplicate_hits(hmm_hits_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse multi-profile hits to one row per protein.

    Groups by target_name (protein ID) and aggregates all matching
    hmm_profile values into a semicolon-separated string. Keeps the
    row with the highest bitscore as the representative.

    Returns:
        DataFrame with one row per unique target_name and a new
        column 'hmm_profiles' (semicolon-separated list).
    """
    if hmm_hits_df.empty:
        return hmm_hits_df

    # Collect all profiles per protein
    profiles_per_protein = (
        hmm_hits_df.groupby("target_name")["hmm_profile"]
        .apply(lambda x: ";".join(sorted(set(x))))
        .rename("hmm_profiles")
    )

    # Keep the best-scoring row per protein
    best_idx = hmm_hits_df.groupby("target_name")["full_sequence_bitscore"].idxmax()
    deduped = hmm_hits_df.loc[best_idx].copy()
    deduped = deduped.merge(profiles_per_protein, on="target_name", how="left")

    n_before = len(hmm_hits_df)
    n_after = len(deduped)
    if n_before != n_after:
        logger.info(
            "Deduplicated HMM hits: %d rows -> %d unique proteins (%d duplicate hits removed)",
            n_before, n_after, n_before - n_after,
        )

    return deduped


def extract_flanking_regions(
    sample: Sample,
    hmm_hits_df: pd.DataFrame,
    outdir: str,
    flank_bp: int = 5000,
) -> tuple[str, list[FlankingRegion]]:
    """Extract flanking DNA for each unique HMM-hit protein in a sample.

    Proteins that match multiple HMM profiles are deduplicated first so
    that each genomic locus produces exactly one flanking region. The
    FlankingRegion records which profiles matched.

    Args:
        sample: Sample with gff_path and fna_path set.
        hmm_hits_df: Filtered HMM hits DataFrame (subset for this sample).
                     Must have columns: target_name, hmm_profile.
        outdir: Directory to write the flanking FASTAs.
        flank_bp: Number of base pairs to extract on each side of the CDS.

    Returns:
        Tuple of (flanking_fasta_path, list of FlankingRegion objects).
        The FASTA contains one record per unique protein, named by hit_id.
    """
    ensure_dir(outdir)

    if hmm_hits_df.empty:
        logger.info("No HMM hits for %s, skipping flanking extraction", sample.sample_id)
        return "", []

    # Deduplicate: one flanking region per protein, regardless of how many
    # HMM profiles matched it
    deduped = _deduplicate_hits(hmm_hits_df)

    # Parse CDS coordinates from the GFF
    cds_map = parse_cds_from_gff(sample.gff_path)

    # Load contig sequences (try GFF first, fall back to standalone FNA)
    sequences = parse_fasta_from_gff(sample.gff_path)
    if not sequences:
        sequences = parse_fasta_file(sample.fna_path)

    if not sequences:
        logger.warning("No sequences found for %s", sample.sample_id)
        return "", []

    regions = []
    fasta_path = os.path.join(outdir, f"{sample.sample_id}_flanking.fasta")

    with open(fasta_path, "w") as fh:
        for _, hit in deduped.iterrows():
            protein_id = hit["target_name"]
            profiles_str = hit["hmm_profiles"]
            profiles_list = profiles_str.split(";")

            # Look up the CDS coordinates
            cds = cds_map.get(protein_id)
            if cds is None:
                logger.debug(
                    "Protein %s not found in GFF for %s, skipping",
                    protein_id, sample.sample_id,
                )
                continue

            contig_seq = sequences.get(cds.contig)
            if contig_seq is None:
                logger.debug(
                    "Contig %s not found in FASTA for %s, skipping",
                    cds.contig, sample.sample_id,
                )
                continue

            contig_len = len(contig_seq)

            # Calculate flanking region (1-based coords, clipped to contig bounds)
            flank_start = max(1, cds.start - flank_bp)
            flank_end = min(contig_len, cds.end + flank_bp)

            # Extract sequence (convert to 0-based for slicing)
            flank_seq = contig_seq[flank_start - 1 : flank_end]

            hit_id = f"{sample.sample_id}__{protein_id}"

            region = FlankingRegion(
                hit_id=hit_id,
                sample_id=sample.sample_id,
                hmm_profiles=profiles_list,
                locus_tag=protein_id,
                contig=cds.contig,
                cds_start=cds.start,
                cds_end=cds.end,
                flank_start=flank_start,
                flank_end=flank_end,
                sequence=flank_seq,
            )
            regions.append(region)

            # Write FASTA record
            fh.write(f">{hit_id} contig={cds.contig} "
                     f"flank={flank_start}-{flank_end} "
                     f"cds={cds.start}-{cds.end} "
                     f"hmm={profiles_str}\n")
            # Wrap sequence at 80 chars
            for i in range(0, len(flank_seq), 80):
                fh.write(flank_seq[i : i + 80] + "\n")

    logger.info(
        "Extracted %d flanking regions for %s (±%d bp) -> %s",
        len(regions), sample.sample_id, flank_bp, fasta_path,
    )
    return fasta_path, regions


def flanking_regions_to_tsv(regions: list[FlankingRegion], output_path: str) -> None:
    """Write flanking region metadata to a TSV for downstream reference.

    Args:
        regions: List of FlankingRegion objects.
        output_path: Where to write the TSV.
    """
    rows = []
    for r in regions:
        rows.append({
            "hit_id": r.hit_id,
            "sample_id": r.sample_id,
            "hmm_profiles": ";".join(r.hmm_profiles),
            "n_profiles": len(r.hmm_profiles),
            "locus_tag": r.locus_tag,
            "contig": r.contig,
            "cds_start": r.cds_start,
            "cds_end": r.cds_end,
            "flank_start": r.flank_start,
            "flank_end": r.flank_end,
            "flank_length": len(r.sequence),
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep="\t", index=False)
    logger.info("Wrote flanking region metadata: %d rows -> %s", len(df), output_path)
