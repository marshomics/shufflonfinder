"""Parse and validate sample input: either raw FASTAs or a pre-annotated sample sheet."""

import csv
import logging
import os
from dataclasses import dataclass, field

logger = logging.getLogger("shufflon-pipeline")


@dataclass
class Sample:
    """Represents one genome sample flowing through the pipeline."""
    sample_id: str
    fna_path: str = ""        # nucleotide FASTA
    faa_path: str = ""        # protein FASTA (Prokka output)
    gff_path: str = ""        # GFF3 annotation (Prokka output)
    needs_prokka: bool = True  # True if Prokka has not been run yet

    def validate(self) -> None:
        """Check that all required file paths exist."""
        for attr in ("fna_path", "faa_path", "gff_path"):
            path = getattr(self, attr)
            if path and not os.path.isfile(path):
                raise FileNotFoundError(
                    f"Sample '{self.sample_id}': {attr} not found: {path}"
                )


def load_sample_sheet(tsv_path: str) -> list[Sample]:
    """Load a TSV sample sheet with columns: sample_id, fna_path, faa_path, gff_path.

    All paths in the sheet must point to existing files.

    Returns:
        List of Sample objects with needs_prokka=False.
    """
    samples = []
    required_cols = {"sample_id", "fna_path", "faa_path", "gff_path"}

    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not required_cols.issubset(set(reader.fieldnames or [])):
            missing = required_cols - set(reader.fieldnames or [])
            raise ValueError(
                f"Sample sheet is missing required columns: {missing}"
            )

        for i, row in enumerate(reader, start=2):
            sample = Sample(
                sample_id=row["sample_id"].strip(),
                fna_path=os.path.abspath(row["fna_path"].strip()),
                faa_path=os.path.abspath(row["faa_path"].strip()),
                gff_path=os.path.abspath(row["gff_path"].strip()),
                needs_prokka=False,
            )
            sample.validate()
            samples.append(sample)

    logger.info("Loaded %d samples from %s", len(samples), tsv_path)
    return samples


def samples_from_fasta_dir(fasta_dir: str) -> list[Sample]:
    """Create Sample objects from a directory of genome FASTA files.

    Recognized extensions: .fasta, .fa, .fna, .fasta.gz, .fa.gz, .fna.gz

    Returns:
        List of Sample objects with needs_prokka=True.
    """
    extensions = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")
    samples = []

    for fname in sorted(os.listdir(fasta_dir)):
        if not any(fname.endswith(ext) for ext in extensions):
            continue
        # Strip all recognized extensions to get sample_id
        sample_id = fname
        for ext in sorted(extensions, key=len, reverse=True):
            if sample_id.endswith(ext):
                sample_id = sample_id[: -len(ext)]
                break
        samples.append(
            Sample(
                sample_id=sample_id,
                fna_path=os.path.abspath(os.path.join(fasta_dir, fname)),
                needs_prokka=True,
            )
        )

    logger.info("Found %d FASTA files in %s", len(samples), fasta_dir)
    return samples


def samples_from_single_fasta(fasta_path: str) -> list[Sample]:
    """Create a single Sample from one FASTA file.

    Returns:
        List containing one Sample with needs_prokka=True.
    """
    fname = os.path.basename(fasta_path)
    extensions = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")
    sample_id = fname
    for ext in sorted(extensions, key=len, reverse=True):
        if sample_id.endswith(ext):
            sample_id = sample_id[: -len(ext)]
            break

    sample = Sample(
        sample_id=sample_id,
        fna_path=os.path.abspath(fasta_path),
        needs_prokka=True,
    )
    logger.info("Single FASTA input: %s", sample.sample_id)
    return [sample]
