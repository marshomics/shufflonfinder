"""Step 1: Run Prokka gene annotation on raw genome FASTAs."""

import logging
import os

from .sample_sheet import Sample
from .utils import ensure_dir, run_cmd

logger = logging.getLogger("shufflonfinder")


def run_prokka(sample: Sample, outdir: str, cpus: int = 4) -> Sample:
    """Annotate a genome FASTA with Prokka and update the Sample paths.

    Args:
        sample: A Sample whose fna_path is set.
        outdir: Directory for Prokka output files.
        cpus: Number of threads for Prokka.

    Returns:
        The same Sample object with faa_path and gff_path populated.
    """
    ensure_dir(outdir)

    cmd = [
        "prokka",
        "--outdir", outdir,
        "--prefix", sample.sample_id,
        "--cpus", str(cpus),
        "--force",
        "--compliant",
        sample.fna_path,
    ]
    run_cmd(cmd, description=f"Prokka: {sample.sample_id}")

    sample.faa_path = os.path.join(outdir, f"{sample.sample_id}.faa")
    sample.gff_path = os.path.join(outdir, f"{sample.sample_id}.gff")
    sample.fna_path = os.path.join(outdir, f"{sample.sample_id}.fna")
    sample.needs_prokka = False

    # Verify outputs exist
    for attr in ("faa_path", "gff_path"):
        path = getattr(sample, attr)
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"Prokka did not produce expected output: {path}"
            )

    logger.info("Prokka complete for %s", sample.sample_id)
    return sample
