"""Step 2: Search protein FASTAs against all HMM profiles in the hmms/ directory."""

import glob
import gzip
import logging
import os
import shutil

import pandas as pd

from .sample_sheet import Sample
from .utils import ensure_dir, run_cmd

logger = logging.getLogger("shufflon-pipeline")

HMMSEARCH_COLUMNS = [
    "target_name", "accession1", "query_name", "accession2",
    "full_sequence_evalue", "full_sequence_bitscore", "full_sequence_bias",
    "best_1_domain_evalue", "best_1_domain_bitscore", "best_1_domain_bias",
    "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc",
    "description",
]


def prepare_hmm_profiles(hmm_dir: str, work_dir: str) -> list[str]:
    """Decompress and prepare all HMM profiles from the hmms/ directory.

    Handles both .hmm and .hmm.gz files. Each profile is decompressed (if
    needed) into work_dir and converted/pressed for hmmsearch.

    Args:
        hmm_dir: Directory containing .hmm or .hmm.gz files.
        work_dir: Working directory for decompressed/converted profiles.

    Returns:
        List of paths to search-ready HMM profile files.
    """
    ensure_dir(work_dir)

    raw_files = (
        glob.glob(os.path.join(hmm_dir, "*.hmm.gz"))
        + glob.glob(os.path.join(hmm_dir, "*.hmm"))
    )
    # Deduplicate (if both .hmm and .hmm.gz exist for the same profile)
    seen_names = set()
    unique_files = []
    for f in raw_files:
        name = os.path.basename(f).replace(".gz", "").replace(".hmm", "")
        if name not in seen_names:
            seen_names.add(name)
            unique_files.append(f)

    if not unique_files:
        raise FileNotFoundError(f"No .hmm or .hmm.gz files found in {hmm_dir}")

    logger.info("Found %d HMM profiles in %s", len(unique_files), hmm_dir)

    ready_profiles = []
    for src in sorted(unique_files):
        basename = os.path.basename(src).replace(".gz", "")
        profile_name = basename.replace(".hmm", "")
        dest = os.path.join(work_dir, basename)

        # Decompress if gzipped
        if src.endswith(".gz"):
            if not os.path.isfile(dest):
                with gzip.open(src, "rb") as gz_in, open(dest, "wb") as f_out:
                    shutil.copyfileobj(gz_in, f_out)
                logger.debug("Decompressed %s -> %s", src, dest)
        else:
            if not os.path.isfile(dest):
                shutil.copy(src, dest)

        # hmmpress for indexed search (skip if already pressed)
        h3i_file = dest + ".h3i"
        if not os.path.isfile(h3i_file):
            try:
                run_cmd(["hmmpress", "-f", dest], description=f"hmmpress {profile_name}")
            except Exception as e:
                logger.warning("hmmpress failed for %s: %s — will try hmmsearch anyway", profile_name, e)

        ready_profiles.append(dest)

    logger.info("Prepared %d HMM profiles for searching", len(ready_profiles))
    return ready_profiles


def run_hmmsearch(
    sample: Sample,
    hmm_profile: str,
    outdir: str,
    cpus: int = 1,
) -> str:
    """Run hmmsearch on a sample's protein FASTA against one HMM profile.

    Args:
        sample: Sample with faa_path set.
        hmm_profile: Path to a single HMM profile file.
        outdir: Directory to write results.
        cpus: Number of threads.

    Returns:
        Path to the tblout results file.
    """
    ensure_dir(outdir)
    profile_name = os.path.basename(hmm_profile).replace(".hmm", "")
    tblout = os.path.join(outdir, f"{sample.sample_id}__{profile_name}.tblout")

    cmd = [
        "hmmsearch",
        "--cpu", str(cpus),
        "--tblout", tblout,
        hmm_profile,
        sample.faa_path,
    ]
    run_cmd(cmd, description=f"hmmsearch {profile_name}: {sample.sample_id}")
    return tblout


def run_hmmsearch_all_profiles(
    sample: Sample,
    profiles: list[str],
    outdir: str,
    cpus: int = 1,
) -> list[str]:
    """Run hmmsearch for a sample against every HMM profile.

    Args:
        sample: Sample with faa_path set.
        profiles: List of HMM profile paths.
        outdir: Directory to write results.
        cpus: Number of threads per search.

    Returns:
        List of tblout file paths (one per profile).
    """
    sample_dir = ensure_dir(os.path.join(outdir, sample.sample_id))
    tblouts = []
    for profile in profiles:
        tblout = run_hmmsearch(sample, profile, sample_dir, cpus=cpus)
        tblouts.append(tblout)
    return tblouts


def parse_hmmsearch_tblout(tblout_path: str, sample_id: str) -> pd.DataFrame:
    """Parse an hmmsearch --tblout file into a DataFrame.

    Args:
        tblout_path: Path to the tblout output.
        sample_id: Sample identifier to add as a column.

    Returns:
        DataFrame with parsed hits (may be empty). Includes an 'hmm_profile'
        column derived from the filename.
    """
    # Extract profile name from filename (e.g., "sample__PF04917.tblout" -> "PF04917")
    basename = os.path.basename(tblout_path)
    if "__" in basename:
        profile_name = basename.split("__", 1)[1].replace(".tblout", "")
    else:
        profile_name = basename.replace(".tblout", "")

    rows = []
    with open(tblout_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split(None, 18)
            if len(parts) >= 19:
                rows.append(parts)

    if not rows:
        return pd.DataFrame(columns=HMMSEARCH_COLUMNS + ["genome", "hmm_profile"])

    df = pd.DataFrame(rows, columns=HMMSEARCH_COLUMNS)
    df["genome"] = sample_id
    df["hmm_profile"] = profile_name
    df["full_sequence_bitscore"] = pd.to_numeric(df["full_sequence_bitscore"], errors="coerce")
    df["full_sequence_evalue"] = pd.to_numeric(df["full_sequence_evalue"], errors="coerce")
    return df


def combine_and_filter_hmmsearch(
    tblout_files: list[tuple[str, str]],
    output_path: str,
    bitscore_threshold: float = 25.0,
) -> pd.DataFrame:
    """Combine hmmsearch results across all samples and profiles, filter by bitscore.

    Args:
        tblout_files: List of (tblout_path, sample_id) tuples.
        output_path: Where to write the combined, filtered TSV.
        bitscore_threshold: Minimum full_sequence_bitscore to keep.

    Returns:
        Filtered DataFrame with columns including 'hmm_profile' and 'genome'.
    """
    frames = []
    for path, sid in tblout_files:
        df = parse_hmmsearch_tblout(path, sid)
        if not df.empty:
            frames.append(df)

    if not frames:
        logger.warning("No HMM hits found across any samples or profiles.")
        combined = pd.DataFrame(columns=HMMSEARCH_COLUMNS + ["genome", "hmm_profile"])
    else:
        combined = pd.concat(frames, ignore_index=True)

    filtered = combined[combined["full_sequence_bitscore"] >= bitscore_threshold].copy()
    filtered.to_csv(output_path, sep="\t", index=False)

    n_profiles = filtered["hmm_profile"].nunique() if not filtered.empty else 0
    logger.info(
        "HMM hits: %d total, %d pass bitscore >= %.1f (across %d profiles)",
        len(combined), len(filtered), bitscore_threshold, n_profiles,
    )
    return filtered
