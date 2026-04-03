#!/usr/bin/env python3
"""
shufflon-pipeline: Annotate shufflon structures in bacterial genomes.

Searches proteins against a library of shufflon-associated HMM profiles,
extracts flanking DNA (±5 kb by default) around each hit, detects inverted
repeats in those flanking regions (via PHAVA), and merges all annotations
into Prokka GFF files with windowed output.

Usage:
    # From raw genome FASTAs (runs Prokka first):
    python bin/run_pipeline.py --input-fasta genomes/ --outdir results/

    # From a single FASTA:
    python bin/run_pipeline.py --input-fasta genome.fna --outdir results/

    # From pre-annotated samples (skip Prokka):
    python bin/run_pipeline.py --sample-sheet samples.tsv --outdir results/

    # Custom HMM directory and flanking distance:
    python bin/run_pipeline.py --input-fasta genomes/ --hmm-dir my_hmms/ --flank-bp 10000 --outdir results/
"""

import argparse
import logging
import os
import sys

# Allow running from the repo root without installing
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, REPO_ROOT)

from lib.utils import setup_logging, ensure_dir, check_tool
from lib.sample_sheet import (
    Sample,
    load_sample_sheet,
    samples_from_fasta_dir,
    samples_from_single_fasta,
)
from lib.step_prokka import run_prokka
from lib.step_hmmsearch import (
    prepare_hmm_profiles,
    run_hmmsearch_all_profiles,
    combine_and_filter_hmmsearch,
)
from lib.step_flanking import (
    extract_flanking_regions,
    flanking_regions_to_tsv,
)
from lib.step_phava import run_phava_on_flanking, combine_ir_tables
from lib.step_gff import (
    hmm_hits_to_gff,
    ir_to_gff,
    merge_gff_into_prokka,
    extract_shufflon_windows,
)

logger = logging.getLogger("shufflon-pipeline")

# Default HMM directory: hmms/ at the repo root
DEFAULT_HMM_DIR = os.path.join(REPO_ROOT, "hmms")


def parse_args():
    p = argparse.ArgumentParser(
        description="Shufflon annotation pipeline for bacterial genomes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Input (mutually exclusive)
    inp = p.add_mutually_exclusive_group(required=True)
    inp.add_argument(
        "--input-fasta",
        help="Path to a genome FASTA file or a directory of FASTA files. "
             "Prokka will be run on these to produce .faa and .gff outputs.",
    )
    inp.add_argument(
        "--sample-sheet",
        help="TSV file with columns: sample_id, fna_path, faa_path, gff_path. "
             "Use this when Prokka has already been run externally.",
    )

    # Required
    p.add_argument(
        "--outdir",
        required=True,
        help="Output directory for all pipeline results.",
    )

    # Optional
    p.add_argument(
        "--hmm-dir",
        default=DEFAULT_HMM_DIR,
        help="Directory containing .hmm or .hmm.gz profile files "
             f"(default: hmms/ in the repo root).",
    )
    p.add_argument("--cpus", type=int, default=4, help="Threads per tool invocation (default: 4).")
    p.add_argument(
        "--bitscore", type=float, default=25.0,
        help="Minimum bitscore for HMM hits (default: 25.0).",
    )
    p.add_argument(
        "--flank-bp", type=int, default=5000,
        help="Base pairs of flanking DNA to extract on each side of an HMM "
             "hit protein for inverted repeat detection (default: 5000).",
    )
    p.add_argument(
        "--window-size", type=int, default=3000,
        help="Window size (bp) for clustering nearby IRs (default: 3000).",
    )
    p.add_argument(
        "--skip-prokka", action="store_true",
        help="Skip Prokka even for FASTA inputs (implies .faa/.gff exist alongside .fna).",
    )
    p.add_argument(
        "-v", "--verbose", action="count", default=1,
        help="Increase verbosity (-v = INFO, -vv = DEBUG).",
    )
    p.add_argument(
        "-q", "--quiet", action="store_true",
        help="Suppress all output except errors.",
    )

    return p.parse_args()


def resolve_samples(args) -> list[Sample]:
    """Build the sample list from CLI arguments."""
    if args.sample_sheet:
        return load_sample_sheet(args.sample_sheet)

    path = args.input_fasta
    if os.path.isdir(path):
        return samples_from_fasta_dir(path)
    elif os.path.isfile(path):
        return samples_from_single_fasta(path)
    else:
        logger.error("Input path does not exist: %s", path)
        sys.exit(1)


def main():
    args = parse_args()
    verbosity = 0 if args.quiet else args.verbose
    setup_logging(verbosity)

    outdir = ensure_dir(args.outdir)
    logger.info("Output directory: %s", outdir)

    # ---- Check external tools ----
    required_tools = ["hmmsearch", "hmmpress", "phava"]
    if not args.sample_sheet and not args.skip_prokka:
        required_tools.append("prokka")

    for tool in required_tools:
        check_tool(tool)

    # ---- Subdirectories ----
    dirs = {
        "prokka":       ensure_dir(os.path.join(outdir, "01_prokka")),
        "hmm":          ensure_dir(os.path.join(outdir, "02_hmmsearch")),
        "hmm_profiles":  ensure_dir(os.path.join(outdir, "02_hmmsearch", "profiles")),
        "hmm_results":  ensure_dir(os.path.join(outdir, "02_hmmsearch", "results")),
        "flanking":     ensure_dir(os.path.join(outdir, "03_flanking")),
        "phava":        ensure_dir(os.path.join(outdir, "04_phava")),
        "gff_hmm":      ensure_dir(os.path.join(outdir, "05_gff", "hmm_hits")),
        "gff_ir":       ensure_dir(os.path.join(outdir, "05_gff", "ir")),
        "gff_merged":   ensure_dir(os.path.join(outdir, "05_gff", "merged")),
        "windows":      ensure_dir(os.path.join(outdir, "06_shufflon_windows")),
    }

    # ==================================================================
    # Step 0: Resolve samples
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 0: Loading samples")
    logger.info("=" * 60)
    samples = resolve_samples(args)
    logger.info("Loaded %d sample(s)", len(samples))

    # ==================================================================
    # Step 1: Prokka (if needed)
    # ==================================================================
    prokka_needed = [s for s in samples if s.needs_prokka and not args.skip_prokka]
    if prokka_needed:
        logger.info("=" * 60)
        logger.info("STEP 1: Running Prokka on %d sample(s)", len(prokka_needed))
        logger.info("=" * 60)
        for sample in prokka_needed:
            run_prokka(sample, dirs["prokka"], cpus=args.cpus)
    else:
        logger.info("STEP 1: Prokka — skipped (all samples pre-annotated)")

    # Validate all samples have required paths
    for s in samples:
        s.validate()

    # ==================================================================
    # Step 2: Prepare HMM profiles + hmmsearch all profiles x all samples
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 2: HMM search (all profiles)")
    logger.info("=" * 60)

    profiles = prepare_hmm_profiles(args.hmm_dir, dirs["hmm_profiles"])
    logger.info("Searching %d HMM profiles against %d samples", len(profiles), len(samples))

    all_tblout_files = []
    for sample in samples:
        tblouts = run_hmmsearch_all_profiles(
            sample, profiles, dirs["hmm_results"], cpus=args.cpus
        )
        for t in tblouts:
            all_tblout_files.append((t, sample.sample_id))

    hits_combined_path = os.path.join(dirs["hmm"], "hmm_hits_combined.tsv")
    hmm_hits_df = combine_and_filter_hmmsearch(
        all_tblout_files, hits_combined_path, bitscore_threshold=args.bitscore
    )

    # ==================================================================
    # Step 3: Extract flanking DNA around each HMM hit
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 3: Extracting ±%d bp flanking DNA around HMM hits", args.flank_bp)
    logger.info("=" * 60)

    all_flanking_regions = []  # (sample, fasta_path, regions)
    all_flanking_flat = []

    for sample in samples:
        sample_hits = hmm_hits_df[hmm_hits_df["genome"] == sample.sample_id]
        fasta_path, regions = extract_flanking_regions(
            sample, sample_hits, dirs["flanking"], flank_bp=args.flank_bp
        )
        all_flanking_regions.append((sample, fasta_path, regions))
        all_flanking_flat.extend(regions)

    # Write combined flanking metadata
    if all_flanking_flat:
        flanking_tsv = os.path.join(dirs["flanking"], "flanking_regions_combined.tsv")
        flanking_regions_to_tsv(all_flanking_flat, flanking_tsv)

    # ==================================================================
    # Step 4: PHAVA inverted repeat detection on flanking regions
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 4: PHAVA inverted repeat detection (flanking regions)")
    logger.info("=" * 60)

    phava_results = []  # (phava_dir, sample_id, flanking_regions)
    for sample, fasta_path, regions in all_flanking_regions:
        if not regions:
            continue
        pdir = run_phava_on_flanking(sample, fasta_path, dirs["phava"], cpus=args.cpus)
        phava_results.append((pdir, sample.sample_id, regions))

    ir_combined_path = os.path.join(dirs["phava"], "IRs_combined_remapped.tsv")
    ir_df = combine_ir_tables(phava_results, ir_combined_path)

    # ==================================================================
    # Step 5: Generate and merge GFF files
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 5: GFF generation and merging")
    logger.info("=" * 60)

    # Generate HMM hit GFFs (per sample, with proper coordinates from GFF)
    hmm_gff_map = hmm_hits_to_gff(hmm_hits_df, samples, dirs["gff_hmm"])

    # Generate IR GFFs (per sample, genome-absolute coordinates)
    ir_gff_map = ir_to_gff(ir_df, dirs["gff_ir"])

    # Merge into Prokka GFFs
    for sample in samples:
        extra_gffs = []
        if sample.sample_id in hmm_gff_map:
            extra_gffs.append(hmm_gff_map[sample.sample_id])
        if sample.sample_id in ir_gff_map:
            extra_gffs.append(ir_gff_map[sample.sample_id])

        merged_path = os.path.join(dirs["gff_merged"], f"{sample.sample_id}.gff")
        merge_gff_into_prokka(sample.gff_path, extra_gffs, merged_path)

    # ==================================================================
    # Step 6: Extract shufflon windows
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 6: Extracting shufflon windows")
    logger.info("=" * 60)

    all_windows = []
    for sample in samples:
        merged_gff = os.path.join(dirs["gff_merged"], f"{sample.sample_id}.gff")
        if not os.path.isfile(merged_gff):
            logger.warning("No merged GFF for %s, skipping window extraction", sample.sample_id)
            continue

        windows = extract_shufflon_windows(
            merged_gff,
            os.path.join(dirs["windows"], sample.sample_id),
            window_size=args.window_size,
        )
        all_windows.extend(windows)

    # ==================================================================
    # Summary
    # ==================================================================
    logger.info("=" * 60)
    logger.info("PIPELINE COMPLETE")
    logger.info("=" * 60)
    logger.info("Samples processed:     %d", len(samples))
    logger.info("HMM profiles searched: %d", len(profiles))
    logger.info("HMM hits (filtered):   %d", len(hmm_hits_df))
    logger.info("Flanking regions:      %d", len(all_flanking_flat))
    logger.info("Inverted repeats:      %d", len(ir_df))
    logger.info("Shufflon windows:      %d", len(all_windows))
    logger.info("Results in:            %s", outdir)


if __name__ == "__main__":
    main()
