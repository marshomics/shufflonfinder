#!/usr/bin/env python3
"""
shufflonfinder: Annotate shufflon structures in bacterial genomes.

Searches proteins against a library of shufflon-associated HMM profiles,
extracts flanking DNA (±5 kb by default) around each hit, detects inverted
repeats in those flanking regions (via EMBOSS einverted), and merges all
annotations into Prokka GFF files with windowed output.

Usage:
    # From raw genome FASTAs (runs Prokka first):
    shufflonfinder --input-fasta genomes/ --outdir results/

    # From a single FASTA:
    shufflonfinder --input-fasta genome.fna --outdir results/

    # From pre-annotated samples (skip Prokka):
    shufflonfinder --sample-sheet samples.tsv --outdir results/

    # Custom HMM directory and flanking distance:
    shufflonfinder --input-fasta genomes/ --hmm-dir my_hmms/ --flank-bp 10000 --outdir results/
"""

import argparse
import importlib.resources
import logging
import os
import shutil
import sys
import tempfile
from collections import defaultdict

import pandas as pd

from .utils import setup_logging, ensure_dir, check_tool
from .sample_sheet import (
    Sample,
    load_sample_sheet,
    samples_from_fasta_dir,
    samples_from_single_fasta,
)
from .step_prokka import run_prokka
from .step_hmmsearch import (
    prepare_hmm_profiles,
    run_hmmsearch_all_profiles,
    combine_and_filter_hmmsearch,
)
from .step_flanking import (
    extract_flanking_regions,
    flanking_regions_to_tsv,
    parse_fasta_from_gff,
)
from .step_phava import (
    detect_inverted_repeats,
    combine_ir_tables,
    filter_ir_table,
    filter_irs_in_hmm_hits,
    filter_shufflon_candidates,
    refine_sfx_sites,
)
from .step_gff import (
    hmm_hits_to_gff,
    ir_to_gff,
    merge_gff_into_prokka,
    extract_shufflon_windows,
    shufflon_windows_to_tsv,
)
from .step_clinker import generate_shufflon_plots
from .step_kofamscan import run_kofamscan, load_ko_table
from .step_kegg_enrichment import run_kegg_enrichment_for_sample

logger = logging.getLogger("shufflonfinder")


def _bundled_hmm_dir() -> str:
    """Return the path to the bundled hmms/ directory inside the package."""
    try:
        # Python 3.9+
        ref = importlib.resources.files("shufflonfinder") / "hmms"
        return str(ref)
    except AttributeError:
        # Fallback for older Python
        return os.path.join(os.path.dirname(__file__), "hmms")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="shufflonfinder",
        description="Annotate shufflon structures in bacterial genomes.",
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
        default=None,
        help="Directory containing .hmm or .hmm.gz profile files "
             "(default: bundled profiles shipped with shufflonfinder).",
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
        "--min-ir-arm-length", type=int, default=13,
        help="Minimum IR arm length in bp to keep. Shufflon sfx "
             "recognition sites are typically 19 bp, but einverted may "
             "report shorter partial alignments for genuine sites "
             "(default: 13).",
    )
    p.add_argument(
        "--max-ir-arm-length", type=int, default=35,
        help="Maximum IR arm length in bp to keep. Shufflon sfx recognition "
             "sites are typically 19 bp; longer arms usually come from "
             "transposon or IS-element IRs, or from einverted extending "
             "across two abutting sites (default: 35).",
    )
    p.add_argument(
        "--min-ir-identity", type=float, default=85.0,
        help="Minimum percent identity between IR arms to keep. "
             "Within-type sfx pairs are near-perfect reverse "
             "complements (≥95%%), but cross-type pairs and "
             "divergent recognition sites may be lower; 85%% "
             "retains genuine sites while excluding most noise "
             "(default: 85.0).",
    )
    p.add_argument(
        "--min-ir-pairs", type=int, default=3,
        help="Minimum number of IR pairs per cluster to qualify as a "
             "shufflon candidate (default: 3).",
    )
    p.add_argument(
        "--cluster-distance", type=int, default=1000,
        help="Maximum gap (bp) between adjacent IR pairs for chaining "
             "into one cluster during shufflon candidate filtering. "
             "Tighter than --window-size to isolate the dense shufflon "
             "core from neighboring transposon or IS-element repeats "
             "(default: 1000).",
    )
    p.add_argument(
        "--min-ir-density", type=float, default=1.0,
        help="Minimum IR pairs per kilobase within a cluster to qualify "
             "as a shufflon candidate. einverted typically detects 3-5 "
             "of the pairwise sfx combinations in a shufflon spanning "
             "~2 kb, giving densities of 1.5-2.5 pairs/kb (default: 1.0).",
    )
    p.add_argument(
        "--skip-prokka", action="store_true",
        help="Skip Prokka even for FASTA inputs (implies .faa/.gff exist alongside .fna).",
    )
    p.add_argument(
        "--ko-profiles-dir",
        default=None,
        help="Path to the KOfamscan profile directory (must contain "
             "profiles/ and ko_list). When provided, KOfamscan is run "
             "on each sample's .faa to assign KO accessions, and KEGG "
             "enrichment analysis is performed on inverton-like windows.",
    )
    p.add_argument(
        "--skip-kofamscan", action="store_true",
        help="Skip KOfamscan even when --ko-profiles-dir is set. Use this "
             "when the sample sheet already has a kofamscan_path column.",
    )
    p.add_argument(
        "-v", "--verbose", action="count", default=1,
        help="Increase verbosity (-v = INFO, -vv = DEBUG).",
    )
    p.add_argument(
        "-q", "--quiet", action="store_true",
        help="Suppress all output except errors.",
    )

    return p.parse_args(argv)


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


def main(argv=None):
    args = parse_args(argv)
    verbosity = 0 if args.quiet else args.verbose
    setup_logging(verbosity)

    # Resolve HMM directory: user-supplied or bundled
    hmm_dir = args.hmm_dir if args.hmm_dir else _bundled_hmm_dir()

    outdir = ensure_dir(args.outdir)
    logger.info("Output directory: %s", outdir)

    # ---- Check external tools ----
    required_tools = ["hmmsearch", "hmmpress", "einverted"]
    if not args.sample_sheet and not args.skip_prokka:
        required_tools.append("prokka")

    # KOfamscan is needed when ko_profiles_dir is set and we're not skipping
    run_kofamscan_step = bool(args.ko_profiles_dir) and not args.skip_kofamscan
    if run_kofamscan_step:
        required_tools.append("exec_annotation")

    for tool in required_tools:
        check_tool(tool)

    # ==================================================================
    # Step 0: Resolve samples
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 0: Loading samples")
    logger.info("=" * 60)
    samples = resolve_samples(args)
    logger.info("Loaded %d sample(s)", len(samples))

    # ---- Per-sample directory trees ----
    # Structure: outdir/<sample_id>/01_prokka/, 02_hmmsearch/, …, 07_shufflon_windows/
    # Shared resources (HMM profiles) stay at outdir level.
    def sample_dirs(sample_id: str) -> dict[str, str]:
        root = ensure_dir(os.path.join(outdir, sample_id))
        win_root = os.path.join(root, "07_shufflon_windows")
        return {
            "root":         root,
            "prokka":       ensure_dir(os.path.join(root, "01_prokka")),
            "hmm":          ensure_dir(os.path.join(root, "02_hmmsearch")),
            "hmm_results":  ensure_dir(os.path.join(root, "02_hmmsearch", "results")),
            "flanking":     ensure_dir(os.path.join(root, "03_flanking")),
            "ir":           ensure_dir(os.path.join(root, "04_inverted_repeats")),
            "shufflon":     ensure_dir(os.path.join(root, "05_shufflon_filter")),
            "gff_hmm":      ensure_dir(os.path.join(root, "06_gff", "hmm_hits")),
            "gff_ir":       ensure_dir(os.path.join(root, "06_gff", "ir")),
            "gff_merged":   ensure_dir(os.path.join(root, "06_gff", "merged")),
            "windows":      ensure_dir(win_root),
            "shufflon_gffs":  ensure_dir(os.path.join(win_root, "shufflon_like", "gffs")),
            "shufflon_plots": ensure_dir(os.path.join(win_root, "shufflon_like", "plots")),
            "inverton_gffs":  ensure_dir(os.path.join(win_root, "inverton_like", "gffs")),
            "inverton_plots": ensure_dir(os.path.join(win_root, "inverton_like", "plots")),
            "kofamscan":    ensure_dir(os.path.join(root, "08_kofamscan")),
            "kegg_enrich":  ensure_dir(os.path.join(root, "09_kegg_enrichment")),
        }

    # Build dirs for each sample up front
    sdirs = {s.sample_id: sample_dirs(s.sample_id) for s in samples}

    # ==================================================================
    # Step 1: Prokka (if needed)
    # ==================================================================
    prokka_needed = [s for s in samples if s.needs_prokka and not args.skip_prokka]
    if prokka_needed:
        logger.info("=" * 60)
        logger.info("STEP 1: Running Prokka on %d sample(s)", len(prokka_needed))
        logger.info("=" * 60)
        for sample in prokka_needed:
            run_prokka(sample, sdirs[sample.sample_id]["prokka"], cpus=args.cpus)
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

    # Prepared profiles go into a temp directory (cleaned up on exit)
    _hmm_tmpdir = tempfile.mkdtemp(prefix="shufflonfinder_hmm_")
    profiles = prepare_hmm_profiles(hmm_dir, _hmm_tmpdir)
    logger.info("Searching %d HMM profiles against %d samples", len(profiles), len(samples))

    all_tblout_files = []
    for sample in samples:
        sd = sdirs[sample.sample_id]
        tblouts = run_hmmsearch_all_profiles(
            sample, profiles, sd["hmm_results"], cpus=args.cpus
        )
        for t in tblouts:
            all_tblout_files.append((t, sample.sample_id))

    # Per-sample combined hits
    hmm_hits_dfs = []
    for sample in samples:
        sd = sdirs[sample.sample_id]
        sample_tblouts = [(t, sid) for t, sid in all_tblout_files if sid == sample.sample_id]
        if sample_tblouts:
            hits_path = os.path.join(sd["hmm"], "hmm_hits.tsv")
            df = combine_and_filter_hmmsearch(
                sample_tblouts, hits_path, bitscore_threshold=args.bitscore
            )
            hmm_hits_dfs.append(df)

    hmm_hits_df = pd.concat(hmm_hits_dfs, ignore_index=True) if hmm_hits_dfs else pd.DataFrame()

    # ==================================================================
    # Step 3: Extract flanking DNA around each HMM hit
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 3: Extracting ±%d bp flanking DNA around HMM hits", args.flank_bp)
    logger.info("=" * 60)

    all_flanking_regions = []  # (sample, fasta_path, regions)
    all_flanking_flat = []

    for sample in samples:
        sd = sdirs[sample.sample_id]
        sample_hits = hmm_hits_df[hmm_hits_df["genome"] == sample.sample_id] if not hmm_hits_df.empty else hmm_hits_df
        fasta_path, regions = extract_flanking_regions(
            sample, sample_hits, sd["flanking"], flank_bp=args.flank_bp
        )
        all_flanking_regions.append((sample, fasta_path, regions))
        all_flanking_flat.extend(regions)
        if regions:
            flanking_tsv = os.path.join(sd["flanking"], "flanking_regions.tsv")
            flanking_regions_to_tsv(regions, flanking_tsv)

    # ==================================================================
    # Step 4: Inverted repeat detection (einverted) on flanking regions
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 4: Inverted repeat detection (flanking regions)")
    logger.info("=" * 60)

    ir_results = []  # (ir_dir, sample_id, flanking_regions)
    for sample, fasta_path, regions in all_flanking_regions:
        if not regions:
            continue
        sd = sdirs[sample.sample_id]
        ir_dir = detect_inverted_repeats(
            sample, fasta_path, sd["ir"], cpus=args.cpus,
        )
        ir_results.append((ir_dir, sample.sample_id, regions))

    # Combine IR tables per sample (remapped to genome coords)
    ir_dfs = []
    for ir_dir, sample_id, regions in ir_results:
        combined_path = os.path.join(sdirs[sample_id]["ir"], "IRs_remapped.tsv")
        df = combine_ir_tables([(ir_dir, sample_id, regions)], combined_path)
        ir_dfs.append(df)
    ir_df = pd.concat(ir_dfs, ignore_index=True) if ir_dfs else pd.DataFrame()

    # Apply IR quality filters (arm length / identity)
    if not ir_df.empty and (args.min_ir_arm_length > 0 or args.max_ir_arm_length > 0 or args.min_ir_identity > 0):
        ir_df = filter_ir_table(
            ir_df,
            min_arm_length=args.min_ir_arm_length,
            max_arm_length=args.max_ir_arm_length,
            min_identity=args.min_ir_identity,
        )
        # Write filtered table per sample
        for sample_id, group in ir_df.groupby("sample_id"):
            filtered_path = os.path.join(sdirs[sample_id]["ir"], "IRs_filtered.tsv")
            group.to_csv(filtered_path, sep="\t", index=False)
        logger.info("Filtered IR table: %d records", len(ir_df))

    # Remove IRs whose arms fall inside a recombinase CDS
    if not ir_df.empty:
        ir_df = filter_irs_in_hmm_hits(ir_df, all_flanking_flat)

    # ==================================================================
    # Step 5: Shufflon candidate filtering (density + CDS overlap)
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 5: Shufflon candidate filtering")
    logger.info("=" * 60)

    # Keep a copy of all quality-filtered IRs before the strict shufflon
    # filter.  IRs that don't pass the shufflon criteria are still
    # candidates for inverton-like classification in step 7.
    quality_ir_df = ir_df.copy() if not ir_df.empty else ir_df

    ir_df = filter_shufflon_candidates(
        ir_df,
        samples,
        cluster_distance=args.cluster_distance,
        min_ir_pairs=args.min_ir_pairs,
        min_ir_density=args.min_ir_density,
        window_size=args.window_size,
    )

    # Motif-based refinement: detect sfx sites missed by einverted
    all_sequences = {}
    for sample in samples:
        if sample.gff_path and os.path.isfile(sample.gff_path):
            seqs = parse_fasta_from_gff(sample.gff_path)
            all_sequences.update(seqs)
    if all_sequences:
        ir_df = refine_sfx_sites(ir_df, all_sequences)

    # Write candidate IRs per sample
    for sample_id, group in ir_df.groupby("sample_id"):
        cand_path = os.path.join(sdirs[sample_id]["shufflon"], "IRs_shufflon_candidates.tsv")
        group.to_csv(cand_path, sep="\t", index=False)
    logger.info("Shufflon candidate IRs: %d records", len(ir_df))

    # Inverton candidates: quality-filtered IRs NOT in shufflon clusters.
    # These lack a cluster_id in the GFF, so extract_shufflon_windows
    # will group them by proximity and apply inverton-specific criteria.
    shufflon_ir_df = ir_df
    if not shufflon_ir_df.empty and not quality_ir_df.empty:
        shufflon_keys = set(zip(
            shufflon_ir_df["sample_id"],
            shufflon_ir_df["IR_Chr"],
            shufflon_ir_df["LeftIRStart"],
            shufflon_ir_df["RightIRStart"],
        ))
        quality_keys = list(zip(
            quality_ir_df["sample_id"],
            quality_ir_df["IR_Chr"],
            quality_ir_df["LeftIRStart"],
            quality_ir_df["RightIRStart"],
        ))
        inverton_mask = [k not in shufflon_keys for k in quality_keys]
        inverton_ir_df = quality_ir_df[inverton_mask].copy()
    else:
        inverton_ir_df = quality_ir_df.copy() if not quality_ir_df.empty else pd.DataFrame()

    logger.info("Inverton candidate IRs: %d records", len(inverton_ir_df))

    # Combined IR DataFrame for GFF generation (shufflon IRs carry
    # cluster_id; inverton IRs do not — this distinction drives the
    # classification logic in extract_shufflon_windows).
    combined_ir_df = pd.concat(
        [shufflon_ir_df, inverton_ir_df], ignore_index=True,
    ) if not shufflon_ir_df.empty or not inverton_ir_df.empty else pd.DataFrame()

    # ==================================================================
    # Step 6: Generate and merge GFF files
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 6: GFF generation and merging")
    logger.info("=" * 60)

    for sample in samples:
        sd = sdirs[sample.sample_id]

        # HMM hit GFF
        sample_hits = hmm_hits_df[hmm_hits_df["genome"] == sample.sample_id] if not hmm_hits_df.empty else hmm_hits_df
        hmm_gff_map = hmm_hits_to_gff(sample_hits, [sample], sd["gff_hmm"])

        # IR GFF (both shufflon + inverton candidates)
        sample_irs = combined_ir_df[combined_ir_df["sample_id"] == sample.sample_id] if not combined_ir_df.empty else combined_ir_df
        ir_gff_map = ir_to_gff(sample_irs, sd["gff_ir"])

        # Merge into Prokka GFF
        extra_gffs = []
        if sample.sample_id in hmm_gff_map:
            extra_gffs.append(hmm_gff_map[sample.sample_id])
        if sample.sample_id in ir_gff_map:
            extra_gffs.append(ir_gff_map[sample.sample_id])

        merged_path = os.path.join(sd["gff_merged"], f"{sample.sample_id}.gff")
        merge_gff_into_prokka(sample.gff_path, extra_gffs, merged_path)

    # ==================================================================
    # Step 7: Extract shufflon windows
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 7: Extracting shufflon windows")
    logger.info("=" * 60)

    all_shufflon_windows = []
    all_inverton_windows = []
    for sample in samples:
        sd = sdirs[sample.sample_id]
        merged_gff = os.path.join(sd["gff_merged"], f"{sample.sample_id}.gff")
        if not os.path.isfile(merged_gff):
            logger.warning("No merged GFF for %s, skipping window extraction", sample.sample_id)
            continue

        shuf_wins, inv_wins = extract_shufflon_windows(
            merged_gff,
            sd["shufflon_gffs"],
            sd["inverton_gffs"],
            sample_id=sample.sample_id,
            window_size=args.window_size,
            min_ir_pairs=args.min_ir_pairs,
        )
        all_shufflon_windows.extend(shuf_wins)
        all_inverton_windows.extend(inv_wins)

        # Summary per sample
        if shuf_wins:
            shufflon_windows_to_tsv(
                shuf_wins,
                os.path.join(sd["windows"], "shufflon_like_summary.tsv"),
            )
        if inv_wins:
            shufflon_windows_to_tsv(
                inv_wins,
                os.path.join(sd["windows"], "inverton_like_summary.tsv"),
            )

    # ==================================================================
    # Step 8: Generate shufflon plots
    # ==================================================================
    logger.info("=" * 60)
    logger.info("STEP 8: Generating shufflon plots")
    logger.info("=" * 60)

    total_plots = 0
    for sample in samples:
        sd = sdirs[sample.sample_id]
        plot_files = generate_shufflon_plots(sd["shufflon_gffs"], sd["shufflon_plots"])
        plot_files += generate_shufflon_plots(sd["inverton_gffs"], sd["inverton_plots"])
        total_plots += len(plot_files)
    logger.info("Generated %d plot file(s)", total_plots)

    # ==================================================================
    # Step 9: KOfamscan annotation
    # ==================================================================
    # Determine whether to run KOfamscan based on CLI args and sample sheet.
    # KOfamscan runs when --ko-profiles-dir is given (and --skip-kofamscan
    # is not set), OR is skipped entirely if every sample already has a
    # kofamscan_path from the sample sheet.
    has_ko_profiles_dir = bool(args.ko_profiles_dir)
    any_sample_has_kofamscan = any(s.kofamscan_path for s in samples)
    ko_enabled = has_ko_profiles_dir or any_sample_has_kofamscan

    # Per-sample KO DataFrames
    ko_dfs: dict[str, pd.DataFrame] = {}

    if ko_enabled:
        logger.info("=" * 60)
        logger.info("STEP 9: KOfamscan annotation")
        logger.info("=" * 60)

        for sample in samples:
            sd = sdirs[sample.sample_id]

            # Decide where the KOfamscan output lives
            ko_path = sample.kofamscan_path
            if ko_path and os.path.isfile(ko_path):
                logger.info(
                    "Using pre-computed KOfamscan output for %s: %s",
                    sample.sample_id, ko_path,
                )
            elif run_kofamscan_step:
                ko_path = run_kofamscan(
                    sample, sd["kofamscan"], args.ko_profiles_dir,
                    cpus=args.cpus,
                )
            else:
                logger.info(
                    "No KOfamscan output for %s and --ko-profiles-dir not set, skipping",
                    sample.sample_id,
                )
                ko_dfs[sample.sample_id] = pd.DataFrame(
                    columns=["sample_id", "gene_id", "ko_accession", "threshold", "score", "e_value"]
                )
                continue

            ko_dfs[sample.sample_id] = load_ko_table(sample, ko_path)

        # ==================================================================
        # Step 10: KEGG enrichment analysis
        # ==================================================================
        logger.info("=" * 60)
        logger.info("STEP 10: KEGG enrichment analysis (inverton-like windows)")
        logger.info("=" * 60)

        # Build per-sample window lists
        sample_shufflon_wins: dict[str, list] = defaultdict(list)
        sample_inverton_wins: dict[str, list] = defaultdict(list)
        for win in all_shufflon_windows:
            sample_shufflon_wins[win.sample_id].append(win)
        for win in all_inverton_windows:
            sample_inverton_wins[win.sample_id].append(win)

        for sample in samples:
            sid = sample.sample_id
            sd = sdirs[sid]
            ko_df = ko_dfs.get(sid, pd.DataFrame())

            if ko_df.empty:
                logger.info("No KO data for %s, skipping enrichment", sid)
                continue

            # Background: all KO accessions in the genome
            background_kos = set(ko_df["ko_accession"].dropna()) - {""}

            s_wins = sample_shufflon_wins.get(sid, [])
            i_wins = sample_inverton_wins.get(sid, [])

            if not s_wins and not i_wins:
                logger.info("No windows for %s, skipping enrichment", sid)
                continue

            run_kegg_enrichment_for_sample(
                sid, s_wins, i_wins, ko_df, background_kos, sd["kegg_enrich"],
            )
    else:
        logger.info("STEPS 9-10: KOfamscan / KEGG enrichment — skipped (no --ko-profiles-dir)")

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
    logger.info("Shufflon candidate IRs:%d", len(shufflon_ir_df))
    logger.info("Inverton candidate IRs:%d", len(inverton_ir_df))
    logger.info("Shufflon-like windows: %d", len(all_shufflon_windows))
    logger.info("Inverton-like windows: %d", len(all_inverton_windows))
    if ko_enabled:
        total_kos = sum(len(df) for df in ko_dfs.values())
        logger.info("KO assignments:        %d", total_kos)
    logger.info("Results in:            %s", outdir)

    # Clean up temporary HMM profile directory
    shutil.rmtree(_hmm_tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
