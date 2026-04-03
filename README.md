# shufflonfinder

Annotate shufflon structures in bacterial genomes. Searches predicted proteins against a library of HMM profiles targeting shufflon-associated recombinases (e.g. Rci), extracts flanking DNA around each hit, detects inverted repeats in those flanking regions using EMBOSS einverted, filters for dense IR clusters where at least one repeat arm overlaps a coding sequence, and produces merged GFF annotations with windowed output for each candidate shufflon locus.


## How it works

The pipeline runs seven steps in sequence:

1. **Prokka** annotates raw genome FASTAs to produce protein sequences (.faa), gene coordinates (.gff), and contigs (.fna). Skipped when you supply pre-annotated inputs via a sample sheet.

2. **HMM search** decompresses all `.hmm.gz` profiles from the bundled `hmms/` directory, prepares them with `hmmpress`, and runs `hmmsearch` for each profile against each sample's proteins. The profiles target shufflon-associated recombinases (Rci family) and related mobile element components. Results are pooled and filtered by a bitscore threshold. Proteins matching multiple profiles are deduplicated so each protein produces one entry (with all matching profiles recorded).

3. **Flanking extraction** maps each hit protein back to its genomic coordinates through the Prokka GFF, then pulls ±5 kb of DNA (configurable) from the genome FASTA on each side of the CDS. The output is a multi-record FASTA plus a metadata TSV tracking the coordinate mappings.

4. **Inverted repeat detection** runs EMBOSS `einverted` at two sensitivity thresholds (51 and 75) on the flanking FASTAs, merges the results with coordinate-based deduplication (preserving distinct pairs in dense clusters), extracts arm sequences, and computes percent identity between arms. IR coordinates are then translated from flanking-region-local back to genome-absolute positions.

5. **Shufflon candidate filtering** clusters nearby IRs on each contig within a configurable window (default 3 kb). Clusters with fewer than `--min-ir-pairs` (default 2) IR pairs are discarded. Remaining clusters must have at least one IR arm overlapping a coding sequence. This selects for the dense repeat architecture characteristic of shufflons while filtering out isolated or intergenic repeats.

6. **GFF generation and merging** converts both HMM hits and the surviving IRs into GFF3 features, then inserts them into the Prokka GFF between the annotation lines and the embedded FASTA section.

7. **Window extraction** identifies CDS features overlapping each shufflon candidate cluster, finds invertible DNA segments between IR arms, predicts ORFs within those segments (to recover small cassette genes that Prokka may have missed), and writes self-contained GFF+FASTA files for each candidate shufflon region. Each window GFF has five annotation tracks: inverted repeats, HMM hits, Prokka CDS, invertible segments, and predicted shufflon cassette ORFs.


## Prerequisites

[Conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/). Everything else is installed automatically.

External tools (installed via conda):

- [Prokka](https://github.com/tseemann/prokka) >= 1.14
- [HMMER](http://hmmer.org/) >= 3.3
- [EMBOSS](http://emboss.sourceforge.net/) >= 6.6 (provides `einverted`)

Python libraries (installed via conda):

- Python >= 3.9
- Biopython >= 1.80
- pandas >= 1.5


## Installation

```bash
git clone https://github.com/your-org/shufflonfinder.git
cd shufflonfinder

conda env create -f environment.yml
conda activate shufflonfinder
```

The `environment.yml` installs all conda packages, pip packages, and shufflonfinder itself (via `pip install -e .`) in one step. After activation, the `shufflonfinder` command is available on your PATH.

Verify everything installed correctly:

```bash
shufflonfinder --help
prokka --version
hmmsearch -h | head -1
einverted --help
```

The 41 HMM profiles ship with the package in `shufflonfinder/hmms/`. No additional downloads are needed.


## Quick start

```bash
# Single genome
shufflonfinder --input-fasta my_genome.fna --outdir results/

# Directory of genomes, 8 threads per tool
shufflonfinder --input-fasta genomes/ --outdir results/ --cpus 8

# Pre-annotated samples (skip Prokka)
shufflonfinder --sample-sheet samples.tsv --outdir results/
```


## Usage

### From raw genome FASTAs

Pass a single FASTA file or a directory of FASTAs. Prokka runs automatically.

```bash
shufflonfinder \
    --input-fasta genomes/ \
    --outdir results/ \
    --cpus 8
```

Recognized FASTA extensions: `.fasta`, `.fa`, `.fna` (and `.gz` versions of each).

### From pre-annotated samples

If Prokka (or a compatible annotator) has already been run, provide a tab-separated sample sheet with these columns:

```
sample_id	fna_path	faa_path	gff_path
genome_001	/data/genome_001.fna	/data/genome_001.faa	/data/genome_001.gff
genome_002	/data/genome_002.fna	/data/genome_002.faa	/data/genome_002.gff
```

```bash
shufflonfinder \
    --sample-sheet samples.tsv \
    --outdir results/
```

All paths in the sample sheet must be absolute or resolvable from the working directory.

### Custom HMM profiles

To use your own profiles instead of the bundled set:

```bash
shufflonfinder \
    --input-fasta genomes/ \
    --hmm-dir /path/to/my_hmms/ \
    --outdir results/
```

The directory can contain `.hmm` or `.hmm.gz` files. Each file should hold one HMM profile.


## Options

```
--input-fasta PATH    Genome FASTA file or directory (mutually exclusive with --sample-sheet)
--sample-sheet PATH   TSV with columns: sample_id, fna_path, faa_path, gff_path
--outdir PATH         Output directory (required)
--hmm-dir PATH        Directory of .hmm/.hmm.gz profiles (default: bundled profiles)
--cpus INT            Threads per tool invocation (default: 4)
--bitscore FLOAT      Minimum HMM bitscore to keep a hit (default: 25.0)
--flank-bp INT        DNA to extract on each side of a hit protein, in bp (default: 5000)
--window-size INT     Max distance in bp for window extraction context (default: 3000)
--min-ir-arm-length INT  Minimum IR arm length in bp to keep (default: 10)
--min-ir-identity FLOAT  Minimum percent identity between IR arms (default: 70.0)
--min-ir-pairs INT    Minimum IR pairs per cluster to qualify as shufflon candidate (default: 3)
--cluster-distance INT  Max gap between IR pairs for chaining into one cluster (default: 1000)
--min-ir-density FLOAT  Minimum IR pairs per kilobase in a cluster (default: 2.0)
--skip-prokka         Skip Prokka even for FASTA inputs
-v                    Verbose output (use -vv for debug)
-q                    Quiet mode (errors only)
```


## Output structure

```
results/
├── 01_prokka/                      # Prokka outputs per sample (when run)
├── 02_hmmsearch/
│   ├── profiles/                   # Decompressed, hmmpress'd profiles
│   ├── results/<sample_id>/        # Per-sample, per-profile tblout files
│   └── hmm_hits_combined.tsv       # All hits above the bitscore threshold
├── 03_flanking/
│   ├── <sample_id>_flanking.fasta  # Flanking DNA around each hit protein
│   └── flanking_regions_combined.tsv
├── 04_inverted_repeats/
│   ├── <sample_id>/                # einverted output + IRs.tsv per sample
│   ├── IRs_combined_remapped.tsv   # All IRs in genome-absolute coordinates
│   └── IRs_combined_filtered.tsv   # After arm-length and identity filters
├── 05_shufflon_filter/
│   └── IRs_shufflon_candidates.tsv # After density + CDS overlap filtering
├── 06_gff/
│   ├── hmm_hits/                   # HMM hit CDS features as GFF
│   ├── ir/                         # IR features as GFF
│   └── merged/                     # Prokka + HMM + IR merged GFFs
└── 07_shufflon_windows/
    ├── shufflon_windows_summary.tsv  # Combined summary table (all samples)
    └── <sample_id>/                  # Per-window GFF+FASTA files
```

### Key output files

`02_hmmsearch/hmm_hits_combined.tsv` lists every protein hit across all samples and profiles. Columns include `target_name` (protein ID), `hmm_profile`, `full_sequence_bitscore`, and `genome` (sample ID).

`03_flanking/flanking_regions_combined.tsv` records the genomic coordinates of each flanking region, the CDS it surrounds, which HMM profiles matched, and the extracted sequence length.

`04_inverted_repeats/IRs_combined_remapped.tsv` contains all detected inverted repeats with coordinates translated back to the original contigs, plus arm sequences and percent identity.

`06_shufflon_windows/shufflon_windows_summary.tsv` is the main results table. Each row represents one CDS feature within a candidate shufflon window. Columns include `sample_id`, `window_id`, `contig`, `window_start`, `window_end`, `window_length_bp`, `n_ir_pairs`, `ir_coords` (compact coordinate string for all IR pairs in the window), `locus_tag` (Prokka ID), `cds_start`, `cds_end`, `strand`, `product` (Prokka annotation), `cds_source`, `is_hmm_hit` (True if this CDS is the gene that triggered the HMM search), `hmm_profiles` (semicolon-separated profiles that matched, empty for non-hit genes), and `gff_path` (path to the per-window GFF+FASTA file).

`07_shufflon_windows/<sample_id>/` contains one GFF+FASTA file per candidate shufflon region. Each GFF has five annotation tracks: `inverted_repeat` features (from einverted), `hmm_hit` features (the recombinase gene with HMM profile info), `CDS` features (overlapping Prokka genes), `invertible_segment` features (the DNA between consecutive IR arms, corresponding to shufflon cassettes), and predicted `CDS` features from shufflonfinder's own ORF finder for small cassette genes that Prokka may have missed (e.g. the PilV C-terminal segments in R64-type shufflons).


## HMM profiles

The bundled profiles (41 total) come from Pfam, PANTHER, TIGRFAM, Gene3D, PIRSF, and other databases, selected for their association with shufflon-specific recombinases (Rci family), site-specific DNA invertases, and related mobile genetic element components. These recombinases catalyze the inversions at sfx recognition sites that define shufflon activity.

A protein counts as a hit if it scores at or above `--bitscore` against any profile. Proteins matching multiple profiles are deduplicated at the flanking extraction step so each genomic locus is scanned for IRs exactly once.


## Inverted repeat filtering

einverted is run at two sensitivity thresholds (51 and 75) with different scoring matrices. Results from both runs are merged using coordinate-based deduplication: a threshold-51 pair is discarded only if a threshold-75 pair exists with all four arm coordinates within 3 bp (i.e. the same detection at different sensitivity). Distinct pairs that happen to share a dense genomic region are preserved. This matters for shufflons, which can pack multiple sfx recognition sites into a few kilobases. A hardcoded minimum of 30 bp between the two arms filters out trivially close pairs.

After detection, two configurable filters are applied:

`--min-ir-arm-length` sets the minimum length (in bp) for both the left and right arms of an IR pair. Short arms are more likely to be noise. The default is 10 bp.

`--min-ir-identity` sets the minimum percent identity between the two arms. Identity is computed from the extracted arm sequences (left arm vs. reverse complement of right arm). The default is 70%.

Both the unfiltered (`IRs_combined_remapped.tsv`) and filtered (`IRs_combined_filtered.tsv`) tables are written to the `04_inverted_repeats/` output directory, so you can always inspect what was removed.

To disable filtering entirely, set both thresholds to zero:

```bash
shufflonfinder --input-fasta genomes/ --outdir results/ --min-ir-arm-length 0 --min-ir-identity 0
```


## Shufflon candidate filtering

After the basic arm-length and identity filters, a second stage selects for the dense inverted repeat clusters characteristic of shufflons. IRs on the same contig are clustered by proximity using `--cluster-distance` (default 1000bp): two IR pairs belong to the same cluster if any of their arms is within this distance. This is deliberately tighter than `--window-size` (used later for window extraction context) to isolate the dense shufflon core from neighboring transposon or IS-element repeats that happen to share the same genomic region.

Clusters must pass three tests:

`--min-ir-pairs` (default 3) sets the minimum number of IR pairs. A shufflon with two invertible segments has at least three recognition sites, which produce at least three detected IR pairs.

`--min-ir-density` (default 2.0 pairs/kb) filters out clusters where the IRs are spread too thinly. Characterized shufflons like R64's pack 7 sfx sites into ~2kb (3.5 pairs/kb). Regions where transposon IRs happen to accumulate near a recombinase gene rarely exceed 1.5 pairs/kb because the individual transposon IRs are spaced further apart.

At least one IR arm in the cluster must overlap a coding sequence. In characterized shufflons, inverted repeats are found within and flanking the invertible gene segments. Clusters where all repeats fall in intergenic regions are discarded.

The unfiltered and basic-filtered tables are in `04_inverted_repeats/`, while the shufflon-filtered table is in `05_shufflon_filter/`, so you can inspect what was removed at each stage.


## Project layout

```
shufflonfinder/
├── shufflonfinder/              # Python package
│   ├── __init__.py
│   ├── cli.py                   # CLI entry point (console_scripts)
│   ├── utils.py                 # Logging, shell commands, path helpers
│   ├── sample_sheet.py          # Input parsing (FASTA, sample sheet)
│   ├── step_prokka.py           # Prokka wrapper
│   ├── step_hmmsearch.py        # Multi-profile HMM search
│   ├── step_flanking.py         # Flanking DNA extraction + deduplication
│   ├── step_phava.py            # einverted IR detection + coordinate remapping
│   ├── step_ir_cds.py           # IR-CDS containment annotation
│   ├── step_gff.py              # GFF generation, merging, window extraction
│   └── hmms/                    # Bundled HMM profiles (.hmm.gz)
├── pyproject.toml               # Package metadata + console_scripts
├── environment.yml              # Conda + pip environment
├── README.md
└── .gitignore
```


## License

MIT
