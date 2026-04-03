# shufflonfinder

Annotate shufflon structures in bacterial genomes. Searches predicted proteins against a library of shufflon-associated HMM profiles, extracts flanking DNA around each hit, detects inverted repeats in those flanking regions using EMBOSS einverted, and produces merged GFF annotations with windowed output for each candidate shufflon locus.


## How it works

The pipeline runs six steps in sequence:

1. **Prokka** annotates raw genome FASTAs to produce protein sequences (.faa), gene coordinates (.gff), and contigs (.fna). Skipped when you supply pre-annotated inputs via a sample sheet.

2. **HMM search** decompresses all `.hmm.gz` profiles from the bundled `hmms/` directory, prepares them with `hmmpress`, and runs `hmmsearch` for each profile against each sample's proteins. Results are pooled and filtered by a bitscore threshold. Proteins matching multiple profiles are deduplicated so each protein produces one entry (with all matching profiles recorded).

3. **Flanking extraction** maps each hit protein back to its genomic coordinates through the Prokka GFF, then pulls ±5 kb of DNA (configurable) from the genome FASTA on each side of the CDS. The output is a multi-record FASTA plus a metadata TSV tracking the coordinate mappings.

4. **Inverted repeat detection** runs EMBOSS `einverted` at two sensitivity thresholds (51 and 75) on the flanking FASTAs, merges the results with overlap deduplication, extracts arm sequences, and computes percent identity between arms. IR coordinates are then translated from flanking-region-local back to genome-absolute positions.

5. **GFF generation and merging** converts both HMM hits and IRs into GFF3 features, then inserts them into the Prokka GFF between the annotation lines and the embedded FASTA section.

6. **Window extraction** clusters nearby IRs within a configurable window (default 3 kb), identifies CDS features overlapping each cluster, and writes self-contained GFF+FASTA files for each candidate shufflon region.


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
--window-size INT     Max distance in bp for clustering nearby IRs (default: 3000)
--min-ir-arm-length INT  Minimum IR arm length in bp to keep (default: 10)
--min-ir-identity FLOAT  Minimum percent identity between IR arms (default: 70.0)
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
│   └── IRs_combined_remapped.tsv   # IRs in genome-absolute coordinates
├── 05_gff/
│   ├── hmm_hits/                   # HMM hit CDS features as GFF
│   ├── ir/                         # IR features as GFF
│   └── merged/                     # Prokka + HMM + IR merged GFFs
└── 06_shufflon_windows/
    ├── shufflon_windows_summary.tsv  # Combined summary table (all samples)
    └── <sample_id>/                  # Per-window GFF+FASTA files
```

### Key output files

`02_hmmsearch/hmm_hits_combined.tsv` lists every protein hit across all samples and profiles. Columns include `target_name` (protein ID), `hmm_profile`, `full_sequence_bitscore`, and `genome` (sample ID).

`03_flanking/flanking_regions_combined.tsv` records the genomic coordinates of each flanking region, the CDS it surrounds, which HMM profiles matched, and the extracted sequence length.

`04_inverted_repeats/IRs_combined_remapped.tsv` contains all detected inverted repeats with coordinates translated back to the original contigs, plus arm sequences and percent identity.

`06_shufflon_windows/shufflon_windows_summary.tsv` is the main results table. Each row represents one CDS feature within a candidate shufflon window. Columns include `sample_id`, `window_id`, `contig`, `window_start`, `window_end`, `window_length_bp`, `n_ir_pairs`, `ir_coords` (compact coordinate string for all IR pairs in the window), `locus_tag` (Prokka ID), `cds_start`, `cds_end`, `strand`, `product` (Prokka annotation), `cds_source`, and `gff_path` (path to the per-window GFF+FASTA file).

`06_shufflon_windows/<sample_id>/` contains one GFF+FASTA file per candidate shufflon region, with IR features, overlapping CDS features, and the corresponding DNA sequence.


## HMM profiles

The bundled profiles (41 total) come from Pfam, PANTHER, TIGRFAM, Gene3D, PIRSF, and other databases, selected for their association with shufflon recombinases, pilus tip adhesins, and related mobile genetic element components.

A protein counts as a hit if it scores at or above `--bitscore` against any profile. Proteins matching multiple profiles are deduplicated at the flanking extraction step so each genomic locus is scanned for IRs exactly once.


## Inverted repeat filtering

einverted is run at two sensitivity thresholds (51 and 75) with different scoring matrices. Results from both runs are merged, with threshold-75 hits taking priority: threshold-51 hits are only kept if they don't overlap any threshold-75 hit on the same contig. A hardcoded minimum of 30 bp between the two arms filters out trivially close pairs.

After detection, two configurable filters are applied:

`--min-ir-arm-length` sets the minimum length (in bp) for both the left and right arms of an IR pair. Short arms are more likely to be noise. The default is 10 bp.

`--min-ir-identity` sets the minimum percent identity between the two arms. Identity is computed from the extracted arm sequences (left arm vs. reverse complement of right arm). The default is 70%.

Both the unfiltered (`IRs_combined_remapped.tsv`) and filtered (`IRs_combined_filtered.tsv`) tables are written to the `04_inverted_repeats/` output directory, so you can always inspect what was removed.

To disable filtering entirely, set both thresholds to zero:

```bash
shufflonfinder --input-fasta genomes/ --outdir results/ --min-ir-arm-length 0 --min-ir-identity 0
```


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
