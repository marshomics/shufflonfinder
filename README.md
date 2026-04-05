# ShufflonFinder

<p align="center">
  <img src="shufflon_pipeline.png" alt="shufflonfinder pipeline overview" width="100%"/>
</p>

Annotate shufflon structures in bacterial genomes. Searches predicted proteins against a library of HMM profiles targeting shufflon-associated recombinases (e.g. Rci), extracts flanking DNA around each hit, detects inverted repeats in those flanking regions using EMBOSS einverted, applies motif-based refinement to recover sfx recognition sites missed by pairwise alignment, filters for dense IR clusters, and produces merged GFF annotations with windowed output and gene-organisation plots for each candidate shufflon locus.


## How it works

The pipeline runs up to ten steps in sequence:

1. **Prokka** annotates raw genome FASTAs to produce protein sequences (.faa), gene coordinates (.gff), and contigs (.fna). Skipped when you supply pre-annotated inputs via a sample sheet.

2. **HMM search** decompresses all `.hmm.gz` profiles from the bundled `hmms/` directory, prepares them with `hmmpress`, and runs `hmmsearch` for each profile against each sample's proteins. The profiles target shufflon-associated recombinases (Rci family) and related mobile element components. Results are pooled and filtered by a bitscore threshold. Proteins matching multiple profiles are deduplicated so each protein produces one entry (with all matching profiles recorded).

3. **Flanking extraction** maps each hit protein back to its genomic coordinates through the Prokka GFF, then pulls ±5 kb of DNA (configurable) from the genome FASTA on each side of the CDS. The output is a multi-record FASTA plus a metadata TSV tracking the coordinate mappings.

4. **Inverted repeat detection** runs EMBOSS `einverted` at two sensitivity thresholds (51 and 75) on the flanking FASTAs, merges the results with coordinate-based deduplication (preserving distinct pairs in dense clusters), extracts arm sequences, and computes percent identity between arms. IR coordinates are then translated from flanking-region-local back to genome-absolute positions.

5. **Shufflon candidate filtering** clusters nearby IRs on each contig using density-based gap chaining with `--cluster-distance` (default 1000 bp). Clusters with fewer than `--min-ir-pairs` (default 3) IR pairs or density below `--min-ir-density` (default 1.0 pairs/kb) are discarded. Remaining clusters must have at least one IR arm overlapping a coding sequence.

6. **Motif-based sfx refinement** derives a core recognition motif from the arms detected in each cluster, then scans the surrounding genomic sequence for additional sfx sites that einverted missed. einverted finds IR pairs by pairwise alignment, so it only detects sfx combinations with sufficient identity to each other; cross-type sfx pairings (e.g. sfxA vs. sfxD') that share the conserved core but diverge in flanking bases are often missed. Extended pairs — einverted artifacts where an arm bridges across two adjacent sfx sites — are identified using a median-based arm-length ratio test (>1.4× the cluster median) and replaced with the individual sites they span. Newly discovered sites that lack a detected partner are stored as unpaired recognition sites.

7. **GFF generation and merging** converts both HMM hits and the surviving IRs into GFF3 features, then inserts them into the Prokka GFF between the annotation lines and the embedded FASTA section.

8. **Window extraction and classification** groups IR features by cluster (shufflon candidates carry a `cluster_id`; remaining quality-filtered IRs are grouped by proximity) and classifies each group into one of two categories:

   * **Shufflon-like** — at least `--min-ir-pairs` (default 3) complete pairs *and* every consecutive pair (sorted by position) is adjacent within 10 bp. This captures the R64/TP114-style organisation where sfx recognition sites sit at cassette boundaries in quick succession.
   * **Inverton-like** — at least one complete IR pair within 2000 bp of an HMM-hit CDS. This captures simpler invertible elements near a recombinase that don't have the tightly packed multi-cassette layout of a shufflon.

   Windows matching neither category are skipped. Each window GFF carries four annotation tracks: `inverted_repeat`, `hmm_hit`, `CDS`, and `invertible_segment`. Shufflon-like and inverton-like windows are written to separate subdirectories under `07_shufflon_windows/`.

9. **Plot generation** produces PNG and SVG gene-organisation figures for each window GFF (both shufflon-like and inverton-like) using [dna_features_viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer). Each plot has three visual tracks: inverted repeats as directional arrows above the annotation line, CDS and recombinase as directional arrows on the main line, and invertible segments as undirectional boxes below. Each IR pair gets a distinct colour (parsed from the feature name, e.g. `inverted_repeat_01_FOR`/`inverted_repeat_01_REV`), so overlapping pairs are visually distinguishable even when drawn on the same line. Other features are colour-coded by category: recombinases (red), invertible segments (orange), CDS containing inverted repeats (teal), and other CDS (grey). Identification is via legend only.

10. **KOfamscan on IR-overlapping CDS** (optional, requires `--ko-profiles-dir`) identifies CDS that overlap at least one inverted repeat arm within each shufflon-like or inverton-like window, extracts those protein sequences from the Prokka `.faa`, and runs [KOfamscan](https://github.com/takaram/kofam_scan) (`exec_annotation`) on the extracted subset. The results are combined into a per-sample table listing each IR-overlapping CDS with its Prokka annotation, KO accession (where assigned), and window category (`shufflon_like` or `inverton_like`). This table is written to `07_shufflon_windows/`.


## Prerequisites

[Conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/). Everything else is installed automatically.

External tools (installed via conda):

- [Prokka](https://github.com/tseemann/prokka) >= 1.14
- [HMMER](http://hmmer.org/) >= 3.3
- [EMBOSS](http://emboss.sourceforge.net/) >= 6.6 (provides `einverted`)

Optional external tools (for KO annotation of IR-containing CDS):

- [KOfamscan](https://github.com/takaram/kofam_scan) — provides `exec_annotation`. Required only when using `--ko-profiles-dir`. Follow the [KOfamscan README](https://github.com/takaram/kofam_scan#usage) to download the KO profiles and ko_list file.

Python libraries (installed via pip):

- Python >= 3.9
- Biopython >= 1.80
- pandas >= 1.5
- [dna_features_viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) >= 3.1
- matplotlib >= 3.5


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
python -c "import dna_features_viewer; print(dna_features_viewer.__version__)"
```

The 41 HMM profiles ship with the package in `shufflonfinder/hmms/`. No additional downloads are needed.

### KOfamscan setup (optional)

KOfamscan is only needed if you want KO annotation of IR-overlapping CDS. Install it and download the profile database:

```bash
# Install KOfamscan (see https://github.com/takaram/kofam_scan)
conda install -c bioconda kofam_scan

# Download KO profiles (this takes a while — ~4 GB)
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
tar xzf profiles.tar.gz
gunzip ko_list.gz

# The resulting directory should contain profiles/ and ko_list:
ls ko_profiles/
# -> ko_list  profiles/
```

Then pass `--ko-profiles-dir ko_profiles/` when running shufflonfinder.


## Quick start

```bash
# Single genome
shufflonfinder --input-fasta my_genome.fna --outdir results/

# Directory of genomes, 8 threads per tool
shufflonfinder --input-fasta genomes/ --outdir results/ --cpus 8

# Pre-annotated samples (skip Prokka)
shufflonfinder --sample-sheet samples.tsv --outdir results/

# With KO annotation of IR-overlapping CDS
shufflonfinder --input-fasta genomes/ --outdir results/ --ko-profiles-dir ko_profiles/
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
--min-ir-arm-length INT   Minimum IR arm length in bp to keep (default: 13)
--max-ir-arm-length INT   Maximum IR arm length in bp to keep (default: 35)
--min-ir-identity FLOAT   Minimum percent identity between IR arms (default: 85.0)
--min-ir-pairs INT        Minimum IR pairs per cluster to qualify as shufflon candidate (default: 3)
--cluster-distance INT    Max gap between IR pairs for chaining into one cluster (default: 1000)
--min-ir-density FLOAT    Minimum IR pairs per kilobase in a cluster (default: 1.0)
--skip-prokka         Skip Prokka even for FASTA inputs
--ko-profiles-dir PATH  KOfamscan profile directory (must contain profiles/ and ko_list)
-v                    Verbose output (use -vv for debug)
-q                    Quiet mode (errors only)
```


## Output structure

Each sample gets its own top-level directory.

```
results/
└── <sample_id>/
    ├── 01_prokka/                             # Prokka outputs (when run)
    │   ├── <sample_id>.faa
    │   ├── <sample_id>.gff
    │   └── ...
    ├── 02_hmmsearch/
    │   ├── results/                           # Per-profile .tblout files
    │   └── hmm_hits.tsv                       # Hits above bitscore threshold
    ├── 03_flanking/
    │   ├── <sample_id>_flanking.fasta         # Flanking DNA FASTA
    │   └── flanking_regions.tsv               # Coordinate mappings
    ├── 04_inverted_repeats/
    │   ├── IRs_remapped.tsv                   # All IRs in genome-absolute coords
    │   └── IRs_filtered.tsv                   # After arm-length + identity filters
    ├── 05_shufflon_filter/
    │   └── IRs_shufflon_candidates.tsv        # After density + CDS + motif refinement
    ├── 06_gff/
    │   ├── hmm_hits/                          # HMM hit features as GFF
    │   ├── ir/                                # IR features as GFF
    │   └── merged/                            # Prokka + HMM + IR merged GFF
    └── 07_shufflon_windows/
        ├── shufflon_like_summary.tsv          # Summary table (shufflon-like windows)
        ├── inverton_like_summary.tsv          # Summary table (inverton-like windows)
        ├── shufflon_like/
        │   ├── gffs/                          # Per-window GFF+FASTA files
        │   │   └── <sample_id>_contig_*_window_*.gff
        │   └── plots/                         # Per-window .png and .svg plots
        │       └── <sample_id>_contig_*_window_*.png
        ├── inverton_like/
        │   ├── gffs/                          # Per-window GFF+FASTA files
        │   │   └── <sample_id>_contig_*_window_*.gff
        │   └── plots/                         # Per-window .png and .svg plots
        │       └── <sample_id>_contig_*_window_*.png
        ├── <sample_id>_ir_cds_ko.tsv          # IR-CDS with KO accessions (when --ko-profiles-dir set)
        ├── <sample_id>_ir_cds.faa             # Extracted IR-CDS proteins
        └── <sample_id>_ir_cds_kofamscan.txt   # Raw KOfamscan output
```

### Key output files

`02_hmmsearch/hmm_hits.tsv` lists every protein hit for the sample across all profiles. Columns include `target_name` (protein ID), `hmm_profile`, `full_sequence_bitscore`, and `genome` (sample ID).

`03_flanking/flanking_regions.tsv` records the genomic coordinates of each flanking region, the CDS it surrounds, which HMM profiles matched, and the extracted sequence length.

`04_inverted_repeats/IRs_remapped.tsv` contains all detected inverted repeats with coordinates translated back to the original contigs, plus arm sequences and percent identity.

`05_shufflon_filter/IRs_shufflon_candidates.tsv` contains the final set of IR pairs after density-based clustering, CDS overlap filtering, and motif-based sfx refinement. Includes `cluster_id` assignments and `unpaired_site` flags for recognition sites discovered by motif search that lack a detected partner arm.

`07_shufflon_windows/shufflon_like_summary.tsv` and `inverton_like_summary.tsv` are the main results tables (one per category). Each row represents one CDS feature within a window. Columns include `sample_id`, `window_id`, `contig`, `window_start`, `window_end`, `window_length_bp`, `n_ir_pairs`, `ir_coords` (compact coordinate string for all IR pairs in the window), `locus_tag` (Prokka ID), `cds_start`, `cds_end`, `strand`, `product` (Prokka annotation), `cds_source`, `is_hmm_hit` (True if this CDS is the gene that triggered the HMM search), `hmm_profiles` (semicolon-separated profiles that matched, empty for non-hit genes), and `gff_path` (path to the per-window GFF+FASTA file).

`07_shufflon_windows/shufflon_like/gffs/` and `inverton_like/gffs/` each contain one GFF+FASTA file per detected window. Each GFF has four annotation tracks: `inverted_repeat` features (from einverted or motif_search, named `inverted_repeat_NN_FOR`/`inverted_repeat_NN_REV`), `hmm_hit` features (the recombinase gene), `CDS` features (overlapping Prokka genes), and `invertible_segment` features (the DNA between consecutive IR arms).

`07_shufflon_windows/shufflon_like/plots/` and `inverton_like/plots/` contain PNG and SVG gene-organisation figures for each window. Each figure has three visual tracks: inverted repeats (each pair in a distinct colour) as directional arrows above the annotation line, CDS and recombinase as directional arrows on the main line, and invertible segments as undirectional boxes below. Features are identified by colour via a shared legend rather than inline labels.

`07_shufflon_windows/<sample_id>_ir_cds_ko.tsv` (present when `--ko-profiles-dir` is set) lists every CDS that overlaps at least one inverted repeat arm within a shufflon-like or inverton-like window. Columns: `sample_id`, `window_id`, `contig`, `locus_tag`, `cds_start`, `cds_end`, `strand`, `product` (Prokka annotation), `ko_accession` (KO ID from KOfamscan, empty when no assignment passes threshold), and `category` (`shufflon_like` or `inverton_like`). KOfamscan is run only on the extracted IR-overlapping proteins, not the full proteome.


## HMM profiles

The bundled profiles (41 total) come from Pfam, PANTHER, TIGRFAM, Gene3D, PIRSF, and other databases, selected for their association with shufflon-specific recombinases (Rci family), site-specific DNA invertases, and related mobile genetic element components. These recombinases catalyse the inversions at sfx recognition sites that define shufflon activity.

A protein counts as a hit if it scores at or above `--bitscore` against any profile. Proteins matching multiple profiles are deduplicated at the flanking extraction step so each genomic locus is scanned for IRs exactly once.


## Inverted repeat filtering

einverted is run at two sensitivity thresholds (51 and 75) with different scoring matrices. Results from both runs are merged using coordinate-based deduplication: a threshold-51 pair is discarded only if a threshold-75 pair exists with all four arm coordinates within 3 bp (i.e. the same detection at different sensitivity). Distinct pairs that happen to share a dense genomic region are preserved. This matters for shufflons, which can pack multiple sfx recognition sites into a few kilobases. A hardcoded minimum of 30 bp between the two arms filters out trivially close pairs.

After detection, three configurable filters are applied:

`--min-ir-arm-length` (default 13 bp) sets the minimum length for both the left and right arms of an IR pair. Shufflon sfx recognition sites are typically 19 bp, but einverted may report shorter partial alignments for genuine sites. The default retains these while filtering noise.

`--max-ir-arm-length` (default 35 bp) sets the maximum arm length. Arms longer than this typically come from transposon or IS-element IRs, or from einverted extending across two abutting sfx sites. The median-based extended pair detection in the motif refinement step handles the latter case more precisely, but this hard ceiling provides an additional filter.

`--min-ir-identity` (default 85.0%) sets the minimum percent identity between the two arms. Within-type sfx pairs are near-perfect reverse complements (≥95%), but cross-type pairs and divergent recognition sites may be lower. The default retains genuine cross-type sfx pairings while filtering most noise.

Both the unfiltered (`IRs_combined_remapped.tsv`) and filtered (`IRs_combined_filtered.tsv`) tables are written to the `04_inverted_repeats/` output directory, so you can always inspect what was removed.


## Shufflon candidate filtering

After the basic arm-length and identity filters, a second stage selects for the dense inverted repeat clusters characteristic of shufflons. IRs on the same contig are clustered by proximity using `--cluster-distance` (default 1000 bp): two IR pairs belong to the same cluster if any of their arms is within this distance. This is deliberately tighter than `--window-size` (used later for window extraction context) to isolate the dense shufflon core from neighbouring transposon or IS-element repeats that happen to share the same genomic region.

Clusters must pass three tests:

`--min-ir-pairs` (default 3) sets the minimum number of IR pairs. A shufflon with two invertible segments has at least three recognition sites, which produce at least three detected IR pairs.

`--min-ir-density` (default 1.0 pairs/kb) filters out clusters where the IRs are spread too thinly. Characterised shufflons like R64's pack 7 sfx sites into ~2 kb (3.5 pairs/kb). The default is set conservatively to accommodate shufflons where einverted detects fewer of the pairwise sfx combinations.

At least one IR arm in the cluster must overlap a coding sequence. In characterised shufflons, inverted repeats are found within and flanking the invertible gene segments. Clusters where all repeats fall in intergenic regions are discarded.


## Motif-based sfx refinement

einverted detects IR pairs by pairwise alignment, which means it only finds sfx site combinations with sufficient mutual identity. In a shufflon with multiple sfx types (e.g. sfxA, sfxB, sfxC, sfxD'), all types share a conserved core sequence (~13 bp in R64) but differ in their flanking bases. Cross-type pairings like sfxA–sfxD' may fall below einverted's detection threshold (~94.7% identity for 19 bp sites differing at one position).

The motif refinement step addresses this by extracting a core motif from the IR arms already detected in each cluster, then scanning the surrounding genomic sequence for additional sites matching that core. Extended pairs — einverted artifacts where a single detected arm bridges across two adjacent sfx sites — are identified when both arms exceed 1.4× the cluster's median arm length, removed from the IR table, and their constituent individual sites are recovered by the motif scan instead. Newly discovered sites are stored as unpaired recognition sites with their partner coordinates set to NA.


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
│   ├── step_phava.py            # einverted IR detection + motif refinement
│   ├── step_ir_cds.py           # IR-CDS containment annotation
│   ├── step_gff.py              # GFF generation, merging, window extraction
│   ├── step_clinker.py          # Gene-organisation plot generation
│   ├── step_kofamscan.py        # KOfamscan on IR-overlapping CDS
│   └── hmms/                    # Bundled HMM profiles (.hmm.gz)
├── pyproject.toml               # Package metadata + console_scripts
├── environment.yml              # Conda + pip environment
├── README.md
└── .gitignore
```


## License

MIT
