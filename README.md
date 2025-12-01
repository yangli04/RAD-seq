# RAD-seq

This repository contains a Snakemake pipeline for preprocessing RAD-seq data: quality control, adapter trimming, rRNA filtering, genome mapping, deduplication, and feature counting.

**Overview**
- **Workflow:** Implemented in `Snakefile`.
- **Configuration:** Provided by `config.yaml` (edit paths there to match your system).
- **Input:** Paired-end FASTQ files expected in `raw_data/` named like `<SAMPLE>_S1_R1_001.fastq.gz` and `<SAMPLE>_S1_R2_001.fastq.gz`.

**Required Software**
The pipeline invokes a number of command-line tools. Install these (preferably via `conda`/`mamba` from `bioconda` and `conda-forge`) and make them available in your `PATH`, or set their absolute paths in `config.yaml`.
- `snakemake` (to run the workflow)
- `hisat-3n` and `hisat-3n-table` (genome/rRNA mapping and table generation)
- `samtools` (BAM processing, indexing, flagstat)
- `umi_tools` (deduplication)
- `cutadapt` (adapter trimming)
- `featureCounts` (from the `subread` package)
- `multiqc` (aggregate QC reports)
- `falco` (used here for FastQC-like reports; the pipeline expects a `falco` binary path in `config.yaml`) 
- `variant motif` (used by the `calculate_table_*` rules; ensure `variant` and the `motif` subcommand are available)

Recommended (example) conda/mamba environment creation:

```bash
# using mamba (faster) if available
mamba create -n radseq -c conda-forge -c bioconda \
	snakemake samtools cutadapt umi_tools subread multiqc hisat-3n hisat3n-table \
	python=3.10 -y
conda activate radseq

# Install falco or other QC tools if not available via conda. Adjust config.yaml to point to the binary.
```

If some packages are unavailable under the given names in your channels, install them from their project pages or place their binaries somewhere on `PATH` and update `config.yaml` accordingly.

**Configuration**
- Edit `config.yaml` to set the paths to binaries and reference files the workflow needs. Keys used in `Snakefile` include:
	- `falco`, `umi_tools`, `cutadapt`, `featureCounts`
	- `ref_genome`, `ref_genome_fa`, `ref_rRNA`, `ref_rRNA_fa`, `chrom_sizes`, `ref_hg_genome_gtf`

Note: `Snakefile` contains `configfile: "workflow/config.yaml"` — if your `config.yaml` is at the repository root, run Snakemake with `--configfile config.yaml` (example below) or move/copy your config to `workflow/config.yaml`.

**How to run**
1. Do a dry-run to check what will run and that inputs exist:

```bash
snakemake -n -s Snakefile --configfile config.yaml
```

2. Run the full workflow using multiple cores (example 16 cores):

```bash
snakemake -s Snakefile --configfile config.yaml --cores 16
```

3. To execute a single target rule (e.g., mapping to genome tables):

```bash
snakemake INTERDIR/hisat3n_table/genome/<SAMPLE>.hisat3n_table.tsv -s Snakefile --configfile config.yaml --cores 8
```

**Notes & Troubleshooting**
- Ensure FASTQ filenames match the sample names listed in `Snakefile` or update the `SAMPLE` list in `Snakefile`.
- If a binary path is set in `config.yaml`, Snakemake will use that path (e.g., `umi_tools` in config). Update these paths when moving the pipeline between systems.
- Logs and intermediate files are written under `preprocess/interdir` and `preprocess/outdir` per the `Snakefile` constants.
- If `hisat-3n` or `hisat-3n-table` are not found, either install them and add to `PATH`, or set their paths in `config.yaml` and adjust the `Snakefile` if necessary.

**Outputs**
- QC reports: `report_qc/` (MultiQC reports)
- Counts: `preprocess/outdir/feature_counts/counts_hg_genome.txt`
- Calculated pileup tables: `preprocess/outdir/calculated_table/{genome|rRNA}/`

**Contact / Changes**
- If you want help adapting this pipeline for new samples or a different reference, open an issue or contact the author.

