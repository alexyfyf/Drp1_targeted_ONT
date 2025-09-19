# Drp1_targeted_ONT

Oxford Nanopore targeted sequencing workflow focused on the Drp1 (Dynamin-related protein 1) locus and related regions. This repo organizes preprocessing, analysis, and results for human and mouse datasets.

## Overview

- Goal: characterize Drp1 using targeted ONT reads, including preprocessing (basecalled reads in, cleaned/demultiplexed/aligned out) and downstream analyses (coverage, isoforms, variants, and QC).
- Species: human and mouse supported in parallel folder structures.

## Repository Structure

```
Drp1_targeted_ONT/
├── README.md
├── preprocess/
│   ├── human/
│   └── mouse/
├── analysis/
│   ├── human/
│   └── mouse/
└── results/
    ├── human/
    └── mouse/
```

- `preprocess/`: scripts and configs for adapter trimming, demultiplexing, filtering, alignment to reference, and basic QC summaries per sample.
- `analysis/`: notebooks and scripts for downstream interpretation (coverage plots, isoform/feature analysis, variant summaries, figures).
- `results/`: generated outputs and figures organized by species and analysis step. Large intermediates should not be committed.

## Environment

Example micromamba environment (preferred) [[memory:3407552]]:

```bash
micromamba create -n drp1_ont -c conda-forge -c bioconda \
  python=3.11 samtools minimap2 bedtools seqkit pandas matplotlib 
micromamba activate drp1_ont
```

Optional tools you may add depending on your workflow: `guppy`, `qcat`/`porechop`, `nanofilt`, `pycoQC`, `fastp`, `htslib`, `r-base`, `bioconductor-*`, `jupyterlab`.

## Data Expectations

- Inputs: raw/basecalled ONT reads (`*.fastq[.gz]`) per sample, and reference genomes/transcripts appropriate to species.
- Outputs: sorted/indexed BAMs, per-sample QC metrics, coverage tracks, isoform/variant tables, and figures saved under `results/<species>/`.

## Typical Workflow

1) Preprocess (per species)

```bash
# inside preprocess/human or preprocess/mouse
# 1. Trim/adapters (if needed)
# 2. Filter/Demultiplex
# 3. Align to reference with minimap2
# 4. Sort/index BAM, compute basic QC
```

2) Analysis (per species)

```bash
# inside analysis/human or analysis/mouse
# - Run notebooks/scripts to generate coverage plots, isoform summaries, variant calls, and figures
# - Save outputs under results/<species>/
```

3) Results

```bash
# figures, tables, and summary reports live under results/human or results/mouse
```

## How to Reproduce

1. Prepare environment (see Environment).
2. Place input FASTQs under an appropriate raw data location (not tracked by git).
3. Run preprocessing for the desired species in `preprocess/<species>/`.
4. Execute analysis notebooks/scripts in `analysis/<species>/`.
5. Find outputs in `results/<species>/`.

## Versioning & Remotes

Set the GitHub remote if not already done:

```bash
git remote set-url origin https://github.com/<your-username>/Drp1_targeted_ONT.git
git push -u origin master
```

## License

TBD.

## Citation

If you use this repository or derived outputs, please cite the appropriate ONT tools and any downstream packages used. Add project-specific citations here when available.
