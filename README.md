# Drp1_targeted_ONT

Targeted ONT sequencing analysis for the Drp1 locus (human/mouse).

- Drp1_targeted_ONT
  - 1, Preprocess targeted ONT reads
  - 2, Isoform and transcript analysis
  - 3, Results and plot generation
  - Citation and Data availability

### 1, Preprocess targeted ONT reads

0. Basecalling/demultiplexing as required by your sequencing run (Dorado, SUP mode without trimming). Place resulting `*.fastq[.gz]` in a raw data location (not tracked).
1. Use `flexiplex` to extract UMI, remove barcode and reorder read header for downstream analysis.
2. Alignment with `minimap2` per species and sample. Use `umi-tools` to deduplicate reads based on UMI, and restricted to Drp1-containing chromosome.
3. (Optional) Merge bam files from 2 runs.
4. Run `bambu` to quantify transcript expression.
5. Filter only Drp1 related transcripts.

### 2, Isoform and transcript analysis

1. Extract raw read length and plot the results.
2. Extract on target rate and plot the results.
3. Generate genomic coverage matrix and plot the results.
4. Re-align to transcriptome and plot transcript coverage.
5. Perform saturation analysis by downsampling the bam files and plot the results.

### 3, Rendered html files

1. Render the qmd file to generate html file.

### Repository layout

```
Drp1_targeted_ONT/
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

### Citation

GEO
PMID
