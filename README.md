# Drp1_targeted_ONT

Targeted ONT sequencing analysis for the Drp1 locus (human/mouse).

- Drp1_targeted_ONT
  - 1, Preprocess targeted ONT reads
  - 2, Isoform and transcript analysis
  - 3, Variant/structural variant and fusion analysis
  - 4, Generate genome browser tracks
  - Citation

This README mirrors the concise, numbered style of the Ki67 ATAC analysis repo for clarity and reproducibility ([Bcell_Ki67_ATAC_JEM](https://github.com/alexyfyf/Bcell_Ki67_ATAC_JEM)).

### 1, Preprocess targeted ONT reads

1. Basecalling/demultiplexing as required by your sequencing run (e.g. Guppy, Dorado). Place resulting `*.fastq[.gz]` in a raw data location (not tracked).
2. Filtering/adapter trimming if needed (e.g. `NanoFilt`, `porechop`).
3. Alignment with minimap2 per species and sample. Example:
   ```bash
   # inside preprocess/human or preprocess/mouse
   REF=reference.fa
   SAMPLE=sample1
   minimap2 -t 16 -ax splice -uf --secondary=no ${REF} ${SAMPLE}.fastq.gz | \
     samtools view -bS - | samtools sort -o ${SAMPLE}.sorted.bam
   samtools index ${SAMPLE}.sorted.bam
   ```
4. QC summaries (e.g. `samtools flagstat`, `seqkit stats`, `pycoQC` as applicable). Record software versions in a `pipeline_info` note.

Environment (micromamba recommended) [[memory:3407552]]:
```bash
micromamba create -n drp1_ont -c conda-forge -c bioconda \
  python=3.11 samtools minimap2 bedtools seqkit \
  pandas matplotlib 
micromamba activate drp1_ont
```

### 2, Isoform and transcript analysis

1. Generate isoform-level annotations from aligned reads (choose one):
   - FLAIR or TALON for long-read isoform calling/annotation
   - StringTie2 (long-read mode) for transcript assembly
2. Summarize Drp1 locus isoforms and expression by species.
3. Export tables and plots to `results/<species>/isoform/`.

### 3, Variant/structural variant and fusion analysis

1. Small variants: `medaka` for ONT polishing or `longshot` for variant calling near Drp1.
2. Structural variants: `sniffles` or `cuteSV` from BAMs.
3. Gene fusions (if applicable): long-read fusion callers (e.g. `JAFFA` long-read mode).
4. Save summary tables/VCFs under `results/<species>/variants/` and figures under `results/<species>/figures/`.

### 4, Generate genome browser tracks

1. Produce coverage bigWig per sample and meta-sample:
   ```bash
   bamCoverage -b ${SAMPLE}.sorted.bam -o ${SAMPLE}.bw -bs 10 --normalizeUsing CPM
   ```
2. Organize tracks by species in `results/<species>/tracks/`.
3. Optionally build a UCSC track hub similar to the style in the Ki67 repo.

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

Add project-specific citation(s) here. Cite ONT tools and any downstream packages used. The numbered section style is adapted from the Ki67 ATAC analysis repo ([Bcell_Ki67_ATAC_JEM](https://github.com/alexyfyf/Bcell_Ki67_ATAC_JEM)).
