#!/usr/bin/env Rscript

library(bambu)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

bam.files <- list.files(args[1], pattern = '_dedup.bam$', full.names = TRUE)

bam.files <- normalizePath(bam.files)
length(bam.files)

if (grepl('human', args[1])) {
    gtf.file <- "/home/users/allstaff/yan.a/vast_scratch/Max_2nd/fastq/merged/ncbi/dedup_bambu/ref/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"
    fa.file <- "/home/users/allstaff/yan.a/vast_scratch/Max_2nd/fastq/merged/ncbi/dedup_bambu/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
} else if (grepl('mouse', args[1])) {
    gtf.file <- "/home/users/allstaff/yan.a/davidson_longread/yan.a/20240326_stVincent_Max_MitocDNAPool/ncbi/dedup/temp_ref/GCF_000001635.27_GRCm39_genomic.gtf"
    fa.file <- "/home/users/allstaff/yan.a/davidson_longread/yan.a/20240326_stVincent_Max_MitocDNAPool/ncbi/dedup/temp_ref/GCF_000001635.27_GRCm39_genomic.fna"
}

bambuAnnotations <- prepareAnnotations(gtf.file)

## run default
dir.create(paste0(args[1],'/default'))

se1 <- bambu(reads = bam.files, 
             annotations = bambuAnnotations, 
             stranded = F, NDR = 0.2,
             genome = fa.file, ncore = 32)

data.frame(mcols(se1)) %>% count(novelTranscript, is.na(readCount))

writeBambuOutput(se1, path = paste0(args[1],'/default'))

saveRDS(se1, paste0(args[1],'/default/se1.rds'))



