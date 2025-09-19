#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Alex Yan
#              yan.a@wehi.edu.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=bambu

# To set a project account for credit charging,
# SBATCH --account=ls25

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12

# Memory usage (MB)
#SBATCH --mem-per-cpu=10000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=7:00:00

# To receive an email when job completes or fails
#SBATCH --mail-user=yan.a@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN

#SBATCH --array=1-11

# Set the file for output (stdout)
# SBATCH --output=stdout

# Set the file for error log (stderr)
# SBATCH --error=stderr

# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name

# SBATCH --partition=genomics
# SBATCH --qos=genomics

# SBATCH --dependency=afterok:10567760
# SBATCH --reservation=highmem

# Command to run a serial job

module purge
module load micromamba/

eval "$(micromamba shell hook --shell bash)"
micromamba activate /vast/projects/lab_davidson/yan.a/conda_env/isoncorrect

minimap2 --version
# 2.26-r1175
module load salmon/1.10.2
salmon --version

module load samtools/1.19.2
samtools --version 

## use raw fasta from output, not SQANTI corrected
num_cores=12
base='Drp1_human'

gtf="/vast/projects/lab_davidson/yan.a/ref/ncbi/human_no_alt/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz"
genome="/vast/projects/lab_davidson/yan.a/ref/ncbi/human_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
chr="chr12"

mkdir -p ${base}_dge
cd ${base}_dge

J=${SLURM_ARRAY_TASK_ID}

raw_reads=$(sed "${J}q;d" ../fq_human.txt)
name=$(basename ${raw_reads} .fastq | awk -F'_' '{print $1}')

paftools.js gff2bed $gtf > ${name}_anno.bed

## for unstraned cDNA

minimap2 -ax splice -I 1000G -t $num_cores --junc-bed ${name}_anno.bed $genome $raw_reads | samtools sort -@ $num_cores -O BAM -o ${name}_sorted.bam

samtools index -@ $num_cores ${name}_sorted.bam

bedtools bamtobed -i ${name}_sorted.bam | sort -k2,2n -k3,3nr > ${name}.bed

# module load umi-tools/1.1.4
micromamba activate /vast/projects/lab_davidson/yan.a/conda_env/umi_tools_env
## umi_tools version 1.1.6

umi_tools dedup -I ${name}_sorted.bam --chrom=${chr} --output-stats=${name} --edit-distance-threshold=2 --method=directional -S ${name}_dedup.bam

cd ..

