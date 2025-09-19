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
#SBATCH --cpus-per-task=16

# Memory usage (MB)
#SBATCH --mem-per-cpu=4000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=8:00:00

# To receive an email when job completes or fails
#SBATCH --mail-user=yan.a@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN

#SBATCH --array=1-56

# Set the file for output (stdout)
# SBATCH --output=stdout

# Set the file for error log (stderr)
# SBATCH --error=stderr

# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name

# SBATCH --partition=genomics
# SBATCH --qos=genomics

# SBATCH --dependency=afterok:10567760
#SBATCH --output=/vast/scratch/users/yan.a/vast_scratch/tx_bam_mouse/%x_%A_%a.out
#SBATCH --error=/vast/scratch/users/yan.a/vast_scratch/tx_bam_mouse/%x_%A_%a.err
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
# base='Drp1_human'
outdir="/vast/scratch/users/yan.a/vast_scratch/tx_bam_mouse"

tx='/vast/projects/lab_davidson/yan.a/ref/ncbi/mouse/GCF_000001635.27_GRCm39_rna.fna.gz' 

# Create minimap2 index if it doesn't exist
# minimap2 -x map-ont -I 1000G -t $num_cores -d $outdir/tx.mmi $tx

tx_index="$outdir/tx.mmi"
input_folder="/vast/projects/davidson_longread/yan.a/20240326_stVincent_Max_MitocDNAPool/ncbi/dedup_bambu/Drp1_mouse_dge/"

J=${SLURM_ARRAY_TASK_ID}

raw_reads=$(sed "${J}q;d" /vast/scratch/users/yan.a/vast_scratch/tx_bam_mouse/bam.txt)

name=$(basename "$raw_reads" .bam)
out_fastq="$outdir/${name}.fastq"
out_bam="${outdir}/${name}_sorted.bam"

# Filter to NC_000082.7 region and convert to FASTQ in one pipe
samtools view -@ $num_cores -b "$raw_reads" NC_000082.7:16129000-16179000 | \
    samtools fastq -@ $num_cores -0 "$out_fastq" -s /dev/null -

# Align using the index
minimap2 -ax map-ont -Y -p 1.0 -N 100 -I 1000G -t $num_cores $tx_index "$out_fastq" | \
    samtools view -@ $num_cores -Sb | \
    samtools sort -@ $num_cores > "$out_bam"

# Index the output BAM
samtools index -@ $num_cores "$out_bam"

# Remove intermediate FASTQ
rm -f "$out_fastq"





