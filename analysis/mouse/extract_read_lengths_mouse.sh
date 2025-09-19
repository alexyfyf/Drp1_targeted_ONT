#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Alex Yan
#              yan.a@wehi.edu.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=lentgh

# To set a project account for credit charging,
# SBATCH --account=ls25

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8

# Memory usage (MB)
#SBATCH --mem-per-cpu=4000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=8:00:00

# To receive an email when job completes or fails
#SBATCH --mail-user=yan.a@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN

#SBATCH --array=1-28

# Set the file for output (stdout)
#SBATCH --output=./length_dist_mouse/slurm-%A_%a.out

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

J=${SLURM_ARRAY_TASK_ID}

mkdir -p length_dist_mouse

raw_reads=$(sed "${J}q;d" fq_mouse.txt)
# base=$(basename $raw_reads _merged.fastq.gz)
if [[ "$raw_reads" =~ ^.*\/([0-9]{8})_.*\/demux\/.*_(barcode[0-9]+)\.fastq\.gz$ ]]; then
    date_part="${BASH_REMATCH[1]}"
    barcode_part="${BASH_REMATCH[2]}"
    base="${date_part}_${barcode_part}"
    echo "Original file path: $raw_reads"
    echo "Extracted name: $base"
else
    echo "Error: Pattern not matched for: $raw_reads"
fi

output_file="length_dist_mouse/${base}_length_distribution.txt"

seqkit fx2tab -l -g -n -i -H -j 8 $raw_reads > $output_file