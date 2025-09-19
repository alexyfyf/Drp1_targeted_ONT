#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Alex Yan
#              yan.a@wehi.edu.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=rysnc

# To set a project account for credit charging,
# SBATCH --account=ls25

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8

# Memory usage (MB)
#SBATCH --mem-per-cpu=2000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=8:00:00

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

J=${SLURM_ARRAY_TASK_ID}

TSO_REV_CMP="TATCAACGCAGAGTACA" # partial sequence, entire TSO AAGCAGTGGTATCAACGCAGAGTACA
UMI_REV_CMP="???YR???YR???YR???V" 
UMI_LEFT_FLANK_REV_CMP="TGGG" 

BC_LEFT_FLANK_FWD="AAGGTTAA"
BC_RIGHT_FLANK_FWD="CAGCACCT"

# Define the exon primer variable using if-else statements
# if [ "$J" -le 8 ]; then
    echo 'Human sample'
    SP="hs"
    EXON_FWD="AAGATGAGTCTCCCGGATTTCAGC"  # human
# elif [ "$J" -gt 8 ]; then
#     echo 'mouse sample'
#     SP="mm"
#     EXON_FWD="AAGATGAGTCTCTCGGATTTCAGCA"  # mouse
# fi

NANOPORE_PRIMER_FWD="TTCGTTCAGTTACGTATTGCT" # partial sequence, entire Naopore primer TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT

BC_SEQ="????????????????????????" # 24 bp

raw_reads=$(sed "${J}q;d" fq.txt)

base=$(basename $raw_reads .fastq.gz)

mkdir -p ${SP}/${base}
cd ${SP}/${base}

# -b "" -k "?" will extract the UMI only, as per the documentation:
#   https://davidsongroup.github.io/flexiplex/#extracting-umis-from-pcr-cdna-ont-data
zcat ${raw_reads} | \
flexiplex -x $TSO_REV_CMP -u $UMI_REV_CMP -x $UMI_LEFT_FLANK_REV_CMP \
  -b "" -k "?" \
  -f 2 -e 1 -p 8 -n "UMI_extraction" | \
flexiplex -x $BC_LEFT_FLANK_FWD -b $BC_SEQ -x $BC_RIGHT_FLANK_FWD -k "?" -f 1 -p 8 -n "BC_extraction"  | \
sed 's/\(\?_#?\)\(_[^#]*\)\(#.*\)/\1\3\2/' > ${base}_fp.fastq
