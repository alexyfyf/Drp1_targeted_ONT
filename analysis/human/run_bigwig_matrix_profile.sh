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
#SBATCH --time=48:00:00

# To receive an email when job completes or fails
#SBATCH --mail-user=yan.a@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN

# SBATCH --array=1-44

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

# Load deeptools
module purge
module load deeptools

# Directories
BAMDIR="Drp1_human_dge"
OUTDIR="plot_coverage"
# GTF="ref/DNM1L.gtf"
GTF="ref/NM_001278464.gtf" # this one contains all exons

mkdir -p "$OUTDIR"

# Step 1: Generate bigwig files from BAMs
for BAM in "$BAMDIR"/*.bam; do
    SAMPLE=$(basename "$BAM" .bam)
    BIGWIG="$OUTDIR/${SAMPLE}.bw"
    if [[ ! -f "$BIGWIG" ]]; then
        echo "Generating bigwig for $BAM..."
        bamCoverage -b "$BAM" -o "$BIGWIG" --normalizeUsing None --binSize 1 -p 16 --region chr12
    else
        echo "$BIGWIG already exists, skipping."
    fi
done

# Step 2: Compute matrix
computeMatrix scale-regions \
    -S "$OUTDIR"/*dedup.bw \
    -R "$GTF" \
    --metagene \
    --beforeRegionStartLength 0 \
    --regionBodyLength 3000 \
    --afterRegionStartLength 0 \
    --skipZeros \
    --binSize 50 \
    -o "$OUTDIR/coverage_matrix_dedup.gz" \
    --outFileNameMatrix "$OUTDIR/coverage_matrix_dedup.tab" \
    -p 16

# Step 3: Plot profile
plotProfile \
    -m "$OUTDIR/coverage_matrix_dedup.gz" \
    -out "$OUTDIR/profile_plot_dedup.pdf" \
    --perGroup \
    --plotType lines \
    --legendLocation best

# Step 2: Compute matrix
computeMatrix scale-regions \
    -S "$OUTDIR"/*sorted.bw \
    -R "$GTF" \
    --metagene \
    --beforeRegionStartLength 0 \
    --regionBodyLength 3000 \
    --afterRegionStartLength 0 \
    --skipZeros \
    --binSize 50 \
    -o "$OUTDIR/coverage_matrix_sorted.gz" \
    --outFileNameMatrix "$OUTDIR/coverage_matrix_sorted.tab" \
    -p 16

# Step 3: Plot profile
plotProfile \
    -m "$OUTDIR/coverage_matrix_sorted.gz" \
    -out "$OUTDIR/profile_plot_sorted.pdf" \
    --perGroup \
    --plotType lines \
    --legendLocation best

echo "All done. Results in $OUTDIR" 