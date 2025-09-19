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
#SBATCH --mem-per-cpu=50000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=48:00:00

# To receive an email when job completes or fails
#SBATCH --mail-user=yan.a@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN

# SBATCH --array=1-11

# Set the file for output (stdout)
# SBATCH --output=stdout

# Set the file for error log (stderr)
# SBATCH --error=stderr

# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name

# SBATCH --partition=long
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

# Activate umi_tools environment for deduplication
micromamba activate /vast/projects/lab_davidson/yan.a/conda_env/umi_tools_env 

# --- User parameters ---
chr="chr12"                      # chromosome to analyse
outdir="saturation_results"   # output folder
# Fractions of total reads to sample (e.g., 0.1 = 10%, 1.0 = 100%)
fractions=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
# fractions=(0.7 0.8 0.9 1.0)
# fractions=(0.9)

# Find all sorted BAM files in Drp1_human_dge/ folder
bam_files=(Drp1_human_dge/*sorted.bam)

# rerun barcode11 due to OOM error
# bam_files=(Drp1_human_dge/barcode11_sorted.bam)

# Check if any BAM files were found
if [ ${#bam_files[@]} -eq 0 ]; then
    echo "Error: No *sorted.bam files found in Drp1_human_dge/ folder"
    exit 1
fi

echo "Found ${#bam_files[@]} BAM file(s): ${bam_files[*]}"

# Create output directory
mkdir -p "${outdir}"

# Process each BAM file
for bam in "${bam_files[@]}"; do
    # Extract base name from BAM file
    name=$(basename "${bam}" _sorted.bam)
    echo "Processing BAM file: ${bam}"
    
    # Header for results CSV for this BAM file
    csv_file="${outdir}/${name}_saturation_curve.csv"
    echo "fraction,total_reads,unique_umi_reads" > "${csv_file}"
    
    # Loop over fractions
    for frac in "${fractions[@]}"; do
      echo "Processing fraction ${frac} for ${name} …"

      # 1. Subsample the BAM on chr12 at fraction=frac
      #    -s takes a FLOAT in [0,1] to randomly subsample reads
      subsampled_bam="${outdir}/${name}_sub_${frac//./}.bam"
      samtools view -b -s "${frac}" -h "${bam}" "${chr}" \
        -o "${subsampled_bam}"

      # Index the subsampled BAM file
      echo "Indexing subsampled BAM file..."
      samtools index "${subsampled_bam}"

      # 2. Deduplicate with umi_tools
      stats_prefix="${outdir}/${name}_stats_${frac//./}"
      dedup_bam="${outdir}/${name}_dedup_${frac//./}.bam"
      
      echo "Running umi_tools dedup on ${subsampled_bam}..."
      umi_tools dedup \
        --method=directional \
        --chrom="${chr}" \
        --edit-distance-threshold=2 \
        --output-stats="${stats_prefix}" \
        -I "${subsampled_bam}" \
        -S "${dedup_bam}"
      
      # Check if umi_tools succeeded
      if [ ! -f "${dedup_bam}" ] || [ ! -s "${dedup_bam}" ]; then
        echo "ERROR: umi_tools dedup failed or produced empty output for ${dedup_bam}"
        unique_umi_reads="ERROR"
      else
        echo "umi_tools dedup completed successfully"
      fi

      # 3. Count total reads in subsampled BAM (pre‐dedup)
      total_reads=$(samtools view -c "${subsampled_bam}")
      echo "Total reads in subsampled BAM: ${total_reads}"

      # 4. Count unique UMI reads (post‐dedup)
      if [ "${unique_umi_reads}" != "ERROR" ]; then
        unique_umi_reads=$(samtools view -c "${dedup_bam}")
        echo "Unique UMI reads: ${unique_umi_reads}"
      fi

      # 5. Append to CSV
      echo "${frac},${total_reads},${unique_umi_reads}" \
        >> "${csv_file}"
        
      # Clean up intermediate files
      rm -f "${subsampled_bam}" "${dedup_bam}"
    done
    
    echo "Done! Results for ${name} in ${csv_file}"
done

echo "All BAM files processed!"
