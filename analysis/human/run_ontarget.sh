#!/bin/bash
# Script to calculate read statistics for ON-target analysis
# Calculates: 1) Total reads from fastq, 2) Aligned reads from bam, 3) Alignments in specific region

# Usage: 
#   Single barcode: ./run_ontarget.sh [barcode_number] [region] [bam_file]
#   All barcodes:   ./run_ontarget.sh all [region] [bam_suffix]
# Example: 
#   ./run_ontarget.sh 01 "chr12:32679301-32745650" barcode01_sorted.bam
#   ./run_ontarget.sh all "chr12:32679301-32745650" sorted.bam

set -e  # Exit on any error

# Function to analyze a single barcode
analyze_barcode() {
    local BARCODE_NUM=$1
    local REGION=$2
    local BAM_FILE=$3
    
    # Set up paths
    local FQ_FILE="/vast/projects/davidson_longread/yan.a/Max_2nd/fastq/merged/hs/barcode${BARCODE_NUM}_merged/barcode${BARCODE_NUM}_merged_fp.fastq"
    local BAM_PATH="./Drp1_human_dge/${BAM_FILE}"
    
    # Check if files exist
    if [ ! -f "${FQ_FILE}" ]; then
        echo "ERROR: FastQ file not found: ${FQ_FILE}"
        return 1
    fi
    
    if [ ! -f "${BAM_PATH}" ]; then
        echo "ERROR: BAM file not found: ${BAM_PATH}"
        return 1
    fi
    
    # Count total reads in fastq file
    local TOTAL_READS=$(wc -l < "${FQ_FILE}")
    TOTAL_READS=$((TOTAL_READS / 4))
    
    # Count primary mapped alignments
    local PRIMARY_ALIGNMENTS=$(samtools view -c -F 256 -F 2048 -F 4 "${BAM_PATH}")
    
    # Count primary mapped alignments in region
    local REGION_ALIGNMENTS=$(samtools view -c -F 256 -F 2048 -F 4 "${BAM_PATH}" ${REGION})
    
    # Calculate rates
    local PRIMARY_ALIGNMENT_RATE="0.00"
    local ON_TARGET_RATE="0.00"
    
    if [ ${TOTAL_READS} -gt 0 ]; then
        PRIMARY_ALIGNMENT_RATE=$(echo "scale=2; ${PRIMARY_ALIGNMENTS} * 100 / ${TOTAL_READS}" | bc -l)
    fi
    
    if [ ${PRIMARY_ALIGNMENTS} -gt 0 ]; then
        ON_TARGET_RATE=$(echo "scale=2; ${REGION_ALIGNMENTS} * 100 / ${PRIMARY_ALIGNMENTS}" | bc -l)
    fi
    
    # Output results
    echo "${BARCODE_NUM},${TOTAL_READS},${PRIMARY_ALIGNMENTS},${REGION_ALIGNMENTS},${PRIMARY_ALIGNMENT_RATE},${ON_TARGET_RATE}"
}

# Function to run analysis for all barcodes
run_all_barcodes() {
    local REGION=$1
    local BAM_SUFFIX=$2
    local OUTPUT_FILE="ontarget_analysis_$(date +%Y%m%d_%H%M%S).csv"
    
    echo "Running ON-target analysis for all barcodes..."
    echo "Region: ${REGION}"
    echo "BAM suffix: ${BAM_SUFFIX}"
    echo "Output file: ${OUTPUT_FILE}"
    echo ""
    
    # Create header
    echo "Barcode,Total_Reads,Primary_Mapped_Alignments,Region_Alignments,Primary_Alignment_Rate(%),On_Target_Rate(%)" > ${OUTPUT_FILE}
    
    # Run analysis for each barcode
    for i in {01..11}; do
        local BAM_FILE="barcode${i}_${BAM_SUFFIX}"
        echo "Processing barcode ${i}..."
        
        if analyze_barcode ${i} "${REGION}" "${BAM_FILE}" >> ${OUTPUT_FILE} 2>/dev/null; then
            echo "  ✓ Completed barcode ${i}"
        else
            echo "  ✗ Failed barcode ${i}"
            echo "${i},ERROR,ERROR,ERROR,ERROR,ERROR" >> ${OUTPUT_FILE}
        fi
    done
    
    echo ""
    echo "Analysis complete! Results saved to: ${OUTPUT_FILE}"
    echo ""
    echo "Summary table:"
    cat ${OUTPUT_FILE}
}

# Check if required arguments are provided
if [ $# -lt 1 ]; then
    echo "Usage:"
    echo "  Single barcode: $0 <barcode_number> [region] [bam_file]"
    echo "  All barcodes:   $0 all [region] [bam_suffix]"
    echo ""
    echo "Examples:"
    echo "  $0 01 \"chr12:32660000-32800000\" barcode01_sorted.bam"
    echo "  $0 all \"chr12:32660000-32800000\" sorted.bam"
    echo ""
    echo "Arguments:"
    echo "  barcode_number: Barcode number (01-11) or 'all' for all barcodes"
    echo "  region: Optional genomic region in format 'chr:start-end' (default: chr12)"
    echo "  bam_file/suffix: BAM file name or suffix (default: barcode{number}_sorted.bam)"
    exit 1
fi

# Load required modules
module purge
module load micromamba/
eval "$(micromamba shell hook --shell bash)"
micromamba activate /vast/projects/lab_davidson/yan.a/conda_env/isoncorrect
module load samtools/1.19.2

# Check if running analysis for all barcodes
if [ "$1" = "all" ]; then
    REGION=${2:-"chr12:32679301-32745650"}
    BAM_SUFFIX=${3:-"sorted.bam"}
    run_all_barcodes "${REGION}" "${BAM_SUFFIX}"
else
    # Single barcode analysis
    BARCODE_NUM=$1
    REGION=${2:-"chr12:32679301-32745650"}
    BAM_FILE=${3:-"barcode${BARCODE_NUM}_sorted.bam"}
    
    # Set up paths
    FQ_FILE="/vast/projects/davidson_longread/yan.a/Max_2nd/fastq/merged/hs/barcode${BARCODE_NUM}_merged/barcode${BARCODE_NUM}_merged_fp.fastq"
    BAM_PATH="./Drp1_human_dge/${BAM_FILE}"
    
    echo "=========================================="
    echo "ON-Target Analysis for Barcode ${BARCODE_NUM}"
    echo "=========================================="
    echo "FastQ file: ${FQ_FILE}"
    echo "BAM file: ${BAM_PATH}"
    echo "Region: ${REGION}"
    echo "=========================================="
    
    # Check if files exist
    if [ ! -f "${FQ_FILE}" ]; then
        echo "ERROR: FastQ file not found: ${FQ_FILE}"
        exit 1
    fi
    
    if [ ! -f "${BAM_PATH}" ]; then
        echo "ERROR: BAM file not found: ${BAM_PATH}"
        exit 1
    fi
    
    echo ""
    echo "1. Counting total reads from FastQ file..."
    echo "----------------------------------------"
    
    # Count total reads in fastq file (divide by 4 since each read has 4 lines)
    TOTAL_READS=$(wc -l < "${FQ_FILE}")
    TOTAL_READS=$((TOTAL_READS / 4))
    echo "Total reads in FastQ: ${TOTAL_READS}"
    
    echo ""
    echo "2. Counting aligned reads from BAM file..."
    echo "----------------------------------------"
    
    # Count primary mapped alignments (exclude secondary, supplementary, and unmapped)
    PRIMARY_ALIGNMENTS=$(samtools view -c -F 256 -F 2048 -F 4 "${BAM_PATH}")
    echo "Primary mapped alignments: ${PRIMARY_ALIGNMENTS}"
    
    echo ""
    echo "3. Counting alignments in specific region..."
    echo "----------------------------------------"
    
    # Count primary mapped alignments in the specified region
    REGION_ALIGNMENTS=$(samtools view -c -F 256 -F 2048 -F 4 "${BAM_PATH}" ${REGION})
    echo "Primary mapped alignments in region ${REGION}: ${REGION_ALIGNMENTS}"
    
    echo ""
    echo "=========================================="
    echo "SUMMARY"
    echo "=========================================="
    echo "1. Total reads from FastQ: ${TOTAL_READS}"
    echo "2. Primary mapped alignments in BAM: ${PRIMARY_ALIGNMENTS}"
    echo "3. Primary mapped alignments in region ${REGION}: ${REGION_ALIGNMENTS}"
    
    # Calculate percentages
    if [ ${TOTAL_READS} -gt 0 ]; then
        PRIMARY_ALIGNMENT_RATE=$(echo "scale=2; ${PRIMARY_ALIGNMENTS} * 100 / ${TOTAL_READS}" | bc -l)
        echo ""
        echo "Primary mapped alignment rate: ${PRIMARY_ALIGNMENT_RATE}%"
    fi
    
    # Calculate on-target rate (alignments in region / total alignments)
    if [ ${PRIMARY_ALIGNMENTS} -gt 0 ]; then
        ON_TARGET_RATE=$(echo "scale=2; ${REGION_ALIGNMENTS} * 100 / ${PRIMARY_ALIGNMENTS}" | bc -l)
        echo "On-target rate (region ${REGION}): ${ON_TARGET_RATE}%"
    fi
    
    echo "=========================================="
fi