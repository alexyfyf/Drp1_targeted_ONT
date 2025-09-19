#!/bin/bash
# Script to calculate read statistics for ON-target analysis
# Calculates: 1) Total reads from fastq, 2) Aligned reads from bam, 3) Alignments in specific region

# Usage: 
#   Single barcode: ./run_ontarget.sh [barcode_number] [region] [bam_suffix] [date]
#   All barcodes:   ./run_ontarget.sh all [region] [bam_suffix]
# Example: 
#   ./run_ontarget.sh 01 sorted.bam 20240326
#   ./run_ontarget.sh 15 sorted.bam 20240328
#   ./run_ontarget.sh all
# Default regions:
#   Human (barcodes 01-08): chr12:32679301-32745650
#   Mouse (barcodes 09-22): NC_000082.7:16130092-16176863

set -e  # Exit on any error

# Function to analyze a single barcode
analyze_barcode() {
    local BARCODE_NUM=$1
    local REGION=$2
    local BAM_FILE=$3
    local DATE=$4
    
    # Determine species and set up paths based on barcode number
    local SPECIES=""
    local FQ_FILE=""
    local BAM_PATH=""
    
    # Remove leading zeros for arithmetic operations
    local BARCODE_NUM_CLEAN=$((10#${BARCODE_NUM}))
    
    # Determine if this is a human or mouse barcode
    if [ ${BARCODE_NUM_CLEAN} -ge 1 ] && [ ${BARCODE_NUM_CLEAN} -le 8 ]; then
        SPECIES="human"
        # Get FQ file path from fq_human.txt based on barcode and date
        if [ "${DATE}" = "20240326" ]; then
            local FQ_LINE_NUM=$((BARCODE_NUM_CLEAN))
        else
            local FQ_LINE_NUM=$((BARCODE_NUM_CLEAN + 8))
        fi
        FQ_FILE=$(sed -n "${FQ_LINE_NUM}p" fq_human.txt)
        BAM_PATH="./Drp1_human_dge/${DATE}_barcode$(printf "%02d" ${BARCODE_NUM_CLEAN})_${BAM_FILE}"
    elif [ ${BARCODE_NUM_CLEAN} -ge 9 ] && [ ${BARCODE_NUM_CLEAN} -le 22 ]; then
        SPECIES="mouse"
        # Get FQ file path from fq_mouse.txt based on barcode and date
        if [ "${DATE}" = "20240326" ]; then
            local FQ_LINE_NUM=$((BARCODE_NUM_CLEAN - 8))
        else
            local FQ_LINE_NUM=$((BARCODE_NUM_CLEAN - 8 + 14))
        fi
        FQ_FILE=$(sed -n "${FQ_LINE_NUM}p" fq_mouse.txt)
        BAM_PATH="./Drp1_mouse_dge/${DATE}_barcode$(printf "%02d" ${BARCODE_NUM_CLEAN})_${BAM_FILE}"
    else
        echo "ERROR: Invalid barcode number ${BARCODE_NUM}. Must be 1-22."
        return 1
    fi
    
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
    local REGION_ALIGNMENTS=$(samtools view -c -F 256 -F 2048 -F 4 "${BAM_PATH}" "${REGION}")
    
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
    echo "${BARCODE_NUM},${DATE},${TOTAL_READS},${PRIMARY_ALIGNMENTS},${REGION_ALIGNMENTS},${PRIMARY_ALIGNMENT_RATE},${ON_TARGET_RATE}"
}

# Function to run analysis for all barcodes
run_all_barcodes() {
    local BAM_SUFFIX=${1:-"sorted.bam"}
    local OUTPUT_FILE="ontarget_analysis_$(date +%Y%m%d_%H%M%S).csv"
    
    echo "Running ON-target analysis for all barcodes..."
    echo "BAM suffix: ${BAM_SUFFIX}"
    echo "Using species-specific regions:"
    echo "  Human samples (barcodes 01-08): chr12:32679301-32745650"
    echo "  Mouse samples (barcodes 09-22): NC_000082.7:16130092-16176863"
    echo "Output file: ${OUTPUT_FILE}"
    echo ""
    
    # Create header
    echo "Barcode,Date,Total_Reads,Primary_Mapped_Alignments,Region_Alignments,Primary_Alignment_Rate(%),On_Target_Rate(%)" > ${OUTPUT_FILE}
    
    # Run analysis for each barcode (human: 01-08, mouse: 09-22)
    for i in {01..22}; do
        # Remove leading zeros for arithmetic operations
        local i_clean=$((10#${i}))
        
        echo "Processing barcode $(printf "%02d" ${i_clean})..."
        
        # Determine species-specific region (ignore REGION parameter for 'all' mode)
        local SPECIES_REGION=""
        if [ ${i_clean} -ge 1 ] && [ ${i_clean} -le 8 ]; then
            SPECIES_REGION="chr12:32679301-32745650"
        elif [ ${i_clean} -ge 9 ] && [ ${i_clean} -le 22 ]; then
            SPECIES_REGION="NC_000082.7:16130092-16176863"
        fi
        
        # Process both dates for each barcode
        for DATE in "20240326" "20240328"; do
            echo "  Processing ${DATE}..."
            
            if analyze_barcode ${i_clean} "${SPECIES_REGION}" "${BAM_SUFFIX}" "${DATE}" >> ${OUTPUT_FILE} 2>/dev/null; then
                echo "    ✓ Completed barcode $(printf "%02d" ${i_clean}) for ${DATE}"
            else
                echo "    ✗ Failed barcode $(printf "%02d" ${i_clean}) for ${DATE}"
                echo "$(printf "%02d" ${i_clean}),${DATE},ERROR,ERROR,ERROR,ERROR,ERROR" >> ${OUTPUT_FILE}
            fi
        done
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
    echo "  Single barcode: $0 <barcode_number> [region] [bam_suffix] [date]"
    echo "  All barcodes:   $0 all [bam_suffix]"
    echo ""
    echo "Examples:"
    echo "  $0 01 sorted.bam 20240326"
    echo "  $0 15 sorted.bam 20240328"
    echo "  $0 all"
    echo ""
    echo "Arguments:"
    echo "  barcode_number: Barcode number (01-22) or 'all' for all barcodes"
    echo "  region: Optional genomic region in format 'chr:start-end' (single barcode only)"
    echo "  bam_suffix: BAM file suffix (default: sorted.bam)"
    echo "  date: Date for single barcode analysis (20240326 or 20240328, default: 20240326)"
    echo ""
    echo "Default regions:"
    echo "  Human samples (barcodes 01-08): chr12:32679301-32745650"
    echo "  Mouse samples (barcodes 09-22): NC_000082.7:16130092-16176863"
    echo ""
    echo "Barcode mapping:"
    echo "  Human samples: barcodes 01-08 (uses fq_human.txt and Drp1_human_dge/)"
    echo "  Mouse samples: barcodes 09-22 (uses fq_mouse.txt and Drp1_mouse_dge/)"
    echo ""
    echo "Note: When using 'all', both 20240326 and 20240328 dates will be processed for each barcode."
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
    BAM_SUFFIX=${2:-"sorted.bam"}
    run_all_barcodes "${BAM_SUFFIX}"
else
    # Single barcode analysis
    BARCODE_NUM=$1
    BAM_SUFFIX=${3:-"sorted.bam"}
    DATE=${4:-"20240326"}
    
    # Determine species and set up paths based on barcode number
    SPECIES=""
    FQ_FILE=""
    BAM_PATH=""
    
    # Remove leading zeros for arithmetic operations
    BARCODE_NUM_CLEAN=$((10#${BARCODE_NUM}))
    
    # Determine if this is a human or mouse barcode and set appropriate default region
    if [ ${BARCODE_NUM_CLEAN} -ge 1 ] && [ ${BARCODE_NUM_CLEAN} -le 8 ]; then
        SPECIES="human"
        REGION=${2:-"chr12:32679301-32745650"}
        # Get FQ file path from fq_human.txt based on barcode and date
        if [ "${DATE}" = "20240326" ]; then
            FQ_LINE_NUM=$((BARCODE_NUM_CLEAN))
        else
            FQ_LINE_NUM=$((BARCODE_NUM_CLEAN + 8))
        fi
        FQ_FILE=$(sed -n "${FQ_LINE_NUM}p" fq_human.txt)
        BAM_PATH="./Drp1_human_dge/${DATE}_barcode$(printf "%02d" ${BARCODE_NUM_CLEAN})_${BAM_SUFFIX}"
    elif [ ${BARCODE_NUM_CLEAN} -ge 9 ] && [ ${BARCODE_NUM_CLEAN} -le 22 ]; then
        SPECIES="mouse"
        REGION=${2:-"NC_000082.7:16130092-16176863"}
        # Get FQ file path from fq_mouse.txt based on barcode and date
        if [ "${DATE}" = "20240326" ]; then
            FQ_LINE_NUM=$((BARCODE_NUM_CLEAN - 8))
        else
            FQ_LINE_NUM=$((BARCODE_NUM_CLEAN - 8 + 14))
        fi
        FQ_FILE=$(sed -n "${FQ_LINE_NUM}p" fq_mouse.txt)
        BAM_PATH="./Drp1_mouse_dge/${DATE}_barcode$(printf "%02d" ${BARCODE_NUM_CLEAN})_${BAM_SUFFIX}"
    else
        echo "ERROR: Invalid barcode number ${BARCODE_NUM}. Must be 1-22."
        exit 1
    fi
    
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