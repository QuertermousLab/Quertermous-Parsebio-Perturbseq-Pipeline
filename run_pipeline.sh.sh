#!/bin/bash

# -----------------------------
# Usage:
# ./run_splitpipe.sh <GENOME_DIR> <CHEMISTRY> <SAMPLE_LIST> <CRISPR_GUIDES> <FASTQ_DIR>
# -----------------------------

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <GENOME_DIR> <CHEMISTRY> <SAMPLE_LIST> <CRISPR_GUIDES> <FASTQ_DIR>"
    exit 1
fi

GENOME_DIR=$1
CHEMISTRY=v2
SAMPLE_LIST=$3
CRISPR_GUIDES=$4
FASTQ_DIR=$5

# -----------------------------
# Automatically detect sample numbers from FASTQ files containing 'WT'
# -----------------------------
SAMPLES=()
for f in ${FASTQ_DIR}/*.fastq.gz; do
    # Extract number after 'WT' in the filename (anywhere in name)
    if [[ $(basename $f) =~ WT([0-9]+) ]]; then
        NUM=${BASH_REMATCH[1]}
        # Avoid duplicates
        if [[ ! " ${SAMPLES[@]} " =~ " ${NUM} " ]]; then
            SAMPLES+=($NUM)
        fi
    fi
done

# Sort sample numbers numerically
IFS=$'\n' SAMPLES=($(sort -n <<<"${SAMPLES[*]}"))
unset IFS

echo "Detected samples: ${SAMPLES[@]}"

# -----------------------------
# Loop over WT and EC samples
# -----------------------------
for i in "${SAMPLES[@]}"; do
    WT_OUT=WT${i}
    EC_OUT=EC${i}
    FQ1=$(ls ${FASTQ_DIR}/*WT${i}*_R1.fastq.gz | head -n 1)
    FQ2=$(ls ${FASTQ_DIR}/*WT${i}*_R2.fastq.gz | head -n 1)
    EC_FQ1=$(ls ${FASTQ_DIR}/*EC${i}*_R1.fastq.gz | head -n 1)
    EC_FQ2=$(ls ${FASTQ_DIR}/*EC${i}*_R2.fastq.gz | head -n 1)

    # Check if files exist
    if [[ ! -f "$FQ1" || ! -f "$FQ2" ]]; then
        echo "Warning: WT FASTQ files for sample $i not found, skipping..."
        continue
    fi
    if [[ ! -f "$EC_FQ1" || ! -f "$EC_FQ2" ]]; then
        echo "Warning: EC FASTQ files for sample $i not found, skipping EC step..."
    fi

    # Run WT split-pipe
    split-pipe --mode all \
        --chemistry $CHEMISTRY \
        --genome_dir $GENOME_DIR \
        --fq1 $FQ1 \
        --fq2 $FQ2 \
        --output_dir $WT_OUT \
        --samp_list $SAMPLE_LIST

    # Run EC split-pipe if files exist
    if [[ -f "$EC_FQ1" && -f "$EC_FQ2" ]]; then
        split-pipe --mode all \
            --chemistry $CHEMISTRY \
            --crispr \
            --crsp_guides $CRISPR_GUIDES \
            --parent_dir ./$WT_OUT \
            --output_dir $EC_OUT \
            --fq1 $EC_FQ1 \
            --fq2 $EC_FQ2
    fi
done

# -----------------------------
# Combine WT samples
# -----------------------------
WT_DIRS=()
for i in "${SAMPLES[@]}"; do
    if [[ -d ./WT${i} ]]; then
        WT_DIRS+=("./WT${i}")
    fi
done
split-pipe --mode comb --sublibraries "${WT_DIRS[@]}" --output_dir TPS_combined

# -----------------------------
# Combine EC samples
# -----------------------------
EC_DIRS=()
for i in "${SAMPLES[@]}"; do
    if [[ -d ./EC${i} ]]; then
        EC_DIRS+=("./EC${i}")
    fi
done
split-pipe --mode comb --parent_dir ./TPS_combined --output_dir ./gRNA_combined --sublibraries "${EC_DIRS[@]}"
