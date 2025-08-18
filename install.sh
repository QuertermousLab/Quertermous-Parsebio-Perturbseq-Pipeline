#!/bin/bash

# -----------------------------
# Usage: ./install_spipe.sh [hg38|mm10]
# Default genome: hg38
# -----------------------------

GENOME=${1:-hg38}  # Default: hg38

# Set download URLs based on genome
if [ "$GENOME" == "hg38" ]; then
    FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.chr.gtf.gz"
    GENOME_NAME="GRCh38"
elif [ "$GENOME" == "mm10" ]; then
    FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
    GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm38.113.chr.gtf.gz"
    GENOME_NAME="mm10"
else
    echo "Unsupported genome: $GENOME"
    exit 1
fi

# -----------------------------
# 1. Create and activate conda environment
# -----------------------------
conda create -n spipe python=3.12.8 -c conda-forge -y
conda activate spipe

# -----------------------------
# 2. Install ParseBiosciences Pipeline
# -----------------------------
unzip -o ParseBiosciences-Pipeline.1.2.1.zip
cd ParseBiosciences-Pipeline.1.2.1
bash ./install_dependencies_conda.sh -i -y
pip install --no-cache-dir ./
cd ..

# -----------------------------
# 3. Install R and R packages
# -----------------------------
conda install -n spipe -c conda-forge r-base=4.3.1 r-reshape r-seurat r-devtools -y

Rscript -e "install.packages('devtools', repos='https://cloud.r-project.org')"
Rscript -e "devtools::install_github('katsevich-lab/sceptre')"

# -----------------------------
# 4. Download reference genome
# -----------------------------
wget -O ${GENOME_NAME}.fa.gz $FASTA_URL
wget -O ${GENOME_NAME}.gtf.gz $GTF_URL

# -----------------------------
# 5. Build reference index
# -----------------------------
split-pipe \
    --mode mkref \
    --genome_name $GENOME_NAME \
    --nthreads 16 \
    --fasta ./${GENOME_NAME}.fa.gz \
    --genes ./${GENOME_NAME}.gtf.gz \
    --output_dir ./${GENOME_NAME}_reference

echo "Installation finished for genome $GENOME!"
