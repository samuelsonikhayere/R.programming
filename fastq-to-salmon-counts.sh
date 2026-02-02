#! /bin/bash
set -euo pipefail
trap 'echo "Error on line $LINENO. Exiting."; exit 1' ERR

#--------------------------------------------------------------------
#Differential Gene Expressions (From Fastq to Salmon Quantification)
#Ikhayere, Samson Samuel (c) November, 2025
#--------------------------------------------------------------------

source ~/miniconda3/etc/profile.d/conda.sh

SECONDS=0

Project_Dir="/mnt/c/Users/USER/Documents/Data Analysis/my_bioinformatics_project" || exit
cd "$Project_Dir" || exit 

#Step 1: RUN FastQC
echo "FastQC Running..."
conda activate fastqc
mkdir -p data/processed 
fastqc data/raw/demo.fastq -o data/processed
conda deactivate
echo "Initial FastQC completed..."

#Step 2: Trim the Data (if neccessary)
echo "Trimmomatic Running..."
conda activate rna-seq 
mkdir -p data/Trimmed
java -jar "trimming/Trimmomatic-0.39/trimmomatic-0.39.jar" SE \
    -threads 4 -phred33 \
    data/raw/demo.fastq \
    data/Trimmed/trimmed_demo.fastq \
    TRAILING:10 \
    LEADING:3 
    
conda deactivate
echo "Data completely trimmed"

#Step 3: Run FastQC (Post-Trimming)
echo "FastQC on Trimmed Data Running..."
conda activate fastqc
fastqc data/Trimmed/trimmed_demo.fastq -o data/Trimmed
conda deactivate
echo "Post-Trim FastQC  completed..."

#Step 4: Build Salmon Index
echo "Setting up Salmon Index"
echo "Downloading the reference Transcriptome file"
conda activate salmon_env
mkdir -p data/salmon_index
cd "data/salmon_index" || exit


if [ ! -f "Homo_sapiens.GRCh38.cdna.all.fa" ]; then
    echo "Downloading the reference Transcriptome file..."
    echo ""
    if [ ! -f "Homo_sapiens.GRCh38.cdna.all.fa.gz" ]; then
        # Download only if compressed file is also missing
        wget "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" 
    fi
    # Attempt to unzip (this handles cases where .gz exists but .fa doesn't)
    gunzip -f "Homo_sapiens.GRCh38.cdna.all.fa.gz" 
fi
echo "reference Transcriptome file prepared \n"


echo "Downloading the genome file for decoy..."
if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
    echo "Downloading the Primary Assemble File..."
    echo ""

    if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" ]; then 
        wget "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    fi
    gunzip "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

fi
echo "genome file for decoy downloaded..."
echo ""

#Concatente the transcriptome and genome data

echo "Concatenating the transcriptome and genome files and extracting decoys..."
cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa > combined.fa

grep "^>" Homo_sapiens.GRCh38.dna.primary_assembly.fa | cut \
    -d " " -f 1 |sed 's/>//g' > decoys.txt

echo "Concatenation and decoy extraction completed..."
echo ""

salmon index -t combined.fa \
    -d decoys.txt \
    -i grch38_salmon_index_decoys \
    --threads 4 \
    -k 31

cd "$Project_Dir"
echo "Salmon Index Built..."
echo ""

#Step 5: Salmon Quantification

echo "Salmon Quant Running..."
conda activate salmon_env
mkdir -p data/quant
salmon quant -i data/salmon_index/grch38_salmon_index_decoys \
    -l A \
    -r data/Trimmed/trimmed_demo.fastq \
    -o data/quant/demo_quant \
    --threads 4 
echo "Quantification completed, file saved in 'data/quant'"
conda deactivate

Duration=$SECONDS
echo "$(($Duration/60)) minutes and $(($Duration%60)) seconds elapsed"
