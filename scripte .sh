#!/bin/bash
#upstream analysis 
#plant samples single end reads for arabidopsis 
#downloding the data 
wget ("path to the data")
#making a directory to work in 
mkdir intiative 
cd intiative 
#uncomprise the data 
for file in *.gz; do     if [ -f "$file" ]; then         echo "Extracting $file...";         gunzip "$file";     else         echo "No .gz files found.";     fi; done
#testing data quality 
mkdir qc before 
for i in *.fq; do fastqc $i ;done
#removing low quality reads or adapters 
 for i *fq; do trimmomatic SE -phred33 $i ${i%.gz}_trimmed.fq.gz ILLUMINACLIP:adapter file -SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36; done
for i in *trimed.fq.gz; do gunzip -f $i; done
#qc again to be sure 
for i in *.fq; do fastqc $i ;done
# indexing the reference genome   

STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /path/to/genome_index/genome1_index \
     --genomeFastaFiles /path/to/genome_fasta_files/genome1.fa \
     --sjdbGTFfile /path/to/annotations/genome1.gtf \
     --sjdbOverhang 99
     
#aligning the samples with the indexing output 
# Define the input and output directories

# Define the path to the STAR genome index
GENOME_DIR="/path/to/genomeDir"

# Define the input and output directories
INPUT_DIR="/path/to/fastq_files"
OUTPUT_DIR="/path/to/output_files"

# List of sample names or FASTQ file prefixes (adjust to match your naming convention)
SAMPLES=("sample1" "sample2" "sample3")

# Loop through each sample and run STAR
for SAMPLE in "${SAMPLES[@]}"; do
    # Define input FASTQ file for single-end reads
    READ="${INPUT_DIR}/${SAMPLE}.fq"

    # Define output prefix
    OUT_PREFIX="${OUTPUT_DIR}/${SAMPLE}_"

    # Run STAR alignment
    STAR --runThreadN 4 \
         --genomeDir ${GENOME_DIR} \
         --readFilesIn ${READ} \
         --outFileNamePrefix ${OUT_PREFIX} \
         --outSAMtype BAM SortedByCoordinate
done
#creat the count matrix 

featureCounts -a /home/m-sharawy/task/ath_annotation.gtf -o /home/m-sharawy/saved/count.txt /home/m-sharawy/yarb/aligned/sample1.bam /home/m-sharawy/yarb/aligned/sample2.bam /home/m-sharawy/yarb/aligned/sample3.bam /home/m-sharawy/yarb/aligned/sample4.bam
#take the output and make downstream analysis 
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2728-2 
#ref for the pipeline #






