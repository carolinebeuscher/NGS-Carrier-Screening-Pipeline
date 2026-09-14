#!/bin/bash
set -euo pipefail

##########################################
##### NGS carrier-screening pipeline #####
##########################################

# Worflow: raw reads -> QC/trimming -> alignment -> variant calling -> annotation & carrier detection for a known pathogenic variant
# Purpose: given an SRA read set, call variants against a single reference chromosome and report carrier status for a specified gene/variant of interest
# Usage: set the variables in the Config section below, then run: bash cmds.sh

### Config ###
# Edit these variables to specify the sample/run you want to process
SRA=SRR710118                     # SRA accession
SAMPLE=HG01890                    # 1000 Genomes sample ID (sequence_read source for $SRA)
REF_CHR=16                        # reference chromosome to align against (UCSC hg19 naming, e.g. 16, X, M)
GENE_CHR=11                       # chromosome containing the gene of interest (GRCh38 1000 Genomes naming)
GENE_NAME=HBB                     # gene of interest
GENE_REGION=5225464-5227071       # gene region (bp) on GENE_CHR, GRCh38 coordinates
VARIANT_ID=rs334                  # known pathogenic variant to screen for carrier status

### Genomic Data Exploration ###

PROJECT_DIR=$HOME/Documents/POP_GEN
mkdir -p $PROJECT_DIR

## SRA ##
SRA_DIR=$PROJECT_DIR/SRA
mkdir -p $SRA_DIR

# Pull from SRA and convert to fastq
prefetch $SRA -O $SRA_DIR

fastq-dump -O $SRA_DIR/$SRA/ \
    --split-files $SRA_DIR/$SRA/$SRA.sra

# Count the number of entries in the fastq file
cat $SRA_DIR/$SRA/${SRA}_1.fastq | grep "@$SRA" | wc -l

# Extract the sequences and write them to a file called $SRA.sequences
cat $SRA_DIR/$SRA/${SRA}_1.fastq | awk '{ if ((NR+2) % 4 == 0) print $0}' > $SRA_DIR/$SRA/$SRA.sequences

# Create a file with the length of each sequence
awk '{ print length }' $SRA_DIR/$SRA/$SRA.sequences > $SRA_DIR/$SRA/$SRA.sequences.length.txt

# Get the average read sequence length
awk '{ total += $1 } END { print total/NR }' $SRA_DIR/$SRA/$SRA.sequences.length.txt

## 1000 Genomes Project ##
GENOMES_DIR=$PROJECT_DIR/1000Genomes
mkdir -p $GENOMES_DIR

wget -P $GENOMES_DIR ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/$SAMPLE/sequence_read/${SRA}_1.filt.fastq.gz
wget -P $GENOMES_DIR ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/$SAMPLE/sequence_read/${SRA}_2.filt.fastq.gz



### Data Access and Quality Control ###

BEFORE_TRIM_DIR=$PROJECT_DIR/DataAccessQC/fastQCBeforeTrimming
mkdir -p $BEFORE_TRIM_DIR

TRIM_DIR=$PROJECT_DIR/DataAccessQC/Trimming
mkdir -p $TRIM_DIR

AFTER_TRIM_DIR=$PROJECT_DIR/DataAccessQC/fastQCAfterTrimming
mkdir -p $AFTER_TRIM_DIR

# Run fastqc for quality control on the downloaded fastq files
fastqc \
    -o $BEFORE_TRIM_DIR/ \
    $GENOMES_DIR/${SRA}_1.filt.fastq.gz

fastqc \
    -o $BEFORE_TRIM_DIR/ \
    $GENOMES_DIR/${SRA}_2.filt.fastq.gz

# Run fastp to trim the reads and remove low quality bases
fastp -i $GENOMES_DIR/${SRA}_1.filt.fastq.gz -o $TRIM_DIR/${SRA}_1.Trimmed.fastq.gz -f 12 -t 15
fastp -i $GENOMES_DIR/${SRA}_2.filt.fastq.gz -o $TRIM_DIR/${SRA}_2.Trimmed.fastq.gz -f 12 -t 15

# Run fastqc again on the trimmed reads to check quality
fastqc \
    -o $AFTER_TRIM_DIR/ \
    $TRIM_DIR/${SRA}_1.Trimmed.fastq.gz

fastqc \
    -o $AFTER_TRIM_DIR/ \
    $TRIM_DIR/${SRA}_2.Trimmed.fastq.gz



### Read Mapping and Variant Calling ###

MAPPING_DIR=$PROJECT_DIR/ReadMappingVariantCalling

RA_DIR=$MAPPING_DIR/ReadAlignment
mkdir -p $RA_DIR

# Download the reference genome (single chromosome, set via REF_CHR above)
wget -O $RA_DIR/chr${REF_CHR}.fa.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr${REF_CHR}.fa.gz
gunzip $RA_DIR/chr${REF_CHR}.fa.gz

# Index the fasta file (-a bwtsw for large/human genomes)
bwa index -a bwtsw $RA_DIR/chr${REF_CHR}.fa

# Align trimmed reads to the reference with bwa mem
# -t: number of threads; output is a SAM file
bwa mem \
    -t 24 \
    $RA_DIR/chr${REF_CHR}.fa \
    $TRIM_DIR/${SRA}_1.Trimmed.fastq.gz \
    $TRIM_DIR/${SRA}_2.Trimmed.fastq.gz > $RA_DIR/$SRA.sam

## Variant Calling
VARIANT_DIR=$MAPPING_DIR/VariantCalling
mkdir -p $VARIANT_DIR

# Convert SAM to BAM
samtools view \
    -bS $RA_DIR/$SRA.sam > $VARIANT_DIR/$SRA.bam

# Sort the BAM file
samtools sort \
    -o $VARIANT_DIR/$SRA.sorted.bam $VARIANT_DIR/$SRA.bam

# Count total sequences (including unmapped)
samtools view -c $VARIANT_DIR/$SRA.sorted.bam

# Count mapped reads only
samtools view -c -F 4 $VARIANT_DIR/$SRA.sorted.bam

samtools flagstat $VARIANT_DIR/$SRA.sorted.bam

samtools stats $VARIANT_DIR/$SRA.sorted.bam | grep '^SN'

# Make a pileup file
samtools mpileup \
    -f $RA_DIR/chr${REF_CHR}.fa $VARIANT_DIR/$SRA.sorted.bam > $VARIANT_DIR/$SRA.mpileup

# Call SNPs using VarScan
# --min-coverage 10:   minimum read depth at position
# --min-reads2 5:      minimum supporting reads for a variant call
# --min-avg-qual 15:   minimum average base quality
# --min-var-freq 0.01: minimum variant allele frequency (1%)
# --p-value 0.99:      p-value threshold for calling
varscan \
    mpileup2snp $VARIANT_DIR/$SRA.mpileup \
    --min-coverage 10 \
    --min-reads2 5 \
    --min-avg-qual 15 \
    --min-var-freq 0.01 \
    --p-value 0.99 > $VARIANT_DIR/${SRA}_snps.vcf

# Call indels using VarScan
varscan \
    mpileup2indel $VARIANT_DIR/$SRA.mpileup \
    --min-coverage 10 \
    --min-reads2 5 \
    --min-avg-qual 15 \
    --min-var-freq 0.01 \
    --p-value 0.99 > $VARIANT_DIR/$SRA.indels.vcf

# Count number of indels called (excluding header lines)
grep -v '^#' $VARIANT_DIR/$SRA.indels.vcf | wc -l

# Combine SNPs and indels into one VCF, headers first
grep '^#' $VARIANT_DIR/${SRA}_snps.vcf > $VARIANT_DIR/${SRA}_variants.vcf
grep -v '^#' $VARIANT_DIR/${SRA}_snps.vcf >> $VARIANT_DIR/${SRA}_variants.vcf
grep -v '^#' $VARIANT_DIR/$SRA.indels.vcf >> $VARIANT_DIR/${SRA}_variants.vcf

# Sort only the data lines (preserve headers at top)
(grep '^#' $VARIANT_DIR/${SRA}_variants.vcf; grep -v '^#' $VARIANT_DIR/${SRA}_variants.vcf | sort -k1,1 -k2,2n) > $VARIANT_DIR/${SRA}_vars.vcf

# Get per-position read counts
varscan \
    readcounts $VARIANT_DIR/$SRA.mpileup \
    --output-file $VARIANT_DIR/$SRA.readcounts

# Filter variants (require 25+ supporting reads)
varscan filter $VARIANT_DIR/${SRA}_snps.vcf \
    --min-reads2 25 \
    --output-file $VARIANT_DIR/$SRA.filtered.vcf

# Extract position + depth for REF_CHR and plot read depth across it
awk -v chr="$REF_CHR" '$1 == "chr"chr || $1 == chr' $VARIANT_DIR/$SRA.mpileup > $VARIANT_DIR/chr${REF_CHR}.pileup
awk '{print $2, $4}' $VARIANT_DIR/chr${REF_CHR}.pileup > $VARIANT_DIR/chr${REF_CHR}_depth.txt

python3 <<PYEOF
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("$VARIANT_DIR/chr${REF_CHR}_depth.txt", sep=' ', names=['Position', 'Depth'])

plt.figure(figsize=(15, 6))
plt.plot(data['Position'], data['Depth'], linewidth=0.5, alpha=0.7)
plt.xlabel('Genomic Position on Chromosome $REF_CHR')
plt.ylabel('Read Depth')
plt.title('Read Depth Distribution Across Chromosome $REF_CHR')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print(f"Mean depth: {data['Depth'].mean():.2f}")
print(f"Median depth: {data['Depth'].median():.2f}")
print(f"Max depth: {data['Depth'].max()}")
print(f"Min depth: {data['Depth'].min()}")
PYEOF



### Variant Annotation and Functional Analysis ###
# Only annotate variants that are clinically important

ANNOTATE_DIR=$PROJECT_DIR/AnnotateVariants
mkdir -p $ANNOTATE_DIR

# Retrieve the 1000 Genomes VCF file for the chromosome carrying the gene of interest (GRCh38)
wget -O $ANNOTATE_DIR/1000Genomes.Chr${GENE_CHR}.GRCh38.vcf.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr${GENE_CHR}_GRCh38.vcf.gz

tabix -p vcf $ANNOTATE_DIR/1000Genomes.Chr${GENE_CHR}.GRCh38.vcf.gz

# Extract the gene region of interest
bcftools view \
    -r chr${GENE_CHR}:${GENE_REGION} \
    -Ov \
    -o $ANNOTATE_DIR/1000Genomes.Chr${GENE_CHR}.${GENE_NAME}.vcf \
    $ANNOTATE_DIR/1000Genomes.Chr${GENE_CHR}.GRCh38.vcf.gz

# Annotate the gene region with snpEff
snpEff -Xmx16g ann GRCh38.105 $ANNOTATE_DIR/1000Genomes.Chr${GENE_CHR}.${GENE_NAME}.vcf > $ANNOTATE_DIR/1000Genomes.Chr${GENE_CHR}.${GENE_NAME}.Annotated.vcf

# Check for the pathogenic variant of interest
grep "$VARIANT_ID" $ANNOTATE_DIR/1000Genomes.Chr${GENE_CHR}.${GENE_NAME}.Annotated.vcf

# How many people in the study carry this variant
bcftools query \
    -f'[%CHROM:%POS %SAMPLE %GT\n]' \
    -i "ID==\"$VARIANT_ID\" && GT=\"alt\"" $ANNOTATE_DIR/1000Genomes.Chr${GENE_CHR}.${GENE_NAME}.Annotated.vcf \
    > $ANNOTATE_DIR/1000Genomes.${VARIANT_ID}.CarrierInfo.txt

wc -l $ANNOTATE_DIR/1000Genomes.${VARIANT_ID}.CarrierInfo.txt
