# NGS Carrier-Screening Pipeline

`cmds.sh` takes a raw SRA read set through QC, alignment, variant calling, and
annotation to report carrier status for a specified gene/variant of interest. Additionally plots read depth.

## Workflow

```
raw reads -> QC/trimming -> alignment -> variant calling -> annotation & carrier detection
```

## Tools used

| Stage | Tool | Purpose |
|---|---|---|
| Data acquisition | [SRA Toolkit](https://github.com/ncbi/sra-tools) (`prefetch`, `fastq-dump`) | Download and convert an SRA read set to FASTQ |
| Data acquisition | `wget` | Pull FASTQ/reference/VCF files from 1000 Genomes and UCSC |
| Quality control | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Read quality reports, before and after trimming |
| Quality control | [fastp](https://github.com/OpenGene/fastp) | Read trimming (adapter/low-quality base removal) |
| Alignment | [BWA](https://github.com/lh3/bwa) (`bwa index`, `bwa mem`) | Index the reference chromosome and align trimmed reads |
| Alignment/variant calling | [SAMtools](https://www.htslib.org/) (`view`, `sort`, `flagstat`, `stats`, `mpileup`) | SAM/BAM conversion, sorting, alignment QC stats, pileup generation |
| Variant calling | [VarScan](https://varscan.sourceforge.net/) (`mpileup2snp`, `mpileup2indel`, `readcounts`, `filter`) | Call SNPs/indels from the pileup and filter by read support |
| Annotation | [tabix](https://www.htslib.org/doc/tabix.html) / [bcftools](https://samtools.github.io/bcftools/) (`view`, `query`) | Index VCFs, extract a gene region, query carrier genotypes |
| Annotation | [snpEff](https://pcingola.github.io/SnpEff/) | Functional annotation of variants in the gene region |
| Visualization | Python (`pandas`, `matplotlib`), invoked inline via `python3 <<PYEOF` | Read-depth plot across `REF_CHR` |

## Config

Edit these at the top of `cmds.sh` before running:

| Variable | Meaning |
|---|---|
| `SRA` | SRA accession to process |
| `SAMPLE` | 1000 Genomes sample ID (FASTQ source for `SRA`) |
| `REF_CHR` | Reference chromosome to align against (UCSC hg19 naming) |
| `GENE_CHR` | Chromosome carrying the gene of interest (GRCh38 1000 Genomes naming) |
| `GENE_NAME` | Gene of interest |
| `GENE_REGION` | Gene region (bp) on `GENE_CHR`, GRCh38 coordinates |
| `VARIANT_ID` | Known pathogenic variant to screen for carrier status |

## Usage

```
bash cmds.sh
```

Requires `python3` with `pandas` and `matplotlib` installed, for the inline read-depth plot.
