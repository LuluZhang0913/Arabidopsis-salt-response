#!/usr/bin/env bash
# Generate a gene-by-cell count matrix from DNBelab C Series scRNA-seq data.
# Run this script from a working directory containing reference/ and fastq/.
# Adjust the paths and sample name below to match the files used for each run.

set -euo pipefail

# Software required: gffread and dnbc4tools (available on PATH).
THREADS=10
SAMPLE_NAME="sample"

# Arabidopsis TAIR10 genome and release 54 annotation.
REF_FASTA="reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
REF_GFF3="reference/Arabidopsis_thaliana.TAIR10.54.gff3"
REF_GTF="reference/Arabidopsis_thaliana.TAIR10.54.gtf"
REF_INDEX="0.reference"

# Paired FASTQ files from the cDNA and Oligo libraries of the SAME sample.
CDNA_R1="fastq/cDNA_R1.fq.gz"
CDNA_R2="fastq/cDNA_R2.fq.gz"
OLIGO_R1="fastq/oligo_R1.fq.gz"
OLIGO_R2="fastq/oligo_R2.fq.gz"  # Verify that this is the mate of OLIGO_R1.

mkdir -p "$REF_INDEX" logs

# 1. Convert the GFF3 annotation to GTF for reference construction.
gffread "$REF_GFF3" -T -o "$REF_GTF"

# 2. Build the dnbc4tools reference in the same directory used by `rna run`.
dnbc4tools rna mkref \
    --fasta "$REF_FASTA" \
    --ingtf "$REF_GTF" \
    --species Arabidopsis_thaliana \
    --genomeDir "$REF_INDEX" \
    --threads "$THREADS" \
    > logs/mkref.stdout.log 2> logs/mkref.stderr.log

# 3. Process paired cDNA and Oligo reads to generate a gene-by-cell matrix.
dnbc4tools rna run \
    --name "$SAMPLE_NAME" \
    --cDNAfastq1 "$CDNA_R1" \
    --cDNAfastq2 "$CDNA_R2" \
    --oligofastq1 "$OLIGO_R1" \
    --oligofastq2 "$OLIGO_R2" \
    --genomeDir "$REF_INDEX" \
    --threads "$THREADS" \
    > "logs/${SAMPLE_NAME}.stdout.log" 2> "logs/${SAMPLE_NAME}.stderr.log"
