# ATAC-seq analysis

This directory contains the scripts used to process and analyze Arabidopsis ATAC-seq data. The scripts are organized according to the main analysis steps. They are annotated records of the analysis and require adjustment and verification before use as an end-to-end workflow.

## Scripts

| Script                             | Description                                                                             |
| ---------------------------------- | --------------------------------------------------------------------------------------- |
| `01.trim_and_align.sh`             | Adapter and quality trimming, read alignment, and BAM file generation.                  |
| `02.mark_duplicates_filter.sh`     | Duplicate marking and alignment filtering.                                              |
| `03.peak_calling_and_frip.sh`      | Peak calling with MACS2 and calculation of the fraction of reads or intervals in peaks. |
| `04.bigwig.sh`                     | Generation of RPKM-normalized bigWig files.                                             |
| `05.diffbind.R`                    | Consensus peak construction and differential accessibility analysis using DiffBind.     |
| `06.annotate_peaks.R`              | Annotation of ATAC-seq peaks using ChIPseeker.                                          |
| `07.annotate_differential_peaks.R` | Annotation of differentially accessible peaks.                                          |

## Input data

The workflow requires:

* Paired-end ATAC-seq FASTQ files.

* The Arabidopsis reference genome.

* The corresponding genome annotation file.

* A sample sheet for DiffBind.

* Additional peak files specified in the scripts.

## Software requirements

The scripts use the following software:

* FastQC v0.12.1

* Trim Galore v0.6.10

* Bowtie2 v2.5.4

* SAMtools v1.22

* BEDTools v2.31.1

* MACS2 v2.2.9.1

* deepTools v3.5.6

* DiffBind v3.18.0

* ChIPseeker v1.42.1

## Running the scripts

The scripts were originally run from different working directories. Before running them:

1. Update the input files and output directories.

2. Check the sample names and paired-end FASTQ file matching.

3. Prepare the reference genome, annotation, and DiffBind sample sheet.

4. Run the scripts in the order shown above, after checking the commands
   and intermediate files.

## Important checks

The scripts preserve the recorded analysis commands. The following details should be checked against the original analysis environment before reuse:

* The recorded SAMtools command marks duplicate reads. Confirm whether duplicates were removed in any additional step.

* Check that the format of the file supplied to MACS2 matches the specified input format.

* Confirm the chromosome names used to remove mitochondrial and plastid
  alignments.

* Check the sample sheet, replicate assignments, and contrast directions used for differential accessibility analysis.

* Confirm the input format and counting unit used for the FRiP calculation.

