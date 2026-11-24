# Arabidopsis-salt-response
# Single-cell and chromatin accessibility analyses of Arabidopsis roots under salt stress

This repository contains the scripts used to analyze single-cell RNA sequencing (scRNA-seq) data from Arabidopsis roots and assay for transposase-accessible chromatin using sequencing (ATAC-seq) data from Arabidopsis seedlings. We used these analyses to investigate cell-type-specific transcriptional responses to salt stress and the contribution of SOS2 to root responses.

## Repository contents

| File or directory       | Description                                                                                        |
| ----------------------- | -------------------------------------------------------------------------------------------------- |
| `01.dnbc4tools.sh`      | Processing of DNBelab C Series scRNA-seq data using DNBC4tools.                                    |
| `02.seurat_scRNA-Seq.r` | Quality control, normalization, dimensionality reduction, and clustering using Seurat and Harmony. |
| `03.muscat_for_DEGs.r`  | Cell-type differential expression analysis using the pseudobulk approach implemented in muscat.    |
| `04.hdWGCNA.r`          | High-dimensional weighted gene co-expression network analysis using hdWGCNA.                       |
| `Augur/`                | Augur analyses comparing treatment or genotype groups across root cell types.                      |
| `ATAC_seq/`             | ATAC-seq read processing, peak calling, differential accessibility analysis, and peak annotation.  |

Additional information about the Augur and ATAC-seq scripts is provided in the README files within the corresponding directories.

## Software requirements

The analyses require R and the software packages used in the individual scripts. Major tools include:

* dnbc4tools v2.1.3

* Seurat v5.2.1

* Harmony v1.2.3

* muscat v1.16.0

* edgeR v4.0.16

* Augur v1.0.3

* hdWGCNA v0.4.8

* Trim Galore v0.6.10

* Bowtie2 v2.5.4

* SAMtools v1.22

* BEDTools v2.31.1

* MACS2 v2.2.9.1

* deepTools v3.5.6

* DiffBind v3.18.0

* ChIPseeker v1.42.1

## Input data

The scRNA-seq and ATAC-seq scripts require the corresponding sequencing data, reference genome, annotations, and intermediate files described in each script. The paths in the scripts reflect the original analysis environment. Please update the input and output paths before running the scripts.

## Data availability

The single-cell RNA-seq (scRNA-seq) and ATAC-seq datasets generated in this study have been deposited in the NCBI Gene Expression Omnibus (GEO) under accession numbers GSE354964 for scRNA-seq and GSE354885 for ATAC-seq.

