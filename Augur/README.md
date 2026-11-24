# Augur analysis

This directory contains the scripts used to compare transcriptional differences across Arabidopsis root cell types using Augur.

## Scripts

| Script                       | Comparison                                                  |
| ---------------------------- | ----------------------------------------------------------- |
| `01.augur_WT_NaCl_vs_MS.r`   | WT under NaCl treatment versus WT under control conditions. |
| `02.augur_sos2_vs_WT_MS.r`   | sos2 versus WT under control conditions.                    |
| `03.augur_sos2_vs_WT_NaCl.r` | sos2 versus WT under NaCl treatment.                        |

## Input

The scripts use an annotated Seurat object containing the following metadata columns:

* `celltype`: annotated cell type of each cell.

* `SampleGroup`: genotype and treatment group of each cell.

The Seurat object should also contain a UMAP reduction for visualization.&#x20;

## Subsampling and Augur analysis

For each comparison, cells are grouped by cell type and sample group. In each repeat, up to 150 cells are randomly selected from each cell-type and sample-group combination. If fewer than 150 cells are available, all available cells in that group are selected. The subsampling procedure is repeated 50 times with a different random seed for each repeat. For each repeat, Augur is run with the following parameters:

```r
subsample_size = 50
n_subsamples = 10
n_threads = 4
```

The mean area under the receiver operating characteristic curve (AUC) is extracted for each cell type in each repeat.

## Output

Each script generates:

* A CSV file containing the AUC results from all 50 repeats.

* A CSV file containing the mean AUC and standard deviation for each
  cell type.

* A PDF file displaying mean cell-type AUC values on the UMAP.

* An RDS file containing intermediate results from the repeated
  subsampling analysis.

## Requirements

The scripts require R and the packages loaded at the beginning of each script, including:

* Augur v1.0.3

* Seurat v5.5.0

* dplyr v1.2.1

* ggplot2 v4.0.3

* RhpcBLASctl v0.23.42

## Notes

* Update the input and output paths before running the scripts.

* A random seed is set for each repeat to support reproducibility.

