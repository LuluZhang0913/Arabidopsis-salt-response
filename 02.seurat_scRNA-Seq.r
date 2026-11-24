# Seurat analysis of WT and sos2 root scRNA-seq data.
# Run this script in the same environment used for the original analysis.

library(Seurat)
library(harmony)

### 1. Read count matrices and inspect quality metrics ----------------------
setwd("/scRNA_data")
data_dirs <- list.files("WT_sos2/")
print(data_dirs)

seurat_list <- list()
for (i in seq_along(data_dirs)) {
  seurat_data <- Read10X(
    data.dir = file.path("WT_sos2", data_dirs[i]),
    gene.column = 1
  )
  seurat_list[[i]] <- CreateSeuratObject(
    counts = seurat_data,
    project = data_dirs[i],
    min.features = 200,
    min.cells = 3
  )
}

# Inspect the metadata for each sample.
for (i in seq_along(seurat_list)) {
  cat("Seurat object", i, "metadata:\n")
  print(head(seurat_list[[i]]@meta.data))
}

# Calculate the fractions of mitochondrial and plastid transcripts.
for (i in seq_along(seurat_list)) {
  seurat_list[[i]][["mt_percent"]] <- PercentageFeatureSet(
    seurat_list[[i]], pattern = "^gene:ATMG"
  )
  seurat_list[[i]][["pt_percent"]] <- PercentageFeatureSet(
    seurat_list[[i]], pattern = "^gene:ATCG"
  )
}

# Plot QC metrics before filtering; each sample is shown on a separate page.
pdf(
  "violin_plot_feature.pdf",
  width = 10, height = 8
)
for (i in seq_along(seurat_list)) {
  violin_plot <- VlnPlot(
    seurat_list[[i]],
    features = c("nFeature_RNA", "nCount_RNA", "mt_percent", "pt_percent"),
    ncol = 3
  )
  print(violin_plot)
}
dev.off()

cat("Cells before QC:\n")
print(sapply(seurat_list, ncol))

### 2. Filter cells and merge samples ----------------------------------------
seurat_list <- lapply(seurat_list, function(obj) {
  subset(
    obj,
    subset = nFeature_RNA > 200 &
      nFeature_RNA < 12000 &
      mt_percent < 5 &
      pt_percent < 5
  )
})

cat("Cells after QC:\n")
print(sapply(seurat_list, ncol))

# add.cell.ids prefixes repeated cell barcodes with their sample names.
seurat_list <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = data_dirs
)

### 3. Normalize, reduce dimensions, and cluster ----------------------------
options(future.globals.maxSize = 1000 * 1024^3)
seurat_list <- NormalizeData(
  seurat_list, normalization.method = "LogNormalize", scale.factor = 10000
)
seurat_list <- FindVariableFeatures(
  seurat_list, nfeatures = 2000, selection.method = "vst"
)
seurat_list <- ScaleData(
  seurat_list, features = VariableFeatures(seurat_list)
)
seurat_list <- RunPCA(
  seurat_list, features = VariableFeatures(seurat_list)
)

# Harmony is used for the UMAP embedding. The neighbor graph below retains
# Seurat's default PCA reduction, as in the original analysis.
seurat_list <- RunHarmony(
  seurat_list,
  reduction = "pca",
  group.by.vars = "orig.ident",
  reduction.save = "harmony"
)
seurat_list <- RunUMAP(seurat_list, reduction = "harmony", dims = 1:40)
seurat_list <- FindNeighbors(seurat_list, dims = 1:40)
seurat_list <- FindClusters(
  seurat_list, resolution = c(seq(0, 0.5, 0.1))
)

# Optional: save the processed Seurat object for downstream analysis.
# saveRDS(seurat_list, "WT_sos2_processed.rds")
