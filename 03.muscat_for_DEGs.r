# Differential expression analysis with muscat
# Load required packages
library(Seurat)
library(muscat)
library(SingleCellExperiment)
library(dplyr)
library(ggpubr)
library(scales)
# NaCl_WT vs MS_WT
# Select the sample groups for this comparison
seurat_list <- readRDS("WT_sos2_processed.rds")
subset_seu <- subset(seurat_list, subset = SampleGroup %in% c("MS_WT", "NaCl_WT"))
table(subset_seu$SampleGroup)  # Number of cells in each group
table(subset_seu$SampleGroup, subset_seu$celltype)  # Number of cells per group and cell type        
subset_seu$SampleGroup <- droplevels(subset_seu$SampleGroup)
table(subset_seu$SampleGroup)  # Number of cells in each group
table(subset_seu$SampleGroup, subset_seu$celltype)  # Number of cells per group and cell type
# Convert the Seurat object to a SingleCellExperiment object
sce <- as.SingleCellExperiment(subset_seu)
# Assign metadata required by muscat
sce$cluster_id <- sce$celltype  # Cell-type annotation
sce$sample_id <- sce$orig.ident  # Sample ID
sce$group_id <- sce$SampleGroup  # Experimental group
# Prepare the object for pseudobulk analysis
sce <- prepSCE(sce, 
               kid = "cluster_id",  # Cell-type column
               gid = "group_id",    # Group column
               sid = "sample_id",   # Sample column
               drop = TRUE)         # Drop unused factor levels
sce$group_id <- factor(sce$group_id, levels = c("MS_WT", "NaCl_WT"))
nk <- length(kids <- levels(sce$cluster_id))
ng <- length(sids <- levels(sce$sample_id))
names(kids) <- kids
names(sids) <- sids
t(table(sce$cluster_id, sce$sample_id))                                  
pb <- aggregateData(sce, 
                    assay = "counts", 
                    fun = "sum", 
                    by = c("cluster_id", "sample_id"))
colData(pb)
pb_mds <- pbMDS(pb)
# Run pseudobulk differential expression analysis with pbDS
pb$group_id
pb$group_id <- gsub("-", "_", pb$group_id)
pb$group_id <- factor(pb$group_id)
pb$group_id
res <- pbDS(pb, verbose = TRUE, filter = "none", method = "edgeR")
tbl <- res$table[[1]]
names(tbl)
# Combine the pbDS results into one data frame for visualization and DEG filtering
result_table <- resDS(sce, res, bind = "row", frq = FALSE, cpm = FALSE)
str(result_table)
# Select DEGs with |logFC| > 1 and p_adj.loc < 0.1
sig_genes <- subset(result_table, 
                    abs(logFC) > 1 &
                    p_adj.loc < 0.1)
# Sort genes by the absolute log fold change
sig_genes <- sig_genes[order(-abs(sig_genes$logFC)), ]
# Organize DEGs by cell type
sig_genes_by_cluster <- sig_genes %>% 
  group_by(cluster_id) %>% 
  arrange(desc(abs(logFC)))
deg_counts <- sig_genes %>% 
  group_by(cluster_id) %>% 
  summarise(deg_count = n()) %>% 
  tibble::deframe()
deg_counts
write.csv(result_table, file = "FC_P_muscat_WT_NaCl_MS_new.csv")
write.csv(sig_genes, file = "FC_1_padj_0.1_muscat_WT_NaCl_MS.csv")
# Separate upregulated and downregulated genes
up_genes <- sig_genes[sig_genes$logFC > 0, ]
down_genes <- sig_genes[sig_genes$logFC < 0, ]
head(down_genes)
head(up_genes)
write.csv(up_genes, file = "FC_1_padj_0.1_up_WT_NaCl_MS.csv")
write.csv(down_genes, file = "FC_1_padj_0.1_down_WT_NaCl_MS.csv")
# Count upregulated and downregulated genes by cell type
up_counts <- as.data.frame(table(up_genes$cluster_id))
down_counts <- as.data.frame(table(down_genes$cluster_id))
# Rename columns
colnames(up_counts) <- c("cluster_id", "up_count")
colnames(down_counts) <- c("cluster_id", "down_count")
# Merge the upregulated and downregulated gene counts
combined_counts <- merge(up_counts, down_counts, by = "cluster_id", all = TRUE)
# Replace missing counts with zero when a cell type has DEGs in only one direction
combined_counts[is.na(combined_counts)] <- 0
# Inspect the results
print(combined_counts)
write.csv(combined_counts, file = "WT_MS_NaCl_DEGs_number_FC_1_padj_0.1.csv")
# Map the number of DEGs per cell type onto the UMAP
a <- as.data.frame(sig_genes %>% group_by(cluster_id) %>% summarise(deg_count = n()))
a
# Extract metadata from the Seurat object
md <- seurat_list@meta.data
# Initialize the deg_count column with NA values
md$deg_count <- NA_real_
# Assign values only to the MS_WT and NaCl_WT groups
# Convert DEG counts to a named vector for matching by cell type
deg_vec <- setNames(a$deg_count, a$cluster_id)
table(md$celltype)  # Check the cell-type labels
# Assign DEG counts using mutate and case_when
md <- md %>%
  # Convert factor levels to character strings
  mutate(celltype_chr = as.character(celltype)) %>%
  # Match cell types to the named count vector
  mutate(
    deg_count = case_when(
      SampleGroup %in% c("MS_WT", "NaCl_WT") ~ deg_vec[celltype_chr],
      TRUE                              ~ NA_real_
    )
  ) %>%
  # Remove the temporary column
  select(-celltype_chr)
# Store the updated metadata in the Seurat object
seurat_list@meta.data <- md
table(is.na(seurat_list@meta.data$deg_count), seurat_list@meta.data$SampleGroup)
# Subset cells from the two comparison groups
sub <- subset(seurat_list, subset = SampleGroup %in% c("MS_WT","NaCl_WT"))
# Plot DEG counts on the UMAP with FeaturePlot
pdf("deg_counts_UMAP_WT_NaCl_MS_FC_1_padj_0.1.pdf")
p <- FeaturePlot(
  object     = sub,
  reduction  = "umap",
  features   = "deg_count",
  cols       = c("white", "#6fb3a8"),  # Set the low and high colors
  min.cutoff = 0,
  max.cutoff = 4500
) +
  ggtitle("MS_WT & NaCl_WT: deg_count on UMAP") +
  theme(plot.title = element_text(hjust = 0.5))
# Set a fixed color scale from white at 0 to #6fb3a8 at 4500
p + scale_colour_gradient(
      low    = "white",
      high   = "#6fb3a8",
      limits = c(0, 4500),
      oob    = squish      # Squish values outside the range to [0, 4500]
    )
print(p)
dev.off()

# NaCl_sos2 vs NaCl_WT
# Select the sample groups for this comparison
subset_seu <- subset(seurat_list, subset = SampleGroup %in% c("NaCl_WT", "NaCl_sos2"))
table(subset_seu$SampleGroup)  # Number of cells in each group
table(subset_seu$SampleGroup, subset_seu$celltype)  # Number of cells per group and cell type        
subset_seu$SampleGroup <- droplevels(subset_seu$SampleGroup)
table(subset_seu$SampleGroup)  # Number of cells in each group
table(subset_seu$SampleGroup, subset_seu$celltype)  # Number of cells per group and cell type
# Convert the Seurat object to a SingleCellExperiment object
sce <- as.SingleCellExperiment(subset_seu)
# Assign metadata required by muscat
sce$cluster_id <- sce$celltype  # Cell-type annotation
sce$sample_id <- sce$orig.ident  # Sample ID
sce$group_id <- sce$SampleGroup  # Experimental group
# Prepare the object for pseudobulk analysis
sce <- prepSCE(sce, 
               kid = "cluster_id",  # Cell-type column
               gid = "group_id",    # Group column
               sid = "sample_id",   # Sample column
               drop = TRUE)         # Drop unused factor levels
sce$group_id <- factor(sce$group_id, levels = c("NaCl_WT", "NaCl_sos2"))
nk <- length(kids <- levels(sce$cluster_id))
ng <- length(sids <- levels(sce$sample_id))
names(kids) <- kids
names(sids) <- sids
t(table(sce$cluster_id, sce$sample_id))                                  
pb <- aggregateData(sce, 
                    assay = "counts", 
                    fun = "sum", 
                    by = c("cluster_id", "sample_id"))
colData(pb)
pb_mds <- pbMDS(pb)
# Run pseudobulk differential expression analysis with pbDS
pb$group_id
pb$group_id <- gsub("-", "_", pb$group_id)
pb$group_id <- factor(pb$group_id)
pb$group_id
res <- pbDS(pb, verbose = TRUE, filter = "none", method = "edgeR")
tbl <- res$table[[1]]
names(tbl)
# Combine the pbDS results into one data frame for visualization and DEG filtering
result_table <- resDS(sce, res, bind = "row", frq = FALSE, cpm = FALSE)
str(result_table)
# Select DEGs with |logFC| > 1 and p_adj.loc < 0.1
sig_genes <- subset(result_table, 
                    abs(logFC) > 1 &
                    p_adj.loc < 0.1)
# Sort genes by the absolute log fold change
sig_genes <- sig_genes[order(-abs(sig_genes$logFC)), ]
# Organize DEGs by cell type
sig_genes_by_cluster <- sig_genes %>% 
  group_by(cluster_id) %>% 
  arrange(desc(abs(logFC)))
deg_counts <- sig_genes %>% 
  group_by(cluster_id) %>% 
  summarise(deg_count = n()) %>% 
  tibble::deframe()
deg_counts
write.csv(result_table, file = "FC_P_muscat_sos2_WT_NaCl_new.csv")
write.csv(sig_genes, file = "FC_1_padj_0.1_muscat_sos2_WT_NaCl.csv")
# Separate upregulated and downregulated genes
up_genes <- sig_genes[sig_genes$logFC > 0, ]
down_genes <- sig_genes[sig_genes$logFC < 0, ]
head(down_genes)
head(up_genes)
write.csv(up_genes, file = "FC_1_padj_0.1_up_sos2_WT_NaCl.csv")
write.csv(down_genes, file = "FC_1_padj_0.1_down_sos2_WT_NaCl.csv")
# Count upregulated and downregulated genes by cell type
up_counts <- as.data.frame(table(up_genes$cluster_id))
down_counts <- as.data.frame(table(down_genes$cluster_id))
# Rename columns
colnames(up_counts) <- c("cluster_id", "up_count")
colnames(down_counts) <- c("cluster_id", "down_count")
# Merge the upregulated and downregulated gene counts
combined_counts <- merge(up_counts, down_counts, by = "cluster_id", all = TRUE)
# Replace missing counts with zero when a cell type has DEGs in only one direction
combined_counts[is.na(combined_counts)] <- 0
# Inspect the results
print(combined_counts)
write.csv(combined_counts, file = "sos2_WT_NaCl_DEGs_number_FC_1_padj_0.1.csv")
# Map the number of DEGs per cell type onto the UMAP
a <- as.data.frame(sig_genes %>% group_by(cluster_id) %>% summarise(deg_count = n()))
a
# Extract metadata from the Seurat object
md <- seurat_list@meta.data
# Initialize the deg_count column with NA values
md$deg_count <- NA_real_
# Convert DEG counts to a named vector for matching by cell type
deg_vec <- setNames(a$deg_count, a$cluster_id)
table(md$celltype)  # Check the cell-type labels
# Assign DEG counts using mutate and case_when
md <- md %>%
  # Convert factor levels to character strings
  mutate(celltype_chr = as.character(celltype)) %>%
  # Match cell types to the named count vector
  mutate(
    deg_count = case_when(
      SampleGroup %in% c("NaCl_WT", "NaCl_sos2") ~ deg_vec[celltype_chr],
      TRUE                              ~ NA_real_
    )
  ) %>%
  # Remove the temporary column
  select(-celltype_chr)
# Store the updated metadata in the Seurat object
seurat_list@meta.data <- md
table(is.na(seurat_list@meta.data$deg_count), seurat_list@meta.data$SampleGroup)
# Subset cells from the two comparison groups
sub <- subset(seurat_list, subset = SampleGroup %in% c("NaCl_WT", "NaCl_sos2"))
# Plot DEG counts on the UMAP with FeaturePlot
pdf("deg_counts_UMAP_sos2_WT_NaCl_FC_1_padj_0.1.pdf")
p <- FeaturePlot(
  object     = sub,
  reduction  = "umap",
  features   = "deg_count",
  cols       = c("white", "#6fb3a8"),  # Set the low and high colors
  min.cutoff = 0,
  max.cutoff = 4500
) +
  ggtitle("NaCl_WT & NaCl_sos2: deg_count on UMAP") +
  theme(plot.title = element_text(hjust = 0.5))
# Set a fixed color scale from white at 0 to #6fb3a8 at 4500
p + scale_colour_gradient(
      low    = "white",
      high   = "#6fb3a8",
      limits = c(0, 4500),
      oob    = squish      # Squish values outside the range to [0, 4500]
    )
print(p)
dev.off()

# MS_sos2 vs MS_WT
# Select the sample groups for this comparison
subset_seu <- subset(seurat_list, subset = SampleGroup %in% c("MS_WT", "MS_sos2"))
table(subset_seu$SampleGroup)  # Number of cells in each group
table(subset_seu$SampleGroup, subset_seu$celltype)  # Number of cells per group and cell type        
subset_seu$SampleGroup <- droplevels(subset_seu$SampleGroup)
table(subset_seu$SampleGroup)  # Number of cells in each group
table(subset_seu$SampleGroup, subset_seu$celltype)  # Number of cells per group and cell type
# Convert the Seurat object to a SingleCellExperiment object
sce <- as.SingleCellExperiment(subset_seu)
# Assign metadata required by muscat
sce$cluster_id <- sce$celltype  # Cell-type annotation
sce$sample_id <- sce$orig.ident  # Sample ID
sce$group_id <- sce$SampleGroup  # Experimental group
# Prepare the object for pseudobulk analysis
sce <- prepSCE(sce, 
               kid = "cluster_id",  # Cell-type column
               gid = "group_id",    # Group column
               sid = "sample_id",   # Sample column
               drop = TRUE)         # Drop unused factor levels
sce$group_id <- factor(sce$group_id, levels = c("MS_WT", "MS_sos2"))
nk <- length(kids <- levels(sce$cluster_id))
ng <- length(sids <- levels(sce$sample_id))
names(kids) <- kids
names(sids) <- sids
t(table(sce$cluster_id, sce$sample_id))                                  
pb <- aggregateData(sce, 
                    assay = "counts", 
                    fun = "sum", 
                    by = c("cluster_id", "sample_id"))
colData(pb)
pb_mds <- pbMDS(pb)
# Run pseudobulk differential expression analysis with pbDS
pb$group_id
pb$group_id <- gsub("-", "_", pb$group_id)
pb$group_id <- factor(pb$group_id)
pb$group_id
res <- pbDS(pb, verbose = TRUE, filter = "none", method = "edgeR")
tbl <- res$table[[1]]
names(tbl)
# Combine the pbDS results into one data frame for visualization and DEG filtering
result_table <- resDS(sce, res, bind = "row", frq = FALSE, cpm = FALSE)
str(result_table)
# Select DEGs with |logFC| > 1 and p_adj.loc < 0.1
sig_genes <- subset(result_table, 
                    abs(logFC) > 1 &
                    p_adj.loc < 0.1)
# Sort genes by the absolute log fold change
sig_genes <- sig_genes[order(-abs(sig_genes$logFC)), ]
# Organize DEGs by cell type
sig_genes_by_cluster <- sig_genes %>% 
  group_by(cluster_id) %>% 
  arrange(desc(abs(logFC)))
deg_counts <- sig_genes %>% 
  group_by(cluster_id) %>% 
  summarise(deg_count = n()) %>% 
  tibble::deframe()
deg_counts
write.csv(result_table, file = "FC_P_muscat_sos2_WT_MS_new.csv")
write.csv(sig_genes, file = "FC_1_padj_0.1_muscat_sos2_WT_MS.csv")
# Separate upregulated and downregulated genes
up_genes <- sig_genes[sig_genes$logFC > 0, ]
down_genes <- sig_genes[sig_genes$logFC < 0, ]
head(down_genes)
head(up_genes)
write.csv(up_genes, file = "FC_1_padj_0.1_up_sos2_WT_MS.csv")
write.csv(down_genes, file = "FC_1_padj_0.1_down_sos2_WT_MS.csv")
# Count upregulated and downregulated genes by cell type
up_counts <- as.data.frame(table(up_genes$cluster_id))
down_counts <- as.data.frame(table(down_genes$cluster_id))
# Rename columns
colnames(up_counts) <- c("cluster_id", "up_count")
colnames(down_counts) <- c("cluster_id", "down_count")
# Merge the upregulated and downregulated gene counts
combined_counts <- merge(up_counts, down_counts, by = "cluster_id", all = TRUE)
# Replace missing counts with zero when a cell type has DEGs in only one direction
combined_counts[is.na(combined_counts)] <- 0
# Inspect the results
print(combined_counts)
write.csv(combined_counts, file = "sos2_WT_MS_DEGs_number_FC_1_padj_0.1.csv")
# Map the number of DEGs per cell type onto the UMAP
a <- as.data.frame(sig_genes %>% group_by(cluster_id) %>% summarise(deg_count = n()))
a
# Extract metadata from the Seurat object
md <- seurat_list@meta.data
# Initialize the deg_count column with NA values
md$deg_count <- NA_real_
# Convert DEG counts to a named vector for matching by cell type
deg_vec <- setNames(a$deg_count, a$cluster_id)
table(md$celltype)  # Check the cell-type labels
# Assign DEG counts using mutate and case_when
md <- md %>%
  # Convert factor levels to character strings
  mutate(celltype_chr = as.character(celltype)) %>%
  # Match cell types to the named count vector
  mutate(
    deg_count = case_when(
      SampleGroup %in% c("MS_WT", "MS_sos2") ~ deg_vec[celltype_chr],
      TRUE                              ~ NA_real_
    )
  ) %>%
  # Remove the temporary column
  select(-celltype_chr)
# Store the updated metadata in the Seurat object
seurat_list@meta.data <- md
table(is.na(seurat_list@meta.data$deg_count), seurat_list@meta.data$SampleGroup)
# Subset cells from the two comparison groups
sub <- subset(seurat_list, subset = SampleGroup %in% c("MS_WT", "MS_sos2"))
# Plot DEG counts on the UMAP with FeaturePlot
pdf("deg_counts_UMAP_sos2_WT_MS_FC_1_padj_0.1.pdf")
p <- FeaturePlot(
  object     = sub,
  reduction  = "umap",
  features   = "deg_count",
  cols       = c("white", "#6fb3a8"),  # Set the low and high colors
  min.cutoff = 0,
  max.cutoff = 4500
) +
  ggtitle("MS_WT & MS_sos2: deg_count on UMAP") +
  theme(plot.title = element_text(hjust = 0.5))
# Set a fixed color scale from white at 0 to #6fb3a8 at 4500
p + scale_colour_gradient(
      low    = "white",
      high   = "#6fb3a8",
      limits = c(0, 4500),
      oob    = squish      # Squish values outside the range to [0, 4500]
    )
print(p)
dev.off()