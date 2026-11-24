library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library(Augur)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(patchwork)
library(viridis)
library(Matrix)
library(sparseMatrixStats)
library(parsnip)
library(recipes)
library(yardstick)
library(scales)
# Update these paths to match your local data and output directories.
setwd("1.WT_NaClvsMS")
seurat_data <- readRDS("WT_sos2_processed.rds")
# Rename the original sample-group labels after loading the Seurat object.
seurat_data$SampleGroup <- as.character(seurat_data$SampleGroup)
seurat_data$SampleGroup <- gsub("WT", "WT", seurat_data$SampleGroup, fixed = TRUE)
seurat_data$SampleGroup <- gsub("sos2", "sos2", seurat_data$SampleGroup, fixed = TRUE)
seurat_data$SampleGroup <- factor(seurat_data$SampleGroup)
seurat_list <- subset(seurat_data, subset = SampleGroup %in% c("MS_WT", "NaCl_WT"))
seurat_list$SampleGroup <- droplevels(seurat_list$SampleGroup)
table(seurat_list$celltype, seurat_list$SampleGroup)
head(seurat_list@assays)
rm(seurat_data)
gc()
if(inherits(seurat_list[["RNA"]], "Assay5")) {
  seurat_list <- JoinLayers(seurat_list)
}
seurat_list$group_for_sampling <- paste(seurat_list$celltype, seurat_list$SampleGroup, sep = "_")
n_repeat <- 50
all_results <- list()
# Sample up to 150 cells per cell type and sample group in each repeat.
max_cells <- 150 
for(i in 1:n_repeat){  
  cat("Running subsampling repeat:", i, "out of", n_repeat, "\n")  
  # Set a reproducible random seed for each repeat.
  set.seed(i + 1000)  
  # Randomly sample cells from each cell type and sample group.
  cells_use <- unlist(
    lapply(
      split(colnames(seurat_list), seurat_list$group_for_sampling),
      function(x){
        sample(x, min(max_cells, length(x))) 
      }
    )
  )
  # Subset the Seurat object to the sampled cells.
  seurat_small <- subset(seurat_list, cells = cells_use)
  # Run Augur with 10 internal subsamples of 50 cells per condition.
  augur_res <- calculate_auc(
    seurat_small,
    cell_type_col = "celltype",
    label_col = "SampleGroup",
    subsample_size = 50,  
    n_subsamples = 10,     
    n_threads = 4
  )
  # Save the mean AUC for each cell type in this repeat.
  auc_df <- augur_res$AUC
  auc_df$repeat_id <- i
  all_results[[i]] <- auc_df
  saveRDS(all_results, "augur_temp_results.rds")
  # Remove temporary objects and free memory before the next repeat.
  rm(seurat_small, augur_res, cells_use, auc_df)
  gc(verbose = FALSE)
}
# Combine AUC results from all 50 repeats.
all_auc <- bind_rows(all_results)
print(all_auc, n = 50)
#-----------------------------
# Calculate the mean AUC and standard deviation for each cell type.
#-----------------------------
summary_auc <- all_auc %>%
  group_by(cell_type) %>%
  summarise(
    mean_AUC = mean(auc),
    sd_AUC   = sd(auc),
    n_repeat = n()
  ) %>%
  arrange(desc(mean_AUC))

# Inspect and export the AUC results.
summary_auc
write.csv(
  all_auc,
  file = "WT_NaClvsMS_Augur_all_repeat50_results.csv",
  row.names = FALSE
)
# Export the summary statistics for plotting.
write.csv(
  summary_auc,
  file = "WT_NaClvsMS_Augur_mean_AUC_summary50.csv",
  row.names = FALSE
)
augur_input <- list(
  AUC = summary_auc %>% 
    rename(cell_type = cell_type, auc = mean_AUC)
)
# Display mean cell-type AUC values on the UMAP.
pdf("WT_NaClvsMS_Augur_umap.pdf")
p <- plot_umap(
  augur = augur_input,            
  sc = seurat_list,               
  mode = "default",               
  reduction = "umap",             
  palette = c("grey90", "#ff8831"), 
  cell_type_col = "celltype"   
)
p <- p + labs(title = "WT NaCl vs MS") +
     theme(plot.title = element_text(hjust = 0.5))
print(p)
dev.off()
