library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library(Augur)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(viridis)
library(Matrix)
library(sparseMatrixStats)
library(parsnip)
library(recipes)
library(yardstick)
library(scales)
setwd("3.sos2vsWT_NaCl")
seurat_data <- readRDS("WT_sos2_processed.rds")
seurat_list <- subset(seurat_data, subset = SampleGroup %in% c("NaCl_WT", "NaCl_sos2"))
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
# 每次最多抽 150 个（如果不足 150 就全抽）
max_cells <- 150 
for(i in 1:n_repeat){  
  cat("Running subsampling repeat:", i, "out of", n_repeat, "\n")  
  # 设置随机种子，确保结果可复现
  set.seed(i + 1000)  
  # 从每个组中随机抽取细胞
  cells_use <- unlist(
    lapply(
      split(colnames(seurat_list), seurat_list$group_for_sampling),
      function(x){
        sample(x, min(max_cells, length(x))) 
      }
    )
  )
  # 提取子集
  seurat_small <- subset(seurat_list, cells = cells_use)
  # 运行 Augur
  augur_res <- calculate_auc(
    seurat_small,
    cell_type_col = "celltype",
    label_col = "SampleGroup",
    subsample_size = 50,  
    n_subsamples = 10,     
    n_threads = 4
  )
  # 5. 提取并记录结果
  auc_df <- augur_res$AUC
  auc_df$repeat_id <- i
  all_results[[i]] <- auc_df
  saveRDS(all_results, "augur_temp_results.rds")
  # 6. 强制清理内存（极其重要！防止变量串线导致结果一模一样）
  rm(seurat_small, augur_res, cells_use, auc_df)
  gc(verbose = FALSE)
}
# 合并 50 次的结果
all_auc <- bind_rows(all_results)
print(all_auc, n = 50)
#-----------------------------
# 计算mean AUC和SD
#-----------------------------
summary_auc <- all_auc %>%
  group_by(cell_type) %>%
  summarise(
    mean_AUC = mean(auc),
    sd_AUC   = sd(auc),
    n_repeat = n()
  ) %>%
  arrange(desc(mean_AUC))

# 查看最终结果
summary_auc
write.csv(
  all_auc,
  file = "NaCl_sos2vsWT_Augur_all_repeat50_results.csv",
  row.names = FALSE
)
##用来绘制小提琴图
write.csv(
  summary_auc,
  file = "NaCl_sos2vsWT_Augur_mean_AUC_summary50.csv",
  row.names = FALSE
)
augur_input <- list(
  AUC = summary_auc %>% 
    rename(cell_type = cell_type, auc = mean_AUC)
)
##绘制UMAP图
pdf("NaCl_sos2vsWT_Augur_umap.pdf")
p <- plot_umap(
  augur = augur_input,            
  sc = seurat_list,               
  mode = "default",               
  reduction = "umap",             
  palette = c("grey90", "#ff8831"), 
  cell_type_col = "celltype"   
)
p <- p + labs(title = "sos2 vs WT NaCl") +
     theme(plot.title = element_text(hjust = 0.5))
print(p)
dev.off()
