## hdWGCNA
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(Seurat)
library(hdWGCNA)
library(WGCNA)
allowWGCNAThreads(nThreads = 20)

# Load the integrated Seurat object and prepare the RNA assay
seurat_obj <- readRDS("/zll_storage/5.scRNA-Seq/2.seurat/1.Col-0_sos2-t/Col-0_sos2-t_cluster.rds")
# Join RNA assay layers and set RNA as the active assay
seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
DefaultAssay(seurat_obj) <- "RNA"
# Set the output directory
setwd("/zll_storage/5.scRNA-Seq/2.seurat/1.Col-0_sos2-t/6.hdWGCNA/1.global")
# Confirm that the existing UMAP embedding has been retained
pdf("umap_check.pdf")
DimPlot(seurat_obj, group.by = "SampleGroup")
dev.off()
# Initialize a global hdWGCNA analysis
seurat_obj <- SetupForWGCNA(
   seurat_obj,
   gene_select = "fraction", # Select network genes by detection frequency
   fraction = 0.05,          # Retain genes detected in at least 5% of all cells
   wgcna_name = "Root_Global"
)
# Construct metacells within each cell-type and experimental-group combination
# Confirm that the object contains a Harmony reduction; otherwise, use an available reduction such as PCA
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("celltype", "SampleGroup"),
  reduction = 'harmony', 
  k = 25,
  max_shared = 10,
  ident.group = "celltype"
)
# Normalize metacell expression profiles
seurat_obj <- NormalizeMetacells(seurat_obj)

# Manually construct the normalized metacell expression matrix for network analysis
# =======================================================
cat("Manually constructing the WGCNA input matrix...\n")

# A. Retrieve the metacell object and join assay layers
m_obj <- GetMetacellObject(seurat_obj)
m_obj[["RNA"]] <- JoinLayers(m_obj[["RNA"]])

# B. Retrieve the selected network genes
target_genes <- seurat_obj@misc$Root_Global$wgcna_genes
cat("Number of genes selected at a detection fraction > 0.05:", length(target_genes), "\n")
# Number of genes selected at a detection fraction > 0.05: 13,800
# C. Extract the normalized expression layer with the Seurat v5 LayerData interface
# The RNA assay data layer contains normalized expression values
expr_sparse <- LayerData(m_obj, assay = "RNA", layer = "data")

# D. Retain only the selected network genes
valid_genes <- intersect(target_genes, rownames(expr_sparse))
cat("Number of selected genes matched in the expression matrix:", length(valid_genes), "\n")
expr_subset <- expr_sparse[valid_genes, ]
# Number of selected genes matched in the expression matrix: 13,800
# E. Transpose the matrix so that metacells are rows and genes are columns
# This matrix size is suitable for in-memory network construction
datExpr <- t(as.matrix(expr_subset))

# F. Replace non-finite values before network construction
datExpr[!is.finite(datExpr)] <- 0

# G. Store the expression matrix in the hdWGCNA experiment slot
seurat_obj@misc$Root_Global$datExpr <- datExpr

# H. Verify the matrix dimensions
cat("WGCNA input matrix constructed successfully. Dimensions:\n")
print(dim(seurat_obj@misc$Root_Global$datExpr))
[1] 24249 13800
# ============================================================================
# Select the soft-thresholding power
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = "signed" # A signed network emphasizes positively correlated gene pairs
)
   Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
1      1  0.81900 20.600          0.938 7450.00  7.51e+03 8130.0
2      2  0.67800  8.840          0.939 4120.00  4.16e+03 4940.0
3      3  0.38800  3.820          0.934 2320.00  2.34e+03 3160.0
4      4  0.10800  1.350          0.960 1340.00  1.33e+03 2120.0
5      5  0.00476 -0.218          0.923  790.00  7.73e+02 1490.0
6      6  0.23100 -1.330          0.854  475.00  4.54e+02 1100.0
7      7  0.49900 -1.830          0.827  291.00  2.71e+02  834.0
8      8  0.72400 -2.150          0.839  183.00  1.64e+02  655.0
9      9  0.84500 -2.270          0.871  117.00  1.00e+02  527.0
10    10  0.90600 -2.310          0.895   76.90  6.22e+01  435.0
11    12  0.94400 -2.170          0.928   35.50  2.47e+01  313.0
12    14  0.93400 -1.960          0.922   18.00  1.02e+01  240.0
13    16  0.92300 -1.750          0.916   10.10  4.30e+00  192.0
14    18  0.90000 -1.610          0.891    6.16  1.88e+00  160.0
15    20  0.88700 -1.480          0.882    4.08  8.46e-01  137.0
16    22  0.88600 -1.390          0.889    2.89  3.92e-01  120.0
17    24  0.89400 -1.310          0.905    2.17  1.84e-01  107.0
18    26  0.91100 -1.250          0.930    1.69  8.85e-02   97.2
19    28  0.90400 -1.200          0.931    1.37  4.35e-02   89.2
20    30  0.89600 -1.160          0.926    1.14  2.18e-02   82.4

# Visualize scale-free topology fit and mean connectivity across tested powers
pdf("1.PlotSoftPowers.pdf")
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol=2)
dev.off()
# Select power 10 based on the topology fit and connectivity profile
# Construct the signed co-expression network and identify modules
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 10,      
  setDatExpr = FALSE,     # Use the manually stored metacell expression matrix
  tom_name = "Root_Global", # Prefix used when saving the TOM matrix
  overwrite_tom = TRUE    # Overwrite an existing TOM with the same name
)
pdf("2.PlotDendrogram.pdf")
PlotDendrogram(seurat_obj, main='Root Global hdWGCNA Dendrogram')
dev.off()

## Export the topological overlap matrix (TOM)
TOM <- GetTOM(seurat_obj)
TOM <- as.data.frame(TOM)
write.csv(TOM,file = "1.TOM.csv")
## Inspect gene-to-module assignments
seurat_obj@misc$Root_Global$wgcna_modules %>% head
          gene_name    module     color
AT1G01470 AT1G01470      grey      grey
AT1G01620 AT1G01620 turquoise turquoise
AT1G01240 AT1G01240    purple    purple
AT1G01350 AT1G01350 royalblue royalblue
AT1G01770 AT1G01770      grey      grey
AT1G01010 AT1G01010     black     black
table(seurat_obj@misc$Root_Global$wgcna_modules$module)
         grey     turquoise        purple     royalblue         black 
         8633           837           177            98           227 
         blue           red           tan       darkred          pink 
          583           237           141            94           213 
    lightcyan    darkorange   lightyellow         brown      darkgrey 
          119            67           100           288            75 
  greenyellow        yellow         green        salmon          cyan 
          154           273           266           133           129 
      skyblue  midnightblue   saddlebrown         white       magenta 
           60           127            52            65           194 
    darkgreen    lightgreen        grey60 darkturquoise        orange 
           82           112           116            79            69 


# Compute module eigengenes (MEs)
# Each ME summarizes the expression pattern of a co-expression module
seurat_obj <- ModuleEigengenes(seurat_obj) # No group.by.vars is supplied, so these are unharmonized MEs

# Compute eigengene-based gene connectivity (kME)
table(seurat_obj$celltype)
# kME measures the correlation between each gene and its module eigengene and is used to identify hub genes
seurat_obj <- ModuleConnectivity(seurat_obj)

# Visualize module eigengenes across the UMAP
# Each plot displays the cell-level eigengene values for one module
pdf("3.ModuleFeaturePlot.pdf")
ModuleFeaturePlot(
  seurat_obj,
  features = "hMEs", # Requires harmonized MEs; use "MEs" if ModuleEigengenes was run without group.by.vars
  order = TRUE
)
dev.off()
# Use a dot plot to compare module eigengene values across cell types
# This view summarizes the relative activity of each module across cell types

# To compare treatment groups within cell types, a combined grouping variable can be added if needed
# Here, cell types are displayed on the x-axis and modules on the y-axis
# Use split.by in a separate plot to compare control and salt-treated samples
# 1. Retrieve module eigengenes
hMEs <- GetMEs(seurat_obj, harmonized=TRUE)
# 2. Add eigengene values to the Seurat metadata for plotting
# This makes the module values available to Seurat plotting functions
seurat_obj <- AddMetaData(seurat_obj, hMEs)
# 3. Prepare module names for plotting and exclude the unassigned grey module
mods <- colnames(hMEs)
mods <- mods[mods != "grey"]
# 4. Generate the dot plot
pdf("4.Dotplot.pdf", width = 12, height = 8)
DotPlot(seurat_obj, features = mods, group.by = 'celltype') + 
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='white', low='blue') +
  ggtitle("Module Eigengenes by Cell Type")
dev.off()

# Add the module eigengenes to metadata
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, hMEs)

# Plot the module of interest by cell type and experimental group
pdf("5.Unknown_blue.pdf")
VlnPlot(
  seurat_obj, 
  features = "blue", 
  group.by = "celltype", 
  split.by = "SampleGroup", # Split the violin plot by the four experimental groups
  pt.size = 0 
)
dev.off()

# Extract the top 25 hub genes from each module
hub_genes <- GetHubGenes(
  seurat_obj,
  n_hubs = 25
)
ModuleNetworkPlot(
   seurat_obj,
   outdir = "modulenetworks",
   n_inner = 10,
   n_outer = 15,
   n_conns = Inf,
   plot_size = c(10,10),
   vertex.label.cex = 1
 )

# Inspect hub genes from the blue module
print(hub_genes[hub_genes$module == "blue", "gene_name"])

# Export the gene-to-module assignment table
gene_module_df <- GetModules(seurat_obj)
write.csv(gene_module_df, "Gene_Module_List.csv")

# Plot a combined network using representative hub genes from each module
options(future.globals.maxSize = 30 * 1024^3)  # Set the maximum future export size to 30 GB
# Plot the hub-gene network; the number of displayed genes can be adjusted
pdf("6.HubGeneNetworkPlot.pdf")
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, 
  n_other=2,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()
# Optionally restrict the visualization to a selected subset of modules
g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)
# Retrieve the list of modules
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
# Plot a hub-gene network for the first five non-grey modules
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 2, 
  n_other= 2,
  edge_prop = 0.75,
  mods = mods[1:5] # only select 5 modules
)
