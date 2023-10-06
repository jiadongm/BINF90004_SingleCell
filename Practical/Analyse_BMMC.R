## Set working directory
setwd("~/Documents/GitHub/BINF90004_SingleCell/Practical")

## Load necessary packages & functions
source("loadPkgs.R")
source("Viz.R")

## Install and load data
# Install data (only need to install once)
if(F){
  AvailableData()
  InstallData("bmcite")
} 
# Update old Seurat obj after updating Seurat pkg
bmmc <- UpdateSeuratObject(object = bmcite) 

# Subset 10% of cells to reduce computation
table(bmmc$celltype.l2)
celltypes <- bmmc$celltype.l2 %>% unique()
toKeepIdx <- c()
set.seed(9048)
for(i in 1:length(celltypes)){
  ctype <- celltypes[i]
  idx <- which(bmmc$celltype.l2 == ctype)
  keepHowMany <- max(round(length(idx)*0.1), 1)
  idx <- sample(idx, keepHowMany)
  toKeepIdx <- c(toKeepIdx, idx)
}
bmmc <- bmmc[,toKeepIdx]

# Simple QC (only delete genes with zero expression in all cells)
feat2keep <- rownames(bmmc)[rowMeans(bmmc[["RNA"]]$counts == 0) < 1]
bmmc <- bmmc[feat2keep,]

# Convert Seurat obj to sce
bmmc <- SingleCellExperiment(
  list(counts = bmmc[["RNA"]]$counts),
  colData = bmmc@meta.data
)

# Scran norm
set.seed(5202056)
bmmc <- computeSumFactors(bmmc, cluster=quickCluster(bmmc))
bmmc <- logNormCounts(bmmc)
if(F){
  temp <- assay(bmmc, "counts")
  temp <- log(t(t(temp)/colSums(temp)) * 1000000 + 1)
  assay(bmmc, "cpm") <- temp
}


### Dimension reduction --------------------------------------------------
# PCA using svds
dat <- assay(bmmc, "logcounts") %>%
  as.matrix() %>%
  t() %>%
  scale(center = T, scale = F)
svd_res <- svds(dat, k = 50)
pc_score <-  svd_res$u %*% diag(svd_res$d)
pc_load <- svd_res$v
dimnames(pc_score) <- list(rownames(dat), paste0("comp", 1:ncol(pc_score)))
dimnames(pc_load) <- list(colnames(dat), paste0("comp", 1:ncol(pc_load)))

# Scree plot (again, not very useful)
totVar <- sum(dat^2)
(svd_res$d^2/totVar) %>% plot()

# UMAP + PCA visualisation 
umap_res <- umap(pc_score %>% as.data.frame() %>% select(comp1:comp50))
umap_score <- umap_res$layout
colnames(umap_score) <- paste0("comp", 1:ncol(umap_score))
p1 <-
  pc_score %>%
  as.data.frame() %>%
  mutate(labs = bmmc$celltype.l2) %>%
  ggplot(aes(x = comp1, y = comp2)) +
  geom_point(aes(colour = labs), alpha = 1, stroke = NA, size = 3)
p2 <-
umap_score %>%
  as.data.frame() %>%
  mutate(labs = bmmc$celltype.l2) %>%
  ggplot(aes(x = comp1, y = comp2)) +
  geom_point(aes(colour = labs), alpha = 1, stroke = NA, size = 3)
ggarrange(p1, p2, common.legend = T)


# matrix plot for PCs, allowing us to see more PCs than the first 2
# (Probably 20 is a good number of PCs?)
plot_dat <- pc_score %>%
  as.data.frame() %>%
  mutate(labs = bmmc$celltype.l2)
outPlots <- matrixPlot(plot_dat, comp_idx = 1:5, groupKey = "labs")
outPlots <- matrixPlot(plot_dat, comp_idx = 6:10, groupKey = "labs")
outPlots <- matrixPlot(plot_dat, comp_idx = 11:15, groupKey = "labs")
outPlots <- matrixPlot(plot_dat, comp_idx = 16:20, groupKey = "labs")
outPlots <- matrixPlot(plot_dat, comp_idx = 21:25, groupKey = "labs")
outPlots <- matrixPlot(plot_dat, comp_idx = 26:30, groupKey = "labs")
outPlots <- matrixPlot(plot_dat, comp_idx = 41:45, groupKey = "labs")
outPlots <- matrixPlot(plot_dat, comp_idx = 45:50, groupKey = "labs")



# Look at specific cell types
plot_dat %>%
  group_by(labs) %>%
  summarise(mean_comp = mean(comp6)) %>%
  print(n=27)
selectedTypes <- c("pDC", "Prog_RBC")
lvls <- c("other", selectedTypes)
manualCols <- c("gray", "red", "blue")
labs <- bmmc$celltype.l2 %>% as.character()
labs[!(labs %in% selectedTypes)] <- "other"
labs <- factor(labs, levels = lvls)
plot_dat <- pc_score %>%
  as.data.frame() %>%
  mutate(labs = labs) %>%
  arrange(labs)
outPlots <- matrixPlot(plot_dat, comp_idx = 1:5, groupKey = "labs", manualCol = manualCols)
outPlots <- matrixPlot(plot_dat, comp_idx = 6:10, groupKey = "labs", manualCol = manualCols)
outPlots <- matrixPlot(plot_dat, comp_idx = 11:15, groupKey = "labs", manualCol = manualCols)
outPlots <- matrixPlot(plot_dat, comp_idx = 16:20, groupKey = "labs", manualCol = manualCols)
outPlots <- matrixPlot(plot_dat, comp_idx = 21:25, groupKey = "labs", manualCol = manualCols)
outPlots <- matrixPlot(plot_dat, comp_idx = 26:30, groupKey = "labs", manualCol = manualCols)
outPlots <- matrixPlot(plot_dat, comp_idx = 41:45, groupKey = "labs", manualCol = manualCols)
outPlots <- matrixPlot(plot_dat, comp_idx = 45:50, groupKey = "labs", manualCol = manualCols)












