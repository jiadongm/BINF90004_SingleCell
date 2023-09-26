## Sources
# https://satijalab.org/seurat/articles/pbmc3k_tutorial
# https://bioconductor.org/books/3.13/OSCA.basic/normalization.html 
# https://satijalab.org/seurat/articles/sctransform_v2_vignette 

## Set working directory
setwd("~/Teaching/BINF")

## Load necessary packages & functions
source("loadPkgs.R")
source("Viz.R")

## Install and load data
# Install data (only need to install once)
if(F){
  AvailableData()
  InstallData("pbmc3k")
} 
# Update old Seurat obj after updating Seurat pkg
pbmc3k <- UpdateSeuratObject(object = pbmc3k) 

## Structure of Seurat object
pbmc3k
pbmc3k@active.assay
pbmc3k@assays$RNA
pbmc3k@assays$RNA$counts
pbmc3k@meta.data

## QC
pbmc3k[["percent.mt"]] <- PercentageFeatureSet(pbmc3k, pattern = "^MT-")
feat2keep <- rownames(pbmc3k)[rowMeans(pbmc3k[["RNA"]]$counts == 0) < 1]
pbmc3k <- subset(pbmc3k, 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5,
                 features = feat2keep)

## Distribution of library size
libSizes <- pbmc3k[["RNA"]]$counts %>% colSums()
plot_dat <- data.frame(libSizes = libSizes)
plot_dat %>% 
  ggplot(aes(x = libSizes)) +
  geom_density(bw = "SJ") +
  geom_point(aes(y = 0), alpha = 0.2, shape = 3)

### Normalisation ----------------------------------------------------------------
## Shifted log
pbmc3k <- NormalizeData(pbmc3k, normalization.method = "LogNormalize", scale.factor = 1)
# Compare to step-by-step log-norm
if(F){
  temp <- pbmc3k[["RNA"]]$counts %>% as.matrix()
  temp <- t(t(temp)/colSums(temp))
  temp <- log(temp + 1) %>% as.sparse()
  temp[1:10,1:20]
  pbmc3k[["RNA"]]$data[1:10, 1:20]
}

# Marginal plots - raw counts
plot_dat <- pbmc3k[["RNA"]]$counts %>% as.matrix() %>% t()
marginMeans <- marginMeanPlot(plot_dat)
geneIdx <- names(marginMeans)[seq(1, 5000, length.out = 16)]
outPlots <- marginPlot(plot_dat %>% data.frame(), featIdx = geneIdx, 
                       plotNrows = 4, plotNcols = 4,
                       omitZeros = F)
outPlots <- marginPlot(plot_dat %>% data.frame(), featIdx = geneIdx, 
                       plotNrows = 4, plotNcols = 4,
                       omitZeros = T)

# Marginal plots - relative counts
plot_dat <- pbmc3k[["RNA"]]$counts %>% as.matrix() %>% t()
plot_dat <- plot_dat/rowSums(plot_dat)
outPlots <- marginPlot(plot_dat %>% data.frame(), featIdx = geneIdx, 
                       plotNrows = 4, plotNcols = 4,
                       omitZeros = T)

# Marginal plots - log-normalised
plot_dat <- pbmc3k[["RNA"]]$data %>% as.matrix() %>% t()
outPlots <- marginPlot(plot_dat %>% data.frame(), featIdx = geneIdx, 
                       plotNrows = 4, plotNcols = 4,
                       omitZeros = T)

# Margin plots - log-normalised but with different shift
pbmc3k <- NormalizeData(pbmc3k, normalization.method = "LogNormalize", scale.factor = 1000000)
plot_dat <- pbmc3k[["RNA"]]$data %>% as.matrix() %>% t()
outPlots <- marginPlot(plot_dat %>% data.frame(), featIdx = geneIdx, 
                       plotNrows = 4, plotNcols = 4,
                       omitZeros = T)
mean(libSizes)
pbmc3k <- NormalizeData(pbmc3k, normalization.method = "LogNormalize", scale.factor = mean(libSizes))
plot_dat <- pbmc3k[["RNA"]]$data %>% as.matrix() %>% t()
outPlots <- marginPlot(plot_dat %>% data.frame(), featIdx = geneIdx, 
                       plotNrows = 4, plotNcols = 4,
                       omitZeros = T)

## Scran
# Convert Seurat obj to sce
pbmc_sce <- as.SingleCellExperiment(pbmc3k)

# Look at sce object
pbmc_sce
pbmc_sce@assays
pbmc_sce@assays$counts
pbmc_sce@assays$logcounts
assay(pbmc_sce, "counts")
pbmc_sce@colData

# Scran norm
set.seed(5202056)
clust_scran <- quickCluster(pbmc_sce)
table(clust_scran)
pbmc_sce <- computeSumFactors(pbmc_sce, cluster=clust_scran)
summary(pbmc_sce$sizeFactor)
pbmc_sce <- logNormCounts(pbmc_sce)
plot_dat <- assay(pbmc_sce, "logcounts") %>% as.matrix() %>% t()
outPlots <- marginPlot(plot_dat %>% data.frame(), featIdx = geneIdx, 
                       plotNrows = 4, plotNcols = 4,
                       omitZeros = T)



### Dimension reduction --------------------------------------------------
# PCA using svds
dat <- assay(pbmc_sce, "logcounts") %>%
  as.matrix() %>%
  t() %>%
  scale(center = T, scale = T)
svd_res <- svds(dat, k = 50)
pc_score <-  svd_res$u %*% diag(svd_res$d)
pc_load <- svd_res$v
dimnames(pc_score) <- list(rownames(dat), paste0("comp", 1:ncol(pc_score)))
dimnames(pc_load) <- list(colnames(dat), paste0("comp", 1:ncol(pc_load)))

pc_score %>%
  as.data.frame() %>%
  mutate(labs = pbmc_sce$seurat_annotations) %>%
  ggplot(aes(x = comp1, y = comp2)) +
  geom_point(aes(colour = labs))

# UMAP
umap_res <- umap(pc_score %>% as.data.frame() %>% select(comp1:comp10))
umap_score <- umap_res$layout
colnames(umap_score) <- paste0("comp", 1:ncol(umap_score))

p <- 
  umap_score %>%
  as.data.frame() %>%
  mutate(labs = pbmc_sce$seurat_annotations) %>%
  ggplot(aes(x = comp1, y = comp2)) +
  geom_point(aes(colour = labs))
p
p +
  scale_y_continuous(limits = c(-10, 11))








