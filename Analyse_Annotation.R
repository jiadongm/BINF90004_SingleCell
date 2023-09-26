## Sources
# https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html 

## Set working directory
setwd("~/Teaching/BINF")

## Load necessary packages & functions
source("loadPkgs.R")
source("Viz.R")

## Prepare query single cell BMMC
if(F){
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
  # Simple QC
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
}
query <- bmmc

## Prepare reference data (stored in celldex)
reference <- MonacoImmuneData()
# Inspect
reference
reference@colData


## Common genes and rank transform
commonGenes <- intersect(rownames(reference), rownames(query))
length(commonGenes)
reference <- reference[commonGenes, ]
query <- query[commonGenes, ]


### SingleR -------------------------------------------------------------------
singleRassay <- "logcounts"
YtrainName <- "label.main"
sr_train <- trainSingleR(assay(reference, singleRassay),
                         labels = colData(reference)[, YtrainName])
sr_re <- classifySingleR(assay(query, singleRassay),
                         sr_train,
                         fine.tune = T)   
# Heatmap of annotation resutls
plotScoreHeatmap(sr_re)
# Store the annotations
query$SingleR_main <- sr_re$labels
YtrainName <- "label.fine"
sr_train <- trainSingleR(assay(reference, singleRassay),
                         labels = colData(reference)[, YtrainName])
sr_re <- classifySingleR(assay(query, singleRassay),
                              sr_train,
                              fine.tune = T)
# Heatmap of annotation resutls
plotScoreHeatmap(sr_re)
# Store the annotations
query$SingleR_fine <- sr_re$labels 

# Compare original labels and SingleR annotations
plotSankey(query$SingleR_main, query$celltype.l1)
plotSankey(query$SingleR_fine, query$celltype.l2)

# Diagnoistics
plotDeltaDistribution(sr_re)
plotScoreDistribution(sr_re)



















