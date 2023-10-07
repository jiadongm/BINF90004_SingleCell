### Alternative datasets for student report
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratData)

## scRNA-seq dataset1: pbmcsca -----------------------------------------------
## This dataset contains cells sequenced using 9 methods (see pbmcsca$Method)
## This allows us to see batch effects and how that may affect analysis.
## You can select cells from 2-3 methods to analyse.
if(F){
  AvailableData()
  InstallData("pbmcsca")
}
data("pbmcsca")
pbmcsca <- UpdateSeuratObject(pbmcsca)
table(pbmcsca$Method)
# Example: create 3 Seurat objects (selection of other sequencing methods is encouraged)
pbmc_10x_v2 <- pbmcsca[,pbmcsca$Method == "10x Chromium (v2)"]
pbmc_10x_v3 <- pbmcsca[,pbmcsca$Method == "10x Chromium (v3)"]
pbmc_combo <- pbmcsca[,pbmcsca$Method %in% c("10x Chromium (v2)", "10x Chromium (v3)")]
# Analyse pbmc_10x_v2, pbmc_10x_v3 individually, and then analyse pbmc_combo.
# This will allow you to see how technical effects (sequencing method) affects data analysis.

## scRNA-seq dataset2: Villani (2017) ------------------------------------
# This is a small but very classic dataset used by Villani (2017) https://pubmed.ncbi.nlm.nih.gov/28428369/
# to study subtypes of human dendritic cells (DCs) and monocytes 
# A possible analysis is to use visualisation methods and SingleR to analyse this dataset.
# It's prbably interesting to pay more attention to DC5, a new DC subtype identified by Villani (2017).
setwd("~/Documents/GitHub/BINF90004_SingleCell/Practical")
Villani <- readRDS("Villani.rds")
Villani
# Note that Villani is an sce object, not Seurat. And it only contains normalised data:
assay(Villani, "data")[1:10,1:10]
# We can renormalise it using scran
# But first, note that Villani contains some missing values, we can impute them by 0.
sum(is.na(assay(Villani, "data"))) # 29 missing values
assay(Villani, "data")[is.na(assay(Villani, "data"))] <- 0
sum(is.na(assay(Villani, "data"))) # no more
library(scran)
Villani <- computeSumFactors(Villani, assay.type = "data",
                             cluster = quickCluster(Villani, assay.type = "data"))
Villani <- logNormCounts(Villani, assay.type = "data")
assay(Villani, "logcounts")[1:10,1:10]



## scRNA-seq dataset3: ifnb -----------------------------------------------
# This dataset contains control and stimulated immune cells (PBMC)
# It would be interesting to explore the difference between activated and non-activated immune cells
# by some of the visualisation skills we've covered.
if(F){
  InstallData("ifnb")
}
data("ifnb")
ifnb <- UpdateSeuratObject(ifnb)
ifnb$stim %>% table()
ifnb$seurat_annotations %>% table()
# How the cell type annotation above was done is explained here: https://satijalab.org/seurat/archive/v3.2/immune_alignment.html








