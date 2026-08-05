################################################################################
# Builds what P2 ships. Run after prep_P1_data.R.
#
#   Rscript analysis/prep_P2_data.R
#
# Deliberately produces only SMALL files. Checkpoint 03 is not shipped as an
# object: it would be checkpoint 02 plus one metadata column, i.e. another
# 79 MB download for ~2 kB of information. Instead we ship the cluster ->
# cell type map, and the P2 sync point applies it to the object students
# already have.
################################################################################

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
})

set.seed(42)
options(future.globals.maxSize = 8 * 1024^3)

OUT <- "/Users/jmao1/Library/CloudStorage/Dropbox/Teaching/BINF90004BinfCaseStudies/BINF_scRNA-seq_2026/dist"
tic <- function() assign(".t0", Sys.time(), envir = .GlobalEnv)
toc <- function(l) cat(sprintf("  [%-22s %5.1f s]\n", paste0(l, "]"),
                               as.numeric(difftime(Sys.time(), .t0, units = "secs"))))

cat("1. Load\n"); tic()
seu <- readRDS(file.path(OUT, "checkpoints", "02_integrated_clustered.rds"))
key <- readRDS(file.path(OUT, "celltype_key.rds"))
toc("readRDS")

## ---- Markers ---------------------------------------------------------------
cat("2. PrepSCTFindMarkers\n"); tic()
seu <- PrepSCTFindMarkers(seu, verbose = FALSE)
toc("PrepSCTFindMarkers")

cat("3. FindAllMarkers\n"); tic()
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.5, verbose = FALSE)
toc("FindAllMarkers")
saveRDS(markers, file.path(OUT, "markers_cache.rds"))
cat("   rows:", nrow(markers), " file:",
    round(file.size(file.path(OUT, "markers_cache.rds")) / 1e6, 2), "MB\n")

## ---- Reference cluster labels ----------------------------------------------
# Majority vote of the published Sim et al. labels within each cluster. This
# stands in for the manual annotation students do in Decision 4 - it is what
# a careful annotation converges on, and it means the composition and DE
# sections downstream run on sensible labels rather than on someone's
# half-finished guess.
cat("4. Cluster labels\n")
ct  <- key[colnames(seu)]
tb  <- table(seu$seurat_clusters, ct)
lab <- setNames(colnames(tb)[apply(tb, 1, which.max)], rownames(tb))
purity <- round(100 * apply(tb, 1, max) / rowSums(tb), 1)

print(data.frame(cluster = names(lab), label = unname(lab),
                 n = as.integer(rowSums(tb)), purity = unname(purity)),
      row.names = FALSE)
cat("\n   mean purity:", round(mean(purity), 1), "%\n")
cat("   clusters:", length(lab), " -> distinct cell types:", length(unique(lab)), "\n")

saveRDS(lab, file.path(OUT, "cluster_labels.rds"))

## ---- Sanity-check the downstream steps run ---------------------------------
seu$celltype <- unname(lab[as.character(seu$seurat_clusters)])

cat("\n5. propeller check\n"); tic()
suppressPackageStartupMessages(library(speckle))
grp <- setNames(as.character(seu$group), seu$sample)[!duplicated(seu$sample)]
pr <- propeller(clusters = seu$celltype, sample = seu$sample,
                group = seu$group)
toc("propeller")
print(pr)

cat("\n6. pseudobulk check\n"); tic()
suppressPackageStartupMessages({library(edgeR); library(limma)})
seu$Celltype <- seu$celltype
seu$Sample   <- seu$sample
DefaultAssay(seu) <- "RNA"
pb <- Seurat2PB(seu, sample = "Sample", cluster = "Celltype")
toc("Seurat2PB")
cat("   pseudobulk samples:", ncol(pb), " genes:", nrow(pb), "\n")
saveRDS(pb, file.path(OUT, "pseudobulk.rds"))
cat("   pseudobulk.rds:",
    round(file.size(file.path(OUT, "pseudobulk.rds")) / 1e6, 2), "MB\n")

cat("\nDone.\n")
