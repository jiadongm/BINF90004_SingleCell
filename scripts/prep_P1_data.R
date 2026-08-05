################################################################################
# Builds everything P1 ships to students.
#
#   Rscript analysis/prep_P1_data.R
#
# Run from the workshop repo root (where data/ lives). Outputs land in
# BINF_scRNA-seq_2026/dist/.
#
# Design notes
# ------------
# * The gene annotation is computed ONCE here and shipped as a plain data
#   frame, so students never install org.Hs.eg.db / AnnotationDbi. That is two
#   fewer Bioconductor packages and a ~100 MB download removed from setup.
#
# * QC metrics (percent.mt, percent.ribo) are computed BEFORE gene filtering,
#   because gene filtering removes the MT- genes that percent.mt depends on.
#   The metrics live in metadata, so Decision 1 still works on the shipped
#   object.
#
# * Cells are downsampled but NOT filtered. Filtering is Decision 1 - if we
#   shipped a pre-filtered object there would be nothing to decide.
################################################################################

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
})

set.seed(42)
options(future.globals.maxSize = 8 * 1024^3)

OUT      <- "/Users/jmao1/Library/CloudStorage/Dropbox/Teaching/BINF90004BinfCaseStudies/BINF_scRNA-seq_2026/dist"

# Pre-QC cell count. Compute is not the constraint (SCTransform on 10.5k took
# 22 s on an M-series Mac); the constraint is the size of the file students
# have to download. 8,000 keeps the rarest published cell type (smooth muscle,
# 0.8% of the data) at ~60 nuclei, which is still enough to cluster.
N_TARGET <- 8000

dir.create(file.path(OUT, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
tic <- function() assign(".t0", Sys.time(), envir = .GlobalEnv)
toc <- function(label) cat(sprintf("  [%-28s %5.1f s]\n", paste0(label, "]"),
                                   as.numeric(difftime(Sys.time(), .t0, units = "secs"))))

## ---- 1. Load ---------------------------------------------------------------
cat("1. Loading raw data\n"); tic()
counts   <- readRDS("data/heart-counts.Rds")
cellinfo <- readRDS("data/cellinfo_updated.Rds")
toc("readRDS")

cat("2. Creating Seurat object\n"); tic()
seu <- CreateSeuratObject(counts = counts, project = "heart",
                          min.cells = 3, min.features = 200)
m <- cellinfo[match(colnames(seu), cellinfo$CellID), ]
stopifnot(all(colnames(seu) == m$CellID))
seu$sample <- m$Sample
seu$group  <- m$Group
seu$sex    <- m$Sex
# Ground-truth labels from Sim et al. Held aside deliberately: if they rode
# along in the shipped object, a student could `table(seu$celltype)` their way
# past the Decision 4 marker exercise in P2. Shipped as a separate key file
# that P2 loads only at the reveal.
celltype_key <- setNames(as.character(m$Celltype), colnames(seu))
rm(counts, cellinfo, m); invisible(gc())
toc("CreateSeuratObject")

## ---- 2. Gene annotation (computed once, shipped) ---------------------------
cat("3. Building gene annotation\n"); tic()
sym <- rownames(seu)
ann <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = sym, keytype = "SYMBOL",
    columns = c("ENTREZID", "GENENAME", "CHR")))
ann <- ann[!duplicated(ann$SYMBOL), ]
ann <- ann[match(sym, ann$SYMBOL), ]
rownames(ann) <- NULL
saveRDS(ann, file.path(OUT, "gene_annotation.rds"))
toc("annotation")
cat("   annotated:", round(100 * mean(!is.na(ann$GENENAME)), 1), "%\n")

## ---- 3. QC metrics (before gene filtering!) --------------------------------
cat("4. QC metrics\n"); tic()
seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
toc("PercentageFeatureSet")

## ---- 4. Gene filtering (plumbing, not a decision) --------------------------
# Same four categories as the workshop's Module 1: mitochondrial (MT- prefix
# and nuclear-encoded), ribosomal, no Entrez ID, and sex chromosomes. Sex
# chromosomes go because sex is unbalanced across developmental groups
# (2M/1F, 2M/1F, 1M/2F), which would otherwise produce spurious DE.
cat("5. Gene filtering\n"); tic()
drop <- unique(c(
    grep("^MT-", sym),
    grep("mitochondrial", ann$GENENAME, ignore.case = TRUE),
    grep("^RP[SL][0-9]", sym),
    grep("ribosomal", ann$GENENAME, ignore.case = TRUE),
    which(is.na(ann$ENTREZID)),
    which(ann$CHR %in% c("X", "Y"))
))
seu <- subset(seu, features = setdiff(sym, sym[drop]))
toc("gene filter")
cat("   genes:", length(sym), "->", nrow(seu), "\n")

## ---- 5. Stratified downsample ----------------------------------------------
# Proportional to sample size, so the per-donor QC variation that Decision 1
# turns on (a3 loses 19.3% at percent.mt < 5, f1 loses 1.5%) is preserved.
cat("6. Downsampling to", N_TARGET, "\n"); tic()
n_by <- table(seu$sample)
take <- round(N_TARGET * n_by / sum(n_by))
keep <- unlist(lapply(names(take), function(s) {
    cells <- colnames(seu)[seu$sample == s]
    sample(cells, min(take[[s]], length(cells)))
}))
seu <- subset(seu, cells = keep)
toc("downsample")
print(table(seu$sample))
# Note the imbalance: adult donors genuinely have far fewer nuclei than foetal
# (a3 = 1,681 vs f2 = 10,948 in the full data). Sampling proportionally keeps
# that, because it is real and it matters for the composition test in P2.

celltype_key <- celltype_key[colnames(seu)]
saveRDS(celltype_key, file.path(OUT, "celltype_key.rds"))

saveRDS(seu, file.path(OUT, "P1_start.rds"))
cat("   P1_start.rds:", round(file.size(file.path(OUT, "P1_start.rds")) / 1e6, 1), "MB\n")

## ---- 6. Apply reference QC thresholds --------------------------------------
# Not saved as a separate checkpoint file: the P1 sync point re-reads
# P1_start.rds and applies these same four lines, which takes ~10 s and saves
# students a ~100 MB download.
cat("7. Reference QC filter\n"); tic()
seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA > 2500 &
                            nCount_RNA < 40000 & percent.mt < 20)
toc("QC filter")
cat("   cells after QC:", ncol(seu), "\n")

## ---- 7. Checkpoint 02: normalised, integrated, clustered -------------------
cat("8. SCTransform\n"); tic()
seu <- SCTransform(seu, vst.flavor = "v2", seed.use = 42, verbose = FALSE)
toc("SCTransform")

cat("9. PCA\n"); tic()
seu <- RunPCA(seu, seed.use = 42, verbose = FALSE)
toc("RunPCA")

N_DIMS <- 20
cat("10. Harmony\n"); tic()
seu <- harmony::RunHarmony(seu, group.by.vars = "sample",
                           reduction.use = "pca", dims.use = 1:N_DIMS,
                           verbose = FALSE)
toc("RunHarmony")

cat("11. UMAP\n"); tic()
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:N_DIMS,
               seed.use = 42, verbose = FALSE)
toc("RunUMAP")

cat("12. Neighbours\n"); tic()
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:N_DIMS, verbose = FALSE)
toc("FindNeighbors")

cat("13. Resolution sweep\n"); tic()
RESOLUTIONS <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8)
for (r in RESOLUTIONS) seu <- FindClusters(seu, resolution = r, verbose = FALSE)
toc("FindClusters x6")
for (r in RESOLUTIONS)
    cat(sprintf("   res %.1f -> %2d clusters\n", r,
                length(unique(seu@meta.data[[paste0("SCT_snn_res.", r)]]))))

Idents(seu) <- "SCT_snn_res.0.4"
seu$seurat_clusters <- Idents(seu)

## ---- 8. Slim down before shipping ------------------------------------------
# The first build of this file was 473 MB, almost all of it two things nobody
# downstream needs:
#   * SCT scale.data - 3,000 x n dense residuals (254 MB). Only DoHeatmap uses
#     it, and P2 can rebuild it for a handful of genes in seconds.
#   * the SNN/NN graphs - only needed to re-run clustering, which P2 does not.
cat("14. Slimming\n")
before_mb <- as.numeric(object.size(seu)) / 1e6
seu[["SCT"]]$scale.data <- matrix(numeric(0), nrow = 0, ncol = 0)
for (g in Graphs(seu)) seu[[g]] <- NULL
cat(sprintf("   in memory: %.0f -> %.0f MB\n", before_mb,
            as.numeric(object.size(seu)) / 1e6))

f <- file.path(OUT, "checkpoints", "02_integrated_clustered.rds")
saveRDS(seu, f)
cat("   02_integrated_clustered.rds:", round(file.size(f) / 1e6, 1), "MB\n")

cat("\nTotal student download:",
    round(sum(file.size(list.files(OUT, recursive = TRUE, full.names = TRUE))) / 1e6, 1),
    "MB\n")
cat("\nDone.\n")
