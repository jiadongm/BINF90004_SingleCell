################################################################################
# The pseudoreplication demonstration.  (~4 min, run once)
#
#   Rscript analysis/make_de_demo.R
#
# What did NOT work, and why it is recorded here
# ----------------------------------------------
# First attempt was the textbook contrast: per-cell Wilcoxon vs pseudobulk on
# adult-vs-foetal cardiomyocytes. Result: 3,288 vs 1,569 genes. Only ~2x, which
# is not a lesson - adult and foetal cardiomyocytes really are enormously
# different, so both methods find plenty.
#
# Second attempt was a power curve, subsampling cells and running both methods:
#
#     n/group   cells   wilcoxon   pseudobulk
#          50     100        330           81
#         100     200        785          689
#         200     400      1,427        1,038
#         333     666      2,175        1,340
#
# Both climb. That is expected and unremarkable - more cells means a better
# estimate of each donor's profile, which helps pseudobulk too. This curve does
# NOT demonstrate pseudoreplication, and presenting it as if it did would be
# teaching a real error with a fake example.
#
# What actually demonstrates it
# -----------------------------
# Compare two donors from the SAME developmental stage. There is no
# developmental contrast to find, so every "significant" gene is donor
# variation - genetics, tissue handling, library prep - being reported as
# biology. Then note that pseudobulk cannot even be run on that comparison
# (n = 1 vs 1), which is the honest answer: you cannot test it. The per-cell
# test manufactures thousands of degrees of freedom from a single donor and
# returns an answer anyway.
#
# All comparisons below use the same number of cells per side, so the counts
# are directly comparable.
################################################################################

suppressPackageStartupMessages({
    library(Seurat); library(edgeR); library(limma); library(dplyr)
})
set.seed(42)
options(future.globals.maxSize = 8 * 1024^3)

OUT   <- "/Users/jmao1/Library/CloudStorage/Dropbox/Teaching/BINF90004BinfCaseStudies/BINF_scRNA-seq_2026/dist"
N_PER <- 300   # cells per side, every comparison

seu <- readRDS(file.path(OUT, "checkpoints", "02_integrated_clustered.rds"))
lab <- readRDS(file.path(OUT, "cluster_labels.rds"))
seu$celltype <- unname(lab[as.character(seu$seurat_clusters)])
cm <- subset(seu, subset = celltype == "Cardiomyocytes")
cat("cardiomyocytes per donor:\n"); print(table(cm$sample))

wilcox_n <- function(obj, id_col, a, b) {
    Idents(obj) <- id_col
    r <- FindMarkers(obj, ident.1 = a, ident.2 = b, assay = "SCT", verbose = FALSE)
    sum(r$p_val_adj < 0.05)
}
take <- function(obj, col, vals, n) {
    cells <- unlist(lapply(vals, function(v) {
        x <- colnames(obj)[obj[[col, drop = TRUE]] == v]
        sample(x, min(n, length(x)))
    }))
    subset(obj, cells = cells)
}

## ---- 1. The real contrast, matched at N_PER per side -----------------------
cat("\n1. Real contrast: adult vs foetal\n")
ad <- colnames(cm)[cm$group == "adult"]; fe <- colnames(cm)[cm$group == "fetal"]
real <- subset(cm, cells = c(sample(ad, min(N_PER, length(ad))),
                             sample(fe, min(N_PER, length(fe)))))
n_real <- wilcox_n(real, "group", "adult", "fetal")
cat("   adult vs fetal:", n_real, "genes\n")

## ---- 2. Same-stage donor pairs: there is nothing to find -------------------
cat("\n2. Same-stage donor pairs (no developmental contrast exists)\n")
pairs <- list(c("f1", "f2"), c("f2", "f3"), c("f1", "f3"))
null_res <- lapply(pairs, function(p) {
    sub <- take(cm, "sample", p, N_PER)
    n   <- wilcox_n(sub, "sample", p[1], p[2])
    cat(sprintf("   %s vs %s: %d genes\n", p[1], p[2], n))
    data.frame(comparison = paste(p, collapse = " vs "),
               kind = "same stage (foetal)", n_sig = n)
})

## ---- 3. Pseudobulk on the real contrast ------------------------------------
cat("\n3. Pseudobulk, adult vs foetal (all cells, 6 donors)\n")
cm2 <- subset(cm, subset = group %in% c("adult", "fetal"))
cm2$Celltype <- "CM"; cm2$Sample <- cm2$sample
DefaultAssay(cm2) <- "RNA"
pb <- Seurat2PB(cm2, sample = "Sample", cluster = "Celltype")
g  <- c(a = "adult", f = "fetal")[sub("[0-9]+$", "", as.character(pb$samples$sample))]
keep <- filterByExpr(pb, group = g)
pb  <- normLibSizes(pb[keep, , keep.lib.sizes = FALSE])
des <- model.matrix(~ 0 + factor(g)); colnames(des) <- levels(factor(g))
fit <- eBayes(contrasts.fit(lmFit(voom(pb, des), des),
                            makeContrasts(adult - fetal, levels = des)))
n_pb <- sum(decideTests(fit) != 0)
cat("   pseudobulk:", n_pb, "genes (n = 3 donors per group)\n")

## ---- Assemble ---------------------------------------------------------------
demo <- bind_rows(
    data.frame(comparison = "adult vs fetal", kind = "real contrast", n_sig = n_real),
    bind_rows(null_res)
)
out <- list(n_per_side = N_PER, wilcoxon = demo,
            pseudobulk_real = n_pb,
            curve = readRDS(file.path(OUT, "pseudorep_curve.rds")))
saveRDS(out, file.path(OUT, "de_demo.rds"))

cat("\n=== summary ===\n"); print(demo)
cat("\nnoise floor: same-stage donor pairs reach ",
    round(100 * max(sapply(null_res, `[[`, "n_sig")) / n_real),
    "% of the 'real' result\n", sep = "")
cat("\nSaved:", file.path(OUT, "de_demo.rds"), "\n")
