################################################################################
# BINF90004 — scRNA-seq session
# Package installation
#
#   source("install_packages.R")
#
# Expect 15-40 minutes on a first run. Several packages compile from source,
# which is slow but normal - long silences are not a crash.
#
# BEFORE YOU RUN THIS: check your R version is 4.3.0 or newer (4.5+ preferred).
# Installing onto an old R is the single most common way this goes wrong, and
# the failure appears much later, in a confusing place.
################################################################################

if (getRversion() < "4.3.0") {
    stop(
        "\n\nYour R is ", getRversion(), ", which is too old.\n",
        "Install R 4.5 or newer from https://cran.r-project.org/ first,\n",
        "restart RStudio, then run this script again.\n",
        "See setup/README.md, section 1.\n\n",
        call. = FALSE
    )
}

message("R ", getRversion(), " - ok\n")

## ---- Stage 1: CRAN -----------------------------------------------------------
cran_pkgs <- c(
    "Seurat",        # core single-cell toolkit (compiles; the slow one)
    "harmony",       # batch correction
    "clustree",      # clustering resolution diagnostics
    "ggplot2", "dplyr", "tidyr", "tibble", "patchwork",
    "RColorBrewer", "pheatmap",
    "BiocManager"    # gateway for stage 2
)

missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1),
                                  quietly = TRUE)]

if (length(missing_cran)) {
    message("Stage 1/2 - installing from CRAN: ",
            paste(missing_cran, collapse = ", "))
    install.packages(missing_cran)
} else {
    message("Stage 1/2 - all CRAN packages already present")
}

## ---- Stage 2: Bioconductor ---------------------------------------------------
# Only four, deliberately. The workshop this is adapted from also uses
# org.Hs.eg.db and AnnotationDbi for gene annotation; we ship the annotation
# pre-computed instead, which removes a large download and a common failure.
bioc_pkgs <- c(
    "glmGamPoi",     # required by SCTransform(vst.flavor = "v2")
    "edgeR",         # pseudobulk aggregation + DE
    "limma",         # linear modelling / voom
    "speckle"        # propeller composition testing
)

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1),
                                  quietly = TRUE)]

if (length(missing_bioc)) {
    message("\nStage 2/2 - installing from Bioconductor: ",
            paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
} else {
    message("\nStage 2/2 - all Bioconductor packages already present")
}

## ---- Report ------------------------------------------------------------------
message("\nInstallation finished. Verifying...\n")

if (file.exists("verify_setup.R")) {
    source("verify_setup.R")
} else {
    message("Now run:  source(\"verify_setup.R\")")
}
