################################################################################
# BINF90004 — scRNA-seq session
# Setup verification
#
#   source("verify_setup.R")
#
# Uses base R only, so it will run even if every workshop package is broken.
# It never stops on an error: whatever the state of your machine, it prints a
# report. Copy the DIAGNOSTIC BLOCK at the end into the LMS thread if anything
# fails.
################################################################################

verify_setup <- function(data_dir = "data") {

    pad <- function(x, n) formatC(x, width = -n, flag = " ")
    rule <- function(ch = "-") cat(strrep(ch, 74), "\n", sep = "")
    problems <- character(0)

    cat("\n"); rule("=")
    cat("BINF90004 scRNA-seq — setup verification\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
    rule("="); cat("\n")

    ## ---- 1. R version -------------------------------------------------------
    # The single most common cause of setup failure is an old R. Bioconductor
    # ties its release to the R version, so an out-of-date R silently installs
    # out-of-date (or no) versions of edgeR, limma, speckle and glmGamPoi.
    rv <- getRversion()
    cat("R version\n")
    cat("  installed : ", as.character(rv), "\n", sep = "")
    cat("  required  : 4.3.0 or newer\n")
    cat("  ideal     : 4.5.0 or newer\n")

    if (rv < "4.3.0") {
        cat("  STATUS    : FAIL - too old, please upgrade before installing anything\n")
        problems <- c(problems, "R is older than 4.3.0")
    } else if (rv < "4.5.0") {
        cat("  STATUS    : OK (older than ideal - upgrade if the install misbehaves)\n")
    } else {
        cat("  STATUS    : OK\n")
    }
    cat("\n")

    ## ---- 2. Packages --------------------------------------------------------
    pkgs <- c(
        # CRAN
        "Seurat", "SeuratObject", "harmony", "clustree",
        "ggplot2", "dplyr", "tidyr", "tibble", "patchwork",
        "RColorBrewer", "pheatmap",
        # Bioconductor
        "glmGamPoi", "edgeR", "limma", "speckle"
    )
    source_of <- c(
        Seurat = "CRAN", SeuratObject = "CRAN", harmony = "CRAN",
        clustree = "CRAN", ggplot2 = "CRAN", dplyr = "CRAN", tidyr = "CRAN",
        tibble = "CRAN", patchwork = "CRAN", RColorBrewer = "CRAN",
        pheatmap = "CRAN", glmGamPoi = "Bioc", edgeR = "Bioc",
        limma = "Bioc", speckle = "Bioc"
    )

    cat("Packages\n")
    cat("  ", pad("package", 16), pad("source", 8), pad("version", 12), "status\n", sep = "")
    rule()

    versions <- character(0)
    for (p in pkgs) {
        # Two separate questions: is it installed, and does it actually load?
        # A package can be present but broken (e.g. compiled against a
        # different R), and only loading it will reveal that.
        ver <- tryCatch(as.character(utils::packageVersion(p)),
                        error = function(e) NA_character_)

        if (is.na(ver)) {
            status <- "MISSING"
            problems <- c(problems, paste0(p, " not installed"))
        } else {
            ok <- suppressWarnings(suppressMessages(tryCatch({
                requireNamespace(p, quietly = TRUE)
            }, error = function(e) FALSE)))
            if (isTRUE(ok)) {
                status <- "ok"
            } else {
                status <- "WILL NOT LOAD"
                problems <- c(problems, paste0(p, " installed but fails to load"))
            }
        }

        versions[p] <- if (is.na(ver)) "-" else ver
        cat("  ", pad(p, 16), pad(source_of[[p]], 8),
            pad(versions[p], 12), status, "\n", sep = "")
    }
    cat("\n")

    ## ---- 3. Data ------------------------------------------------------------
    # Sizes are approximate; we only check the file is present and not a
    # truncated download (a failed curl often leaves a small HTML error page).
    # Exactly the files the two vignettes read. `gene_annotation.rds`,
    # `pseudobulk.rds` and `pseudorep_curve.rds` are build artefacts from the
    # prep scripts - they should not go in the student bundle.
    want <- c(`P1_start.rds` = 54,
              `celltype_key.rds` = 0.05,
              `cluster_labels.rds` = 0.0002,
              `markers_cache.rds` = 0.6,
              `de_demo.rds` = 0.0003,
              `checkpoints/02_integrated_clustered.rds` = 154)

    cat("Data files (looking in '", data_dir, "')\n", sep = "")
    if (!dir.exists(data_dir)) {
        cat("  directory not found\n")
        problems <- c(problems, paste0("data directory '", data_dir, "' not found"))
    } else {
        for (f in names(want)) {
            path <- file.path(data_dir, f)
            if (!file.exists(path)) {
                cat("  ", pad(f, 26), "MISSING\n", sep = "")
                problems <- c(problems, paste0(f, " not downloaded"))
            } else {
                mb <- file.size(path) / 1e6
                # Anything under half the expected size is a broken download.
                if (mb < want[[f]] * 0.5) {
                    cat("  ", pad(f, 26),
                        sprintf("TRUNCATED (%.1f MB, expected ~%.1f MB)\n", mb, want[[f]]),
                        sep = "")
                    problems <- c(problems, paste0(f, " looks truncated - re-download"))
                } else {
                    cat("  ", pad(f, 26), sprintf("ok (%.1f MB)\n", mb), sep = "")
                }
            }
        }
    }
    cat("\n")

    ## ---- 4. Verdict ---------------------------------------------------------
    rule("=")
    if (length(problems) == 0) {
        cat("ALL CHECKS PASSED - you are ready for the session.\n")
        rule("=")
        return(invisible(TRUE))
    }

    cat(length(problems), " problem(s) found:\n\n", sep = "")
    for (p in problems) cat("  - ", p, "\n", sep = "")
    cat("\nWhat to do:\n")
    cat("  1. Re-read setup/README.md - most problems are a stale R version.\n")
    cat("  2. Re-run setup/install_packages.R.\n")
    cat("  3. Still stuck? Post the block below in the LMS thread.\n")
    cat("     Do NOT spend more than 30 minutes on this alone. You can follow\n")
    cat("     the whole session on the rendered website without a working R.\n\n")

    rule("=")
    cat("DIAGNOSTIC BLOCK - copy everything between the lines\n")
    rule("=")
    cat("R           : ", as.character(rv), "\n", sep = "")
    cat("Platform    : ", R.version$platform, "\n", sep = "")
    cat("OS          : ", utils::sessionInfo()$running, "\n", sep = "")
    cat("Library     : ", paste(.libPaths(), collapse = " | "), "\n", sep = "")
    cat("Packages    : ",
        paste(names(versions), versions, sep = "=", collapse = ", "), "\n", sep = "")
    cat("Problems    : ", paste(problems, collapse = "; "), "\n", sep = "")
    rule("=")

    invisible(FALSE)
}

# Run immediately when sourced, looking for data/ next to this script's parent.
verify_setup()
