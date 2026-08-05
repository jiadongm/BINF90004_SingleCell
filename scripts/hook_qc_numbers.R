################################################################################
# Generates the real numbers behind the AI-hook QC slide (PLAN.md §3.1 beat 4).
#
# Point of this script: the hook must be built on what this dataset actually
# does, not on a plausible-sounding story. The first version of the plan
# asserted that a naive `percent.mt < 5` filter would preferentially delete
# cardiomyocytes, because cardiac tissue is mitochondria-rich. That is a
# reasonable heuristic and it is WRONG here - see the results below.
#
# Run from the workshop repo root (where data/ lives).
################################################################################

suppressPackageStartupMessages(library(Matrix))

counts   <- readRDS("data/heart-counts.Rds")
cellinfo <- readRDS("data/cellinfo_updated.Rds")

## Mirror CreateSeuratObject(min.cells = 3, min.features = 200) --------------
keep_g <- Matrix::rowSums(counts > 0) >= 3
keep_c <- Matrix::colSums(counts > 0) >= 200
x <- counts[keep_g, keep_c]

mt   <- grep("^MT-", rownames(x))          # all 13 are present in this file
pmt  <- Matrix::colSums(x[mt, , drop = FALSE]) / Matrix::colSums(x) * 100

meta <- cellinfo[match(names(pmt), cellinfo$CellID), ]
df   <- data.frame(pmt = as.numeric(pmt), ct = meta$Celltype,
                   grp = meta$Group, samp = meta$Sample)
df   <- df[!is.na(df$ct), ]

## 1. The distribution is low, because these are nuclei ---------------------
# Median 0.25%. Nuclei preps wash away the cytoplasm, and mitochondria go
# with it. This is why the cardiac-tissue intuition fails.
print(summary(df$pmt))

## 2. What each threshold costs ---------------------------------------------
for (t in c(1, 5, 10, 20)) {
    cat(sprintf("percent.mt < %-2d  removes %5.1f%% of nuclei\n",
                t, 100 * mean(df$pmt >= t)))
}
# < 5  -> 6.5%   (the pbmc3k default an agent will reach for)
# < 20 -> 1.0%   (what the workshop uses)
# ~3,000 nuclei separate the two.

## 3. The loss is not uniform across cell types -----------------------------
tot  <- table(df$ct)
lost <- table(df$ct[df$pmt >= 5])
res  <- data.frame(celltype = names(tot),
                   n        = as.integer(tot),
                   pct_lost = round(100 * as.integer(lost[names(tot)]) /
                                        as.integer(tot), 1))
print(res[order(-res$pct_lost), ], row.names = FALSE)
# Endothelial 15.6%, Immune 13.4%  vs  Cardiomyocytes 6.0%.
# The rarest populations are hit hardest - exactly the ones you have least
# power to detect. Cardiomyocytes have the LOWEST median MT (0.16%), the
# opposite of the tissue-level intuition.

## 4. And it is mostly a per-donor effect -----------------------------------
print(round(tapply(df$pmt, df$samp, function(z) 100 * mean(z >= 5)), 1))
# a3 = 19.3%   f1 = 1.5%   -> a 13x spread across donors.
# One fixed global threshold lands completely differently on each sample.
# With n = 3 per group, a single bad donor can move a group's composition,
# which is noise injected straight into the propeller test in P2.

## 5. Within-cell-type by group (for completeness) --------------------------
print(tapply(df$pmt, list(df$ct, df$grp), function(z) round(100 * mean(z >= 5), 1)))
# Group differences exist but are inconsistent in direction - do NOT claim a
# clean "adults lose more" story on a slide. The donor-level story is the
# defensible one.
