

# Single-cell RNA-seq — BINF90004 Bioinformatics Case Studies

Guest lecture and practical for [BINF90004 Bioinformatics Case Studies](https://handbook.unimelb.edu.au/subjects/binf90004), Master of Bioinformatics, University of Melbourne.

### 📖 [Read the materials online →](https://jiadongm.github.io/BINF90004_SingleCell/)

Both practicals and both lecture decks, fully rendered with all code, output and figures. **If your R installation will not work, you can follow the entire session there** — the decisions you are assessed on are answered by thinking, not by running code. Come anyway and sit next to someone whose laptop works.

------------------------------------------------------------------------


## Setup


### 1. Check your R version first

This is the step people skip, and it causes most failures. An old R silently installs old versions of the analysis packages, and the problem surfaces much later as a confusing error.

In RStudio:

``` r
getRversion()
```

| Result           | What to do                                          |
|------------------|-----------------------------------------------------|
| 4.5.0 or newer   | Good. Continue.                                     |
| 4.3.0 – 4.4.x    | Will probably work. Upgrading is still recommended. |
| Older than 4.3.0 | **Stop.** Install a current R before going further. |

To upgrade: get R from <https://cran.r-project.org/>, install it, then **fully restart RStudio** — it does not switch R versions until you do. Update RStudio too (*Help → Check for Updates*) if you have not in the last year.

**Build tools.** A few packages compile from source. If you do not already have a compiler:

- **Windows** — [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version
- **macOS** — run `xcode-select --install` in Terminal
- **Linux** — `build-essential` plus your distribution's R development headers

### 2. Install the packages

Clone or download this repository, open the `.Rproj` or set your working directory to the repo root, then:

``` r
source("setup/install_packages.R")
```

Fourteen packages in two stages, CRAN then Bioconductor. Expect **15–40 minutes** on a fresh machine.

**Long silences are normal.** Compilation prints nothing for minutes at a time. Walls of red text during compilation are usually warnings, not errors — what matters is the verification report in step 4, not the noise along the way.

<details>

<summary>What gets installed, and why</summary>

| Source | Package | Used for |
|------------------------|------------------------|------------------------|
| CRAN | **Seurat** | core single-cell toolkit |
| CRAN | **harmony** | batch correction across donors |
| CRAN | **clustree** | choosing a clustering resolution |
| CRAN | ggplot2, dplyr, tidyr, tibble, patchwork, RColorBrewer, pheatmap | plotting and data handling |
| Bioconductor | **glmGamPoi** | required by `SCTransform(vst.flavor = "v2")` |
| Bioconductor | **edgeR** | pseudobulk aggregation |
| Bioconductor | **limma** | differential expression (voom) |
| Bioconductor | **speckle** | cell type composition testing (`propeller`) |

Deliberately short. Anything that could be pre-computed and shipped with the data instead of installed, has been — including the gene annotation, which removes a large Bioconductor download and a common point of failure.

</details>

### 3. Download the data

The data bundle is on Zenodo:

[**https://doi.org/10.5281/zenodo.21803240**](https://doi.org/10.5281/zenodo.21803240){.uri}

Download **`binf90004_scrnaseq_data_2026.zip`** and unzip it into the root of this repository, so your working folder looks like this:

```         
BINF90004_SingleCell/
├── setup/
├── practicals/
└── data/
    ├── P1_start.rds
    ├── celltype_key.rds
    ├── cluster_labels.rds
    ├── markers_cache.rds
    ├── de_demo.rds
    └── checkpoints/
        └── 02_integrated_clustered.rds
```

About **200 MB**. **Download it at home, not in the classroom** — thirty people pulling 200 MB through lecture-theatre wifi does not work.

Keep the `checkpoints/` folder. It is what lets you rejoin the practical if you fall behind, and there are two points in the session where everyone loads from it regardless of where they have got to.

Several slow, choice-free steps have already been done for you in this bundle (initial filtering, gene annotation, QC metrics, downsampling), so you are not waiting on your laptop during the session.

### 4. Verify

``` r
source("setup/verify_setup.R")
```

You want `ALL CHECKS PASSED`. If not, the script prints a numbered list of problems and a **diagnostic block** at the bottom. Copy that whole block into the LMS thread — it contains everything needed to diagnose the problem and saves a long back-and-forth.

**Then stop.** Do not spend your evening on it. Follow the session on the [rendered materials](https://jiadongm.github.io/BINF90004_SingleCell/) instead — you will not be disadvantaged.


------------------------------------------------------------------------

## The data

Single-**nucleus** RNA-seq of human heart tissue across three developmental stages, from Sim et al. (2021):

| Group  | Donors     | Age         |
|--------|------------|-------------|
| foetal | f1, f2, f3 | 19–20 weeks |
| young  | y1, y2, y3 | 4–14 years  |
| adult  | a1, a2, a3 | 35–42 years |

Nine donors, three per stage. That number — not the number of nuclei — turns out to be the most important number in the analysis.

Note this is **snRNA-seq**, not scRNA-seq: nuclei were isolated rather than whole cells, because adult heart muscle does not dissociate gently. That distinction changes how quality control should be done, and you will meet the consequences within the first twenty minutes.

> Sim CB, Phipson B, Ziemann M, et al. **Sex-Specific Control of Human Heart Maturation by the Progesterone Receptor.** *Circulation*. 2021;143(10):1614–1628. [doi:10.1161/CIRCULATIONAHA.120.051921](https://doi.org/10.1161/CIRCULATIONAHA.120.051921)

Raw counts are publicly available at [Zenodo record 18237749](https://zenodo.org/records/18237749) (CC-BY-4.0). You do not need them for this session — the bundle in step 3 is derived from them.

------------------------------------------------------------------------

## Acknowledgements

The practicals are adapted from the [Advanced Methods for Single-Cell RNA-seq Analysis](https://phipsonlab.github.io/single_cell_workshop/) workshop, quality control, integration, by Belinda Phipson, Melody Jin and Jiadong Mao. 

## License

[MIT](LICENSE)
