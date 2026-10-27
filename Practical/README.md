# Alternative datasets for the report

This folder is **not** part of the session practicals — those are in [`practicals/`](../practicals/).
It holds the dataset options for your report, carried over from earlier years.

[`alternative_datasets.R`](alternative_datasets.R) loads three datasets and suggests an angle for each:

| Dataset | Source | What it is good for |
|---|---|---|
| `pbmcsca` | SeuratData | PBMCs sequenced with 9 methods — technical/batch effects |
| Villani (2017) | [`Villani.rds`](Villani.rds), this folder | Dendritic cell and monocyte subtypes (look at DC5) |
| `ifnb` | SeuratData | Control vs interferon-stimulated PBMCs |

`Villani.rds` is a `SingleCellExperiment`, not a Seurat object. Run the script from the repo
root or from this folder; it finds the file either way.

**SeuratData is not installed by `setup/install_packages.R`.** It is GitHub-only:

``` r
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")
```

The first `InstallData()` for each dataset downloads it, so do that on a good connection.
