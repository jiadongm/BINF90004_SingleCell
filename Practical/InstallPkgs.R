# For Ubuntu: install these linux packages using ubuntu first
# sudo apt-get install libxml2-dev
# sudo apt-get install libssl-dev
# sudo apt-get install libfontconfig1-dev
# sudo apt-get install libharfbuzz-dev libfribidi-dev
# sudo apt-get install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
# sudo apt-get install
# Prep for devtools
install.packages("curl")
install.packages("usethis")
install.packages("pkgdown")
# R tools
install.packages("devtools")
install.packages("BiocManager")
## Le Cao Lab
devtools::install_github('meiosis97/Sincast@main',subdir = 'pkg')
devtools::install_github("mixOmicsTeam/mixOmics")
## Seurat family
remotes::install_github("satijalab/seurat", "seurat5")
remotes::install_github("satijalab/seurat-data", "seurat5")
remotes::install_github("satijalab/azimuth", "seurat5")
remotes::install_github("satijalab/seurat-wrappers", "seurat5")
remotes::install_github("stuart-lab/signac", "seurat5")
remotes::install_github("mojaveazure/seurat-disk")
## Other comput bio 
BiocManager::install("scran")
BiocManager::install("scuttle")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scater")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
ArchR::installExtraPackages()
## TiddyVerse universe
install.packages("TiddyVerse")
install.packages("markdown")
## Maths and Stats
install.packages("pracma")
install.packages("ks")
install.packages("rARPACK")
install.packages("glmnet")
install.packages("igraph")
install.packages("network")
install.packages("networkD3")
install.packages("philentropy") # distance metrics
## Parallel
install.packages("doParallel")
install.packages("foreach")




