# For linux OS: might neet to install these linux packages
# sudo apt-get install libxml2-dev
# sudo apt-get install libssl-dev
# sudo apt-get install libfontconfig1-dev
# sudo apt-get install libharfbuzz-dev libfribidi-dev
# sudo apt-get install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
# Prep for devtools (you do not need these if you have already had devtools installed)
if(F){
  install.packages("curl")
  install.packages("usethis")
  install.packages("pkgdown")
}
# R tools (you most likely have already installed these)
install.packages("devtools")
install.packages("BiocManager")
## Seurat family 
remotes::install_github("satijalab/seurat", "seurat5")
remotes::install_github("satijalab/seurat-data", "seurat5")
remotes::install_github("satijalab/seurat-wrappers", "seurat5")
## Other comput bio 
BiocManager::install("scran")
BiocManager::install("scuttle")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scater")
## TiddyVerse universe
install.packages("dplyr", "patchwork", "ggplot2", "gridExtra",
                 "ggpubr", "tidyr", "tibble")
## Maths and Stats
install.packages("rARPACK")
install.packages("glmnet")
install.packages("networkD3")




