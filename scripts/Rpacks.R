# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("pcaExplorer")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("BiocStyle")

install.packages("shinyBS")
# install.packages("dqshiny")  # No longer available? Get from archive.
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/dqshiny/dqshiny_0.0.4.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
install.packages("remotes")
remotes::install_github("daqana/dqshiny")
install.packages("ggseqlogo")

BiocManager::install("DEP")
install.packages("plotly")

install.packages("sna")
install.packages("network")
packageurl <- "https://cran.r-project.org/src/contrib/Archive/geomnet/geomnet_0.3.1.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
install.packages("ggpubr")
install.packages("ggvenn")
install.packages("ggdendro")
install.packages("heatmap3")