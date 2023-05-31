if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repo = "http://cran.rstudio.com/")

if (!require("SingleCellExperiment", quietly = TRUE))
    BiocManager::install("SingleCellExperiment")

if (!require("BayesSpace", quietly = TRUE))
    BiocManager::install("BayesSpace")

if (!require("rhdf5", quietly = TRUE))
    BiocManager::install("rhdf5")


#library(SingleCellExperiment)
#library(ggplot2)
#library(BayesSpace)
# https://edward130603.github.io/BayesSpace/articles/BayesSpace.html
#library(rhdf5)
