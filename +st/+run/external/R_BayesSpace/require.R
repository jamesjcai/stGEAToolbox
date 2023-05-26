if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BayesSpace")

library(SingleCellExperiment)
#library(ggplot2)
library(BayesSpace)
# https://edward130603.github.io/BayesSpace/articles/BayesSpace.html
library(rhdf5)
