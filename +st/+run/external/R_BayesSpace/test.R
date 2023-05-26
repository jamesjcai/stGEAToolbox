


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BayesSpace")

library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)

# https://edward130603.github.io/BayesSpace/articles/BayesSpace.html

melanoma <- getRDS(dataset="2018_thrane_melanoma", sample="ST_mel1_rep2")
set.seed(102)
melanoma <- spatialPreprocess(melanoma, platform="ST", 
                              n.PCs=7, n.HVGs=2000, log.normalize=FALSE)

melanoma <- qTune(melanoma, qs=seq(2, 10), platform="ST", d=7)
qPlot(melanoma)

set.seed(149)
melanoma <- spatialCluster(melanoma, q=5, platform="ST", d=7,
                           init.method="mclust", model="t", gamma=2,
                           nrep=1000, burn.in=100,
                           save.chain=TRUE)
head(colData(melanoma))

clusterPlot(melanoma)


