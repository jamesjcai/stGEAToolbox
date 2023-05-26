library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
# setwd("C:\\Users\\jcai\\Documents\\GitHub\\spatial_transcriptomics\\external\\R_BayesSpace")
# setwd("U:\\GitHub\\spatial_transcriptomics\\external\\R_BayesSpace")

library(Seurat)
library(Matrix)
counts <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
geneList <- make.unique(counts[,1])[-1]
# bcList <- countsmatrix[1,][-1]
countMatrix <- Matrix(as.matrix(counts[-1,-1]))
rownames(countMatrix) <- geneList
# colnames(countMatrix) <- bcList

colData <- read.csv("positions.csv", header=FALSE)
colnames(colData) <- c("spot", "in_tissue", "row", "col", "imagerow", "imagecol")
rownames(colData) <- colData$spot
# colData <- colData[colData$in_tissue > 0, ]

colnames(countMatrix) <- colData$spot

#sce <- SingleCellExperiment(assays=list(counts=as(countMatrix, "dgCMatrix")),
#                            rowData=rownames(countMatrix),
#                            colData=colData)

sce <- SingleCellExperiment(assays=list(counts=countMatrix),
                            rowData=rownames(countMatrix),
                            colData=colData)

libsizes <- colSums(countMatrix)
size.factors <- libsizes/mean(libsizes)
logcounts(sce) <- log2(t(t(countMatrix)/size.factors) + 1)

#metadata(sce)$BayesSpace.data <- list()
#metadata(sce)$BayesSpace.data$platform <- "Visium"
#metadata(sce)$BayesSpace.data$is.enhanced <- FALSE


sce <- spatialPreprocess(sce, platform="Visium", 
                              n.PCs=7, n.HVGs=2000, log.normalize=FALSE)

#sce <- qTune(sce, qs=seq(2, 10), platform="Visium", d=7)
#qPlot(sce)


sce <- spatialCluster(sce, q=6, platform="Visium", d=7,
                           init.method="mclust", model="t", gamma=2,
                           nrep=1000, burn.in=100,
                           save.chain=FALSE)

sce.enhanced <- spatialEnhance(sce, q=6, platform="Visium", d=7,
                                model="t", gamma=2,
                                jitter_prior=0.3, jitter_scale=3.5,
                                nrep=1000, burn.in=100,
                                save.chain=TRUE)

write.csv(colData(sce.enhanced), "positions_enhanced.csv");
