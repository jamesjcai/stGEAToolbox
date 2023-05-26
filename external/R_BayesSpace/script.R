library(SingleCellExperiment)
library(BayesSpace)
library(Seurat)
library(Matrix)
library(rhdf5)

# library(ggplot2)
# setwd("C:\\Users\\jcai\\Documents\\GitHub\\spatial_transcriptomics\\external\\R_BayesSpace")
# setwd("U:\\GitHub\\spatial_transcriptomics\\external\\R_BayesSpace")

usingmat<-TRUE

if (usingmat) {
    counts <- h5read(file = "input.mat", name = "/X")
    countMatrix <- Matrix(as.matrix(counts))
    rownames(countMatrix) <- paste0("G", seq_len(nrow(countMatrix)))
    # colnames(countMatrix) <- paste0("C", seq_len(ncol(countMatrix)))    
} else {
    counts <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
    geneList <- make.unique(counts[,1])[-1]
    countMatrix <- Matrix(as.matrix(counts[-1,-1]))
    rownames(countMatrix) <- geneList
    # bcList <- countsmatrix[1,][-1]
    # colnames(countMatrix) <- bcList
}

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

sce <- qTune(sce, qs=seq(2, 10), platform="Visium", d=7)
qPlot(sce)


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


markers <- c("PMEL", "CD2", "CD19", "COL1A1")
sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                     feature_names=markers,
                                     nrounds=0)
X<-logcounts(sce.enhanced)[markers, 1:5]
h5createFile("output.h5")
h5write(as.matrix(X), "output.h5","/X")

# rowData(sce.enhanced)[markers, ]
# featurePlot(sce.enhanced, "PMEL")

