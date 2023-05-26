setwd("C:\\Users\\jcai\\Documents\\GitHub\\spatial_transcriptomics\\external\\R_BayesSpace")


library(Seurat)
library(Matrix)
countsmatrix <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
geneList <- make.unique(countsmatrix[,1])[-1]
bcList <- countsmatrix[1,][-1]
countMatrix <- Matrix(as.matrix(countsmatrix[-1,-1]))
rownames(countMatrix) <- geneList
colnames(countMatrix) <- bcList
sce1 <- CreateSeuratObject(countMatrix)
# sce1 <- NormalizeData(countMatrix)

library(SingleCellExperiment)
sce2 <- SingleCellExperiment(assays=list(counts=as(countMatrix, "dgCMatrix")),
                            rowData=geneList,
                            colData=bcList)


sce2 <- SingleCellExperiment(assays=list(counts=as(countMatrix, "dgCMatrix")),
                            rowData=rownames(countMatrix),
                            colData=colnames(countMatrix))




counts <- read.csv("in.csv",
                   row.names=1, check.names=F, stringsAsFactors=FALSE))

assays=list(counts=counts),

sce2 <- SingleCellExperiment(assays=list(counts=as(counts, "dgCMatrix")),
                            rowData=geneList,
                            colData=bcList)

library(Matrix)
rowData <- read.csv("path/to/rowData.csv", stringsAsFactors=FALSE)
colData <- read.csv("path/to/colData.csv", stringsAsFactors=FALSE, row.names=1)
counts <- read.csv("path/to/counts.csv.gz",
                   row.names=1, check.names=F, stringsAsFactors=FALSE))


sce <- SingleCellExperiment(assays=list(counts=as(counts, "dgCMatrix")),
                            rowData=rowData,
                            colData=colData)