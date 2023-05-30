library(SingleCellExperiment)
library(SpaTalk)
#setwd("C:\\Users\\jcai\\Documents\\GitHub\\spatial_transcriptomics\\external\\R_SpaTalk")
#setwd("U:\\GitHub\\spatial_transcriptomics\\external\\R_SpaTalk")
setwd("D:\\GitHub\\stGEAToolbox\\+st\\+run\\external\\R_SpaTalk")

library(Seurat)
library(Matrix)


counts <- h5read(file = "input.mat", name = "/X")
countMatrix <- Matrix(as.matrix(counts))
rownames(countMatrix) <- paste0("G", seq_len(nrow(countMatrix)))
colnames(countMatrix) <- paste0("C", seq_len(ncol(countMatrix)))


#countsmatrix <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
#geneList <- make.unique(countsmatrix[,1])[-1]
#bcList <- countsmatrix[1,][-1]
#countMatrix <- Matrix(as.matrix(countsmatrix[-1,-1]))
#rownames(countMatrix) <- geneList
#colnames(countMatrix) <- bcList


starmap_data <- rev_gene(data = as.matrix(countMatrix),
                         data_type = "count",
                         species = "Mouse",
                         geneinfo = rownames(countMatrix))
