#setwd("C:\\Users\\jcai\\Documents\\GitHub\\spatial_transcriptomics\\external\\R_BayesSpace")
setwd("U:\\GitHub\\spatial_transcriptomics\\example_data\\DD073R_A1_processed")
# setwd("C:\\Users\\jcai\\Documents\\GitHub\\spatial_transcriptomics\\example_data\\DD073R_A1_processed")

colData <- read.csv(file.path("spatial", "tissue_positions_list.csv"), header=FALSE)
colnames(colData) <- c("spot", "in_tissue", "row", "col", "imagerow", "imagecol")
rownames(colData) <- colData$spot
colData <- colData[colData$in_tissue > 0, ]


library(Seurat)
library(Matrix)
library(hdf5r)
counts <- Read10X_h5("filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)



library(SingleCellExperiment)
counts <- counts[, rownames(colData)]

sce <- SingleCellExperiment(assays=list(counts=counts),                                
                        rowData=rownames(counts),
                        colData=colData)

libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(sce) <- log2(t(t(counts)/size.factors) + 1)


library(ggplot2)
library(BayesSpace)                                
sce <- spatialPreprocess(sce, platform="Visium", n.PCs=7, n.HVGs=2000, log.normalize=FALSE)


sce <- qTune(sce, qs=seq(2, 10), platform="Visium", d=7)
qPlot(sce)

# a <- exampleSCE()

set.seed(149)
sce <- spatialCluster(sce, q=6, platform="Visium", d=7,
                           init.method="mclust", model="t", gamma=2,
                           nrep=1000, burn.in=100,
                           save.chain=TRUE)
head(colData(sce))

clusterPlot(sce)


clusterPlot(sce, palette=c("purple", "red", "blue", "yellow","gray","green"), color="black") +
  theme_bw() +
  xlab("Column") +
  ylab("Row") +
  labs(fill="BayesSpace\ncluster", title="Spatial clustering of ST_mel1_rep2")

sce.enhanced <- spatialEnhance(sce, q=6, platform="Visium", d=7,
                                model="t", gamma=2,
                                jitter_prior=0.3, jitter_scale=3.5,
                                nrep=1000, burn.in=100,
                                save.chain=TRUE)

clusterPlot(sce.enhanced)


markers <- c("PMEL", "CD2", "CD19", "COL1A1", "ROR1")

sce.enhanced <- enhanceFeatures(sce.enhanced, sce, feature_names=markers,
                                     nrounds=0)