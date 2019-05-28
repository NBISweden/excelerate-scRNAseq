#'---
#'title: "Trajectory inference on random data"
#'author: "Jules GILET"
#'output:
#'  html_document:
#'    keep_md: true
#'---

# a synthetic expression matrix is generated (see the introduction lecture on day 1)
# the extraDistr package is required to model zero-inflated negative binomial
# distributions
# we simulate a sparse matrix of 1000 cells with 25000 genes

#' ## Data generation and loading

# emat <- Matrix::Matrix(data=extraDistr::rzinb(25000*1000, 50, 0.95, 0.75), nrow=25000, ncol=1000, sparse=TRUE)
# saveRDS(emat, file='random_sparse_emat.rds')

# for simplicity a rds file of the random matrix is provided

emat <- readRDS('random_sparse_emat.rds')
dim(emat)
emat[1:5,1:5]

#' ## Clustering of the articifial dataset

# we use Seurat to partition data into clusters
library(Seurat)

seu <- CreateSeuratObject(as.matrix(emat), min.features=500, min.cells=5)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- FindNeighbors(seu)
seu <- FindClusters(seu)
seu <- RunUMAP(seu, dims=1:20)
DimPlot(seu, reduction='umap')
head(seu@meta.data)
meta <- seu@meta.data

#' ## Monocle/DDRTree model

# creation of a cell dataset object for monocle TI
library(monocle)

pDat <- data.frame(cell=colnames(emat), celltype=meta$seurat_clusters, stringsAsFactors=FALSE)
fDat <- data.frame(gene_short_name=rownames(emat), stringsAsFactors=FALSE)
head(fDat)
head(pDat)
rownames(fDat) <- rownames(emat)
rownames(pDat) <- colnames(emat)
library(monocle)
cds <- newCellDataSet(as.matrix(emat), phenoData=Biobase::AnnotatedDataFrame(pDat), featureData=Biobase::AnnotatedDataFrame(fDat))
cds

# notmalisation and detection of 'genes' that are expressed in many 'cells'
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr=0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 250))
length(expressed_genes)

# differential expression testing: it is possible to find some DE 'genes' from random data
diff_test_res <- differentialGeneTest(cds[ expressed_genes, ], fullModelFormulaStr="~ celltype")
ordering_genes <- row.names(subset(diff_test_res, pval < 0.05))
length(ordering_genes)

# inference of the DDRTree model according to these ordering genes
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method='DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by="celltype")

# the model articially create a couple of branches with group of cell clusters
# on the tips

#' ## Diffusion map model

# as a comparison we will use diffusion map on random data
library(scran)
library(destiny)

sce <- SingleCellExperiment(assays=list(counts=emat, logcounts=log1p(emat)))
sce <- normalize(sce)

dm <- DiffusionMap(sce)
plot(dm)

# There is no clear structure in the data, even in the firsts eigen vectors (diffusion components)

sessionInfo()
