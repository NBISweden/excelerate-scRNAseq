#'---
#'title: "Trajectory inference"
#'author: "Jules GILET"
#'
#'---

#' ## Datasets

# Transcriptional trajectories will be inferred from data by Nestorowa, Hamey et al. (Blood, 2016):
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5305050/
# The dataset consists in 1600 hematopoietic stem and progenitor cells from mouse bone marrow
# (sequenced using the SMARTseq2 technology). Using flow cytometry and index sorting, 12 HSPC of different phenotypes
# (about 10 cells each) have been included in the dataset, and will be used in this lab as a biological prior for the
# identification of the root and the branches in the transcriptional trajectory models.

# bash:
# 
# mkdir data && cd data
# wget -nH -nd    http://blood.stemcells.cam.ac.uk/data/nestorowa_corrected_log2_transformed_counts.txt.gz \
# http://blood.stemcells.cam.ac.uk/data/nestorowa_corrected_population_annotation.txt.gz
# wget -nH -nd ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81682/suppl/GSE81682%5FHTSeq%5Fcounts%2Etxt%2Egz
# 
# gunzip *
# cd ..

#' ======

#' ## Part I - Monocle2/DDRtree
# Inference done with Monocle2/DDRtree available via Bioconductor

library(monocle)
library(biomaRt)

# data loading

# The authors provide an expression matrix that has been filtered (highly expressed genes, high quality cells)
# scaled and log-normalized. An annotation table is also provided, with each cell type labelled according to
# the immunophenotyping done by flow cytometry.

lognorm <- t(read.table('data/nestorowa_corrected_log2_transformed_counts.txt', sep=" ", header=TRUE))
annotable <- read.table('data/nestorowa_corrected_population_annotation.txt')

# To infer a trajectory with Monocle2/DDRtree, using non-normalized UMI-based counts is highly recommended,
# as Monocle2 will scale and normalize the data internaly and is especting data ditributed according to a negative binomial.

# The count matrix has been downloaded and will be used for Monocle2.

counts <- read.table('data/GSE81682_HTSeq_counts.txt', sep="\t", header=TRUE, row.names='ID')

counts[1:5,1:5]
dim(counts)

lognorm[1:5,1:5]
dim(lognorm)

# Note that the count matrix is not filtered, and genes are labelled
# according to ensembl gene IDs. We will first filter the matrix according to the authors choices
# and we will map the gene official symbols.

# we filter the counts to keep only high quality cells

counts <- counts[ , colnames(lognorm) ]
dim(counts)

# we create a annotation data frame to label the cell types
# as defined byu the authors
pDat <- data.frame(cell=colnames(counts), celltype='undefined', stringsAsFactors=FALSE)
rownames(pDat) <- pDat$cell
pDat[ rownames(anno), 2] <- as.character(anno$celltype)
head(pDat)

# we create a feature annotation data frame
# that will contain gene informations and matching symbols and ID
# we now annotate the genes IDs in the counts matrix with biomaRt
mart <- biomaRt::useDataset("mmusculus_gene_ensembl", biomaRt::useMart("ensembl"))
genes_table <- biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name"), values=rownames(counts), mart=mart)
rownames(genes_table) <- genes_table$ensembl_gene_id
head(genes_table)
fDat <- genes_table[ rownames(counts), ]
# to be consistent with Monocle naming conventions
colnames(fDat) <- c('ensembl_gene_id', 'gene_short_name')
head(fDat)

# we can now use this table to filter the genes
# in the counts matrix that are highly expressed
# according to the quality filters used by the authors

fDat <- fDat[ fDat$gene_short_name %in% rownames(lognorm), ]

# and we finally keep in the counts matrix only these genes

counts <- counts[ rownames(fDat), ]
dim(counts)
dim(fDat)
dim(pDat)

# we build a cell dataset object
# in an appropriate format for monocle
# default method for modeling the expression values is VGAM::negbinomial.size()
# and adapted to counts

cds <- newCellDataSet(as.matrix(counts), phenoData=Biobase::AnnotatedDataFrame(pDat), featureData=Biobase::AnnotatedDataFrame(fDat))
cds

# the monocle cds object is done
# and ready for trajectory inference
dir.create('monocle', showWarnings=FALSE)
saveRDS(cds, 'monocle/cds_hematopoiesis.rds')

# Monocle2 preprocess
# normalization and scaling
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# no need for further filtering, the expression matrix has already be filtered
# we find the genes that are expressed
cds <- detectGenes(cds, min_expr=0.1)
print(head(fData(cds)))

# we identify genes that are expressed in at least 10 cells
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
length(expressed_genes)
print(head(pData(cds)))

# Identification of the ordering genes by differential testing
# ie. genes that are presumed to be important in the differentiation
# process captured in the sample. We used the cell types identified by the authors
# to define the ordering genes by DE testing. Alternatively, a classical approach
# consist in clustering the cells, then identify markers genes per clusters.

diff_test_res <- differentialGeneTest(cds[ expressed_genes, ], fullModelFormulaStr="~ celltype")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
length(ordering_genes)

# we mark the genes that will be used for the ordering 
cds <- setOrderingFilter(cds, ordering_genes)

# We use the DDRTree algorithm to infer a trajectory with potential branching
# points.
cds <- reduceDimension(cds, max_components = 2, method='DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by="celltype")

# changing the cell color
cell_colors <- c('lightblue','blue','red','black','orange','yellow','turquoise','lightgrey')
plot_cell_trajectory(cds, color_by="celltype") + scale_color_manual(values=cell_colors)

# The most imamture HSC in the sample express E-Slam
# we will define the root of this model according to this subset of cells
table(pData(cds)$State, pData(cds)$celltype)[,"ESLAM"]

# The state 4 defines the root
# We define the root in the model
cds <- orderCells(cds, root_state = 4)
# The pseudotime is now defined by the distance to the root
plot_cell_trajectory(cds, color_by = "Pseudotime")

# Differential expression testing per branch
# We look at the genes that are differentially expressed
# according to the pseudotime model
diff_test_res <- differentialGeneTest(cds[ ordering_genes, ], fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(cds[ sig_gene_names[1:50], ], num_clusters = 3, cores=4, show_rownames=TRUE)

# Differential expression per branch 
# let's look for genes involved in the erythroid pathway
BEAM_res <- BEAM(cds, branch_point = 1, cores = 4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)

plot_genes_branched_heatmap(cds[row.names(BEAM_res)[1:50]], branch_point = 1, num_clusters = 3, cores=4, use_gene_short_name=TRUE, show_rownames=TRUE)

# There is a clear separation between genes that are involved in the erythroid
# differentiation (eg. Gata1) on the left (cell fate1)
# with genes involved in the leukocyte differentiation (eg. Sell, Ccl9)

plot_genes_branched_pseudotime(cds[row.names(BEAM_res)[1:5]], branch_point = 1, color_by = "celltype", ncol = 1)  + scale_color_manual(values=cell_colors)


#
#' ======

#' ## Part II - Diffusion map

# Analysis and inference done with the destiny package available via Bioconductor

# Trajectory inference by diffusion map an diffusion pseudotime

library(destiny)
library(ggplot2)

# data loading
# we now will directly use the filtered, scaled, log-normalised expression matrix
# provided by the author of the article
lognorm <- t(read.table('data/nestorowa_corrected_log2_transformed_counts.txt', sep=" ", header=TRUE))
lognorm[1:5,1:5]

# We load the annotation of cell types that has been defined
# using flow cytometry and index sorting.
# The cell subsets (final differentiation stages) will be used
# to validate the trajectory model
anno <- read.table('data/nestorowa_corrected_population_annotation.txt')
pDat <- data.frame(cell=colnames(lognorm), celltype='undefined', stringsAsFactors=FALSE)
rownames(pDat) <- pDat$cell
pDat[ rownames(anno), 2] <- as.character(anno$celltype)

# we build an expression set object for an easier integration with destiny
eset <- Biobase::ExpressionSet(lognorm, phenoData=Biobase::AnnotatedDataFrame(pDat))
eset

# the expression set is ready for inference with destiny
dir.create('destiny', showWarnings=FALSE)
saveRDS(cds, 'destiny/eset_hematopoiesis.rds')

# diffusion map
dmap <- DiffusionMap(eset)
# the process takes less than 60 seconds

# We look at the global model

plot.DiffusionMap(dmap)
plot.DiffusionMap(dmap, dims=c(1,2))
plot.DiffusionMap(dmap, dims=c(2,3))
plot.DiffusionMap(dmap, dims=c(1,3))

# components 1-2 describe well the branching process between erythroid (red) and myeloid/lymphoid (white) lineages

# we use ggplot2 to have a better rendering and project the cell labels
# as defined by flow cytometry experiment and index sorting

qplot(DC1, DC2, data=dmap, colour=celltype) + scale_color_manual(values=c('lightblue','brown','red','black','orange','yellow','blue','lightgrey')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Pseudotime inference by diffusion:
# the transcriptional distance between cells is estimated by
# random walk along a neighborhood graph.
# The resulting "transcriptional" transition probability
# between cells is used to infer a pseudo-time scale of the differentiation process.

# we first define a root cell (origin) for the model
# we find the index of a ESLAM positive cells
which(anno$celltype=="ESLAM")

# we use this cell as a starting point

dpt <- DPT(dmap, tips=19)
plot(dpt)

# We can project the level of expression of known marker genes on
# the trajectory model

# Procr / Endothelial protein C is a marker of HSC subsets
plot(dpt, col_by='Procr', pal=viridis::magma)

# Gata1 is a key TF of the erythroid lineage
plot(dpt, col_by='Gata1', pal=viridis::magma)

# cathepsin G is a marker of neutrophils
plot(dpt, col_by='Ctsg', pal=viridis::magma)

sessionInfo()
