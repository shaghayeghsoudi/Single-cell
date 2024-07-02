
rm(list = ls())
#Tutorial of single-cell RNA-ATAC multiomic

## load required libraries
library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(dbplyr)
library(hdf5r)
library(AnnotationHub)


################################
######## tutorial data 1 #######
################################

## Step1. Load the data and create the Seurat object
counts_NC8 <- Read10X_h5("data/NC8/filtered_feature_bc_matrix.h5")
counts_H9 <- Read10X_h5("data/H9/filtered_feature_bc_matrix.h5")

library(AnnotationHub)
ah <- AnnotationHub()
ensdbs <- query(ah, c("EnsDb.Hsapiens"))
ensdb_id <- ensdbs$ah_id[grep(paste0(" 98 EnsDb"), ensdbs$title)]
ensdb <- ensdbs[[ensdb_id]]
seqlevelsStyle(ensdb) <- "UCSC"
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotations) <- "hg38"










################################
######## tutorial data 2 #######
################################

#Analyzing PBMC scATAC-seq
# source: https://stuartlab.org/signac/articles/pbmc_vignette



#Pre-processing workflow
#When pre-processing chromatin data, Signac uses information from two related input files, both of which can be created using CellRanger:#

root<-("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Seurat_signac_tutorial_data")

counts <- Read10X_h5(filename = paste(root,"/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5", sep = ""))
metadata <- read.csv(file = paste(root,"/atac_v1_pbmc_10k_singlecell.csv", sep = ""),header = TRUE,  row.names = 1)
                   
#creating a Seurat object

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = paste(root,"/atac_v1_pbmc_10k_fragments.tsv.gz", sep = ""),
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(counts = chrom_assay,assay = "peaks",meta.data = metadata)
pbmc
  
granges(pbmc)

### adding gene annotations to the pbmc object for the human genome
## extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg19"
# add the gene information to the object
Annotation(pbmc) <- annotations


### Computing QC Metrics

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments


#### visualization :
#The relationship between variables stored in the object metadata can be visualized using the DensityScatter() function

pdf(file = paste (root,"/atac_v1_pbmc_DensityScatter.pdf", sep = ""), height = 8, width = 14)
den_scatter<-DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
print(den_scatter)
dev.off()

## TSS plot
pdf(file = paste (root,"/atac_v1_pbmc_TSS.pdf", sep = ""), height = 8, width = 10)
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 3, 'High', 'Low')
tss<-TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
print(tss)
dev.off()

## plot fragment length
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')


##plot the distribution of each QC metric separately using a violin plot:

pdf(file = paste (root,"/atac_v1_pbmc_QCplot.pdf", sep = ""), height = 8, width = 10)
qc<-VlnPlot(
  object = pbmc,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  alpha= 0.1,
  ncol = 5
)
print(qc)
dev.off()


###remove cells that are outliers for these QC metrics.

pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 &
    #pct_reads_in_peaks > 15 &
    #blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)
pbmc

#Normalization and linear dimensional reduction


pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

DepthCor(pbmc)

# Non-linear dimension reduction and clustering

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

pdf(file = paste (root,"/atac_v1_pbmc_DimPlot.pdf", sep = ""), height = 8, width = 10)
dim<-DimPlot(object = pbmc, label = TRUE) + NoLegend()
print(dim)
dev.off()


#Create a gene activity matrix

gene.activities <- GeneActivity(pbmc)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)


## Integrating with scRNA-seq data
# Load the pre-processed scRNA-seq data for PBMCs

## To construct the object: https://github.com/satijalab/Integration2019/blob/master/preprocessing_scripts/pbmc_10k_v3.R
pbmc_rna <- readRDS(paste(root,"/pbmc_10k_v3.rds", sep = ""))
pbmc_rna <- UpdateSeuratObject(pbmc_rna)
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)


predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)


pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2


# replace each label with its most likely prediction
for(i in levels(pbmc)) {
  cells_to_reid <- WhichCells(pbmc, idents = i)
  newid <- names(which.max(table(pbmc$predicted.id[cells_to_reid])))
  Idents(pbmc, cells = cells_to_reid) <- newid
}


#Find differentially accessible peaks between cell types

DefaultAssay(pbmc) <- 'peaks'

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
head(da_peaks)


plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14+ Monocytes")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

fc <- FoldChange(pbmc, ident.1 = "CD4 Naive", ident.2 = "CD14+ Monocytes")
# order by fold change
fc <- fc[order(fc$avg_log2FC, decreasing = TRUE), ]
head(fc)

open_cd4naive <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)
head(closest_genes_cd4naive)


head(closest_genes_cd14mono)

## Plotting genomic regions

# set plotting order
levels(pbmc) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 effector","Double negative T cell","NK dim", "NK bright", "pre-B cell",'B cell progenitor',"pDC","CD14+ Monocytes",'CD16+ Monocytes')

CoveragePlot(
  object = pbmc,
  region = rownames(da_peaks)[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)


###############################################
###### Joint analysis of RNA and ATCseq #######
###############################################

# source: https://stuartlab.org/signac/articles/pbmc_multiomic

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

root<-"~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024/Seurat_signac_tutorial_data/joint_RNA_ATAC"

# load the RNA and ATAC data
counts <- Read10X_h5(paste(root,"/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5", sep = ""))
fragpath <- paste(root,"/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz", sep = "")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)


### Quality control

DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)


DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)



# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1800 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
pbmc


#Gene expression data processing
#normalize the gene expression data using SCTransform, and reduce the dimensionality using PCA
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)

#DNA accessibility data processing
DefaultAssay(pbmc) <- "ATAC"    ## feature selection
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)  ## feature selection
pbmc <- RunTFIDF(pbmc)  ## Normalization
pbmc <- RunSVD(pbmc)  ## Linear dimension reduction

## Clustering and further dimensionality reduction ##
### non-linear dimension reduction with UMAP for visualization (in the git tutorial)
p1<-ElbowPlot(pbmc,ndims = 30, reduction = "lsi")
p2<-DepthCor(pbmc, n= 30)
p2


##### Tutorial 1
RunUMAP(pbmc , reduction = "lsi", dims= 2:30, reduction.name = "umap_atac",reduction.key= "UMAPATAC_")  ## umap
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:34, k.param=13)
pbmc <- FindClusters(object = pbmc, verbose = FALSE)


p1 <- DimPlot(object = pbmc, label = TRUE, dims = c(2, 3), reduction = "lsi") +
NoLegend() +
ggtitle('SVD')
p2 <- DimPlot(object = pbmc, label = TRUE) +
NoLegend() +
ggtitle('UMAP')

p2 <- FeaturePlot(b,
c("PDGFB","CLDN5","PDGFRB","COL5A1"),
reduction = "umap_atac") & NoAxes() & NoLegend()

p1 | p2


#####
### integrated samples ###
RunUMAP(pbmc , reduction = "lsi", dims= 2:30, reduction.name = "umap_atac",reduction.key= "UMAPATAC_")  ## umap
p1 <- DimPlot(object = pbmc, label = TRUE, dims = c(2, 3), reduction = "lsi") +


