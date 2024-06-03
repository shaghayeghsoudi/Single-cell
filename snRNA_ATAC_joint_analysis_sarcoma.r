
rm(list = ls())
### Load required libraries
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(ggplot2)
library(patchwork)



#Pre-processing workflow
#When pre-processing chromatin data, Signac uses information from two related input files, both of which can be created using CellRanger:#

root<-("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024")

counts <- Read10X_h5(filename = paste(root,"/med399-1/filtered_feature_bc_matrix.h5", sep = ""))
fragpath <- paste(root,"/med399-1/atac_fragments.tsv.gz", sep = "")
#metadata <- read.csv(file = paste(root,"/med399-1/atac_v1_pbmc_10k_singlecell.csv", sep = ""),header = TRUE,  row.names = 1)
                   
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

# pbmc
#An object of class Seurat 
#78268 features across 1977 samples within 2 assays 
#Active assay: RNA (36601 features, 0 variable features)
# 1 layer present: counts
# 1 other assay present: ATAC


### Quality control
DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

# generate density plot
pdf(file = paste(root,"scatter_quality_plot_SRC3991-1.pdf", sep = ""), width = 10, height= 8)
denscatter<-DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
print(denscatter)
dev.off()


# box plots
pdf(file = paste(root,"violin_quality_plot_SRC3991-1.pdf", sep = ""), width = 10, height= 5)
violin_qual<-VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
print(violin_qual)
dev.off()


# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 50000 &
    nCount_ATAC > 1000 &
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

p2 <- FeaturePlot(pbmc,
c("PDGFB","CLDN5","PDGFRB","CD99"),
reduction.name = "umap_atac") & NoAxes() & NoLegend()

#p1 | p2


#####
### integrated samples ###
RunUMAP(pbmc , reduction = "lsi", dims= 2:30, reduction.name = "umap_atac",reduction.key= "UMAPATAC_")  ## umap
p1 <- DimPlot(object = pbmc, label = TRUE, dims = c(2, 3), reduction = "lsi") +


