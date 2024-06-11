
rm(list = ls())
### Load required libraries
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(ggplot2)
library(patchwork)
library(simspec)  ## for data integration
library("harmony")



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
c("PDGFB","CLDN5","PDGFRB","CD99"),alpha = 0.2)
reduction.name = "umap_atac") & NoAxes() & NoLegend()

#p1 | p2

### save data here
save(pbmc, file="pbmc.test.Rda")

#####
### integrated samples ###
RunUMAP(pbmc , reduction = "lsi", dims= 2:30, reduction.name = "umap_atac",reduction.key= "UMAPATAC_")  ## umap
p1 <- DimPlot(object = pbmc, label = TRUE, dims = c(2, 3), reduction = "lsi") +


##########################
##### joint samples ######
##########################
root<-("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024")


counts_399_1 <- Read10X_h5(filename = paste(root,"/med399-1/filtered_feature_bc_matrix.h5", sep = ""))
counts_399_2 <- Read10X_h5(filename = paste(root,"/med399-2/filtered_feature_bc_matrix.h5", sep = ""))


ATAC_frag399_1<-"~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024/med399-1/atac_fragments.tsv.gz"
ATAC_frag399_2<-"~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024/med399-2/atac_fragments.tsv.gz"



#### get gene annotations from hg38 assembly 
ah <- AnnotationHub()
ensdbs <- query(ah, c("EnsDb.Hsapiens"))
ensdb_id <- ensdbs$ah_id[grep(paste0(" 98 EnsDb"), ensdbs$title)]
ensdb <- ensdbs[[ensdb_id]]
seqlevelsStyle(ensdb) <- "UCSC"
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotations) <- "hg38"


# create a Seurat object containing the RNA adata
seurat_SRC399_1 <- CreateSeuratObject(counts = counts_399_1$`Gene Expression`,
                                 assay = "RNA",
                                 project = "SRC399-1")



seurat_SRC399_2 <- CreateSeuratObject(counts = counts_399_2$`Gene Expression`,
                                      assay = "RNA",
                                      project = "SRC399-2")


seurat_SRC399_1[['ATAC']] <- CreateChromatinAssay(counts = counts_399_1$`Peaks`,
                                             annotation = annotations,fragments =ATAC_frag399_1,
                                             sep = c(":", "-"),
                                             genome = 'hg38')
                                             
                                               
                                             

seurat_SRC399_2[['ATAC']] <- CreateChromatinAssay(counts = counts_399_2$`Peaks`,
                                                  annotation = annotations,fragments =ATAC_frag399_2,
                                                  sep = c(":", "-"),
                                                  genome = 'hg38')


### merge for samples (works for RNA)
seurat_src <- merge(seurat_SRC399_1, seurat_SRC399_2)      

## merge ATACseq assays 


peaks <- reduce(unlist(as(c(seurat_SRC399_1@assays$ATAC@ranges,
                            seurat_SRC399_2@assays$ATAC@ranges),
                          "GRangesList")))
peakwidths <- width(peaks)
peaks <- peaks[peakwidths < 20000 & peakwidths > 7]


counts_atac_merged <- FeatureMatrix(seurat_src@assays$ATAC@fragments,
                                    features = peaks,
                                    cells = colnames(seurat_src))

seurat_src[['ATAC']] <- CreateChromatinAssay(counts_atac_merged,
                                         fragments = seurat_src@assays$ATAC@fragments,
                                         annotation =
                                           seurat_src@assays$ATAC@annotation,
                                         sep = c(":",
                                                 "-"),
                                         genome = "hg38")



###############
### Step 2 ####
###############
### Quality control


seurat_src <- PercentageFeatureSet(seurat_src, pattern = "^MT-", col.name = "percent.mt",
                               assay = "RNA")
seurat_src <- NucleosomeSignal(seurat_src, assay = "ATAC")
seurat_src <- TSSEnrichment(seurat_src, assay = "ATAC")

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024/VlnPlot_quality_combined_SRC_samples.pdf",height = 6, width =12 )
quality_plot<-VlnPlot(seurat_src,
        features = c("nFeature_RNA",
                     "percent.mt",
                     "nFeature_ATAC",
                     "TSS.enrichment",
                     "nucleosome_signal"),
                      ncol = 5,
                      pt.size = 0)
print(quality_plot)
dev.off()

#### filter our low quality cells 

seurat_src <- subset(seurat_src,
                 subset = nFeature_RNA > 100 &   ### adjust thresholds
                   nFeature_RNA < 10000 &
                   percent.mt < 30 &
                   nFeature_ATAC > 100 &
                   nFeature_ATAC < 20000 &
                   TSS.enrichment > 1 &
                   nucleosome_signal < 3
)


#################################
##### Step 3 ####################
## Analysis on the RNA assay ####
#################################

DefaultAssay(seurat_src) <- "RNA"

seurat_src <- NormalizeData(seurat_src) %>%   ### Normalizing
  FindVariableFeatures(nfeatures = 3000) %>%  ### Finding top vaiables
  #CellCycleScoring(s.features = cc.genes.updated.2019$s.genes,
                   #g2m.features = cc.genes.updated.2019$g2m.genes) %>%
  ScaleData() %>%         ### Scaling data
  RunPCA(npcs = 50) %>%   ### Running PCA
  RunUMAP(dims = 1:20, reduction.name = "umap_rna", reduction.key = "UMAPRNA_")    ###UMAP embedding

## checking batch effect
p1 <- DimPlot(seurat_src, group.by = "orig.ident", reduction = "umap_rna",pt.size =1.5) 
p1


p2 <- FeaturePlot(seurat_src,
c("IGFBP7", "PREX2", "CALCRL", "ADAMTS9", "CYYR1", "SPP1", "VEGFC"," NFIB"), alpha = 1, pt.size =1)
#reduction = "umap_css_rna") & NoAxes() & NoLegend()
pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024/RNA_feature_plot_combined_SRC_samples_batch_effect_relaxed.pdf",height = 14, width =19)
both<-p1 + p2
print(both)
dev.off()

#########################################
#### Visualzie top variable feature #####
#### OPTIONAL ####
pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024/LabelPoint_plot_snRNA_top_features.pdf",height = 9, width =14)
top_features <- head(VariableFeatures(seurat_src), 20)
plot1 <- VariableFeaturePlot(seurat_src)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
top_vars<-plot1 + plot2
print(top_vars)
dev.off()


###########################################################
### performing data integration to fix batch effect(RNA) ###

#seurat_src <- cluster_sim_spectrum(seurat_src,
#label_tag = "orig.ident",
#cluster_resolution = 0.6,
#reduction.name = "css_rna",
#reduction.key = "CSSRNA_")

### Using Harmony 
seurat_src <- RunHarmony(seurat_src, group.by.vars = "orig.ident",plot_converge=TRUE,max_iter = 10, reduction = "pca",reduction.save = "harmony_rna",assay.use="RNA")


## NOTE: To make sure our Harmony integration is reflected in the data visualization, we still need to generate a UMAP derived from these harmony embeddings instead of PCs:
harmonized_seurat <- RunUMAP(seurat_src, reduction = "harmony_rna", assay = "RNA", dims = 1:40)


seurat_src <- RunUMAP(seurat_src,
    reduction = "harmony",
    dims = 1:ncol(Embeddings(seurat_src,"harmony")))


### save the object
p1 <- DimPlot(seurat_src, group.by = "orig.ident", reduction = "umap",pt.size =1.5) 

p2 <- FeaturePlot(seurat_src,
c("PTPRC","COL5A1","CD4","CD68","CD99","CTSK","COL1A1","FBLN1","TOP2A"),
reduction = "umap",alpha = 1, pt.size =1)


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024/RNA_feature_plot_combined_SRC_samples_Harmony.pdf",height = 18, width =24)
both<-p1 + p2
print(both)
dev.off()

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024/P1_RNA_feature_plot_combined_SRC_samples_Harmony.pdf",height = 18, width =24)
p2 <- FeaturePlot(seurat_src,
c("PTPRC","COL5A1","CD4","CD68","CD99","CTSK","COL1A1","FBLN1","TOP2A"),
reduction = "umap",alpha = 1, pt.size =1)
print(p2)
dev.off()

#p2 <- FeaturePlot(seurat_src,
#c("PREX2","CD99", "CAV1","KDSR","CALCRL", "CYYR1", "SNTG2", "SPP1", "RAMP3", "ADAMTS9", "IGFBP5", "TPO" ),
#reduction = "umap",alpha = 1, pt.size =1)
#print(p2)


## clustering and cluster marker identification

seurat_src <- FindNeighbors(seurat_src,
reduction = "harmony",
dims = 1:ncol(Embeddings(seurat_src,"harmony"))) %>%
FindClusters(resolution = 0.2)


#DE_cl_rna <- presto::wilcoxauc(seurat_src, "RNA_snn_res.0.2")
#top_markers <- DE_cl_rna %>%
#filter(logFC > log(1.2) &
#auc > 0.7 &
#padj < 0.05 &
#pct_in - pct_out > 30 &
#pct_out < 30) %>%
#group_by(group) %>%
#top_n(1, wt = auc)


p1 <- DimPlot(seurat_src,
group.by="RNA_snn_res.0.2",
reduction="umap", label=T,pt.size =1.9, alpha = 0.7) 

p2 <- FeaturePlot(seurat_src,
c("PTPRC","COL5A1","CD4","CD68","CD99","CTSK","COL1A1","FBLN1","TOP2A"),
order = T,
ncol=3) 

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024/clustered_RNA_feature_plot_combined_SRC_samples_Harmony.pdf",height = 18, width =24)
both_clustered<-(p1 | p2) + patchwork::plot_layout(widths = c(2,3))
print(both_clustered)
dev.off()


# You may also want to save the object
#saveRDS(seurat, file="integrated_seurat.rds")

##########################################
### Step 4. Analysis on the ATAC assay ###
##########################################

#Feature selection
DefaultAssay(seurat_src) <- "ATAC"
seurat_src <- FindTopFeatures(seurat_src, min.cutoff = 50)

## Normalization
seurat_src <- RunTFIDF(seurat_src)

## linear dimenssion reduction ##
seurat_src <- RunSVD(seurat_src, n = 50)


## Non-linear dimension reduction with UMAP for visualization ##
p1 <- ElbowPlot(seurat_src, ndims = 30, reduction="lsi")
p2 <- DepthCor(seurat_src, n = 30)
p1 + p2


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024/ATAC_ellbow_depthCor.pdf",height = 6, width =7)
atac_both_raw<- p1+ p2
print(atac_both_raw)
dev.off()
#2nd to the 30th SVD components to generate the UMAP embedding of the ATAC assay

seurat_src <- RunUMAP(seurat_src,
reduction = "lsi",
dims = 2:30,
reduction.name = "umap_atac",
reduction.key = "UMAPATAC_")

p1 <- DimPlot(seurat_src,
group.by = "orig.ident",
reduction = "umap_atac",pt.size =1.5) 

p2 <- FeaturePlot(seurat_src,
c("PTPRC","COL5A1","CD4","CD68","CD99","CTSK","COL1A1","FBLN1","TOP2A"),
reduction = "umap_atac")
#reduction = "umap_atac") 


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024/ATAC_feature_plot_combined_SRC_samples.pdf",height = 14, width =19)
both_ATAC<-p1 | p2
print(both_ATAC)
dev.off()


### Data integration of the ATAC assay using Harmony
seurat_src <- RunHarmony(seurat_src,
    group.by.vars = "orig.ident",
    reduction = "lsi",   ### it should be SVD embedding ("Isi") instead of PCA embedding
    dims.use = 2:30,
    max.iter.harmony = 50,
    reduction.save = "harmony_atac")


seurat_src <- RunUMAP(seurat_src,
    reduction = "harmony",
    dims = 2:30,
    reduction.name = "umap_atac",
    reduction.key = "umapatac_")


#hm.integrated <- RunHarmony(seurat_src, group.by.vars ="orig.ident", reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)
#hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')
#DimPlot(hm.integrated, group.by = 'dataset', pt.size = 0.1)

p1<-DimPlot(seurat_src, group.by = "orig.ident",reduction = "umap_atac",pt.size =1.5)



p2 <- FeaturePlot(seurat_src,
    c("PTPRC","COL5A1","CD4","CD68","CD99","CTSK","COL1A1","FBLN1","TOP2A"),
    reduction = "umap_atac",pt.size =1.2)

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/Medgenome_multiome10X_March2024/ATAC_Integrated_harmony_feature_plot_combined_SRC_samples.pdf",height = 14, width =19)
both_ATAC_hm<-p1 + p2
print(both_ATAC_hm)
dev.off()

#################################################################################
### Section 2. Bi-modal integrative analysis of the RNA-ATAC scMultiome data ###
#################################################################################

#Step 1. Weighted nearest neighbor analysis
seurat_src <- FindMultiModalNeighbors(seurat_src,
        reduction.list = list("harmony", "harmony"),
        dims.list =
        list(1:ncol(Embeddings(seurat_src,"harmony")),
        1:ncol(Embeddings(seurat_src,"harmony"))),
        modality.weight.name =
        c("RNA.weight","ATAC.weight"),
        verbose = TRUE)