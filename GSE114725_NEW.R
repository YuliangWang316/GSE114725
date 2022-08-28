library(patchwork)
library(dplyr)
library(Seurat)
library(data.table)

# BRCA1<-fread("c:/Users/xjmik/Downloads/GSE114725_rna_imputed.csv/imputed_corrected.csv",sep = ",")
BRCA2<-fread("c:/Users/xjmik/Downloads/raw_corrected.csv",sep = ",")

BRCA_metadata<-BRCA2[,c(1:5)]
BRCA2<-BRCA2[,-1]
BRCA2<-BRCA2[,-1]
BRCA2<-BRCA2[,-1]
BRCA2<-BRCA2[,-1]
BRCA2<-BRCA2[,-1]

BRCA2<-as.data.frame(t(BRCA2))
BRCA_metadata<-as.data.frame(BRCA_metadata)
rownames(BRCA_metadata)<-BRCA_metadata$cellid
colnames(BRCA2)<-rownames(BRCA_metadata)

BRCA_metadata_Tumor<-BRCA_metadata[which(as.character(BRCA_metadata$tissue) == "TUMOR"),]
BRCA_metadata_Blood<-BRCA_metadata[which(as.character(BRCA_metadata$tissue) == "BLOOD"),]

BRCA_Tumor.data<-BRCA2[,rownames(BRCA_metadata_Tumor)]
BRCA_Blood.data<-BRCA2[,rownames(BRCA_metadata_Blood)]

BRCA_Tumor <- CreateSeuratObject(counts = BRCA_Tumor.data, project = "IMMUNE_BRCA_Tumor", min.cells = 3,meta.data = BRCA_metadata_Tumor)
BRCA_Tumor$type <- "BRCA_Tumor"
BRCA_Tumor <- subset(BRCA_Tumor, subset = nFeature_RNA > 200)
BRCA_Tumor <- NormalizeData(BRCA_Tumor, verbose = FALSE)
BRCA_Tumor <- FindVariableFeatures(BRCA_Tumor, selection.method = "vst", nfeatures = 2000)

BRCA_Blood <- CreateSeuratObject(counts = BRCA_Blood.data, project = "IMMUNE_BRCA_Blood", min.cells = 3,meta.data = BRCA_metadata_Blood)
BRCA_Blood$type <- "BRCA_Blood"
BRCA_Blood <- subset(BRCA_Blood, subset = nFeature_RNA > 200)
BRCA_Blood <- NormalizeData(BRCA_Blood, verbose = FALSE)
BRCA_Blood <- FindVariableFeatures(BRCA_Blood, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(BRCA_Blood,BRCA_Tumor), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "type")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
p1
p2
DimPlot(immune.combined, reduction = "umap", split.by = "type")

FeaturePlot(immune.combined,features = "FOXP3")
VlnPlot(immune.combined,features = "FOXP3",sort = TRUE,pt.size = 0,assay = "RNA")
Treg<-subset(immune.combined,idents = "12")
Idents(Treg)<-Treg@meta.data$tissue
VlnPlot(Treg,features = "JMJD1C",sort = TRUE,pt.size = 0)
library(ggpubr)
VlnPlot(Treg,features = "JMJD1C",sort = TRUE,pt.size = 0)+stat_compare_means()
JMJD1C<-FetchData(Treg,vars = "JMJD1C")
JMJD1C_metadata<-data.frame(rownames(Treg@meta.data),Treg@meta.data$tissue)
JMJD1C_new<-cbind(JMJD1C,JMJD1C_metadata[,2])
write.table(JMJD1C_new,file = "GSE114725.txt",sep = "\t")
