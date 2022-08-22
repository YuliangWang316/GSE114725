library(patchwork)
library(dplyr)
library(Seurat)
library(data.table)

# BRCA1<-fread("c:/Users/xjmik/Downloads/GSE114725_rna_imputed.csv/imputed_corrected.csv",sep = ",")
BRCA2<-fread("c:/Users/xjmik/Downloads/GSE114725_rna_raw.csv/raw_corrected.csv",sep = ",")

BRCA_metadata<-BRCA2[,c(1:5)]
BRCA2<-BRCA2[,-1]
BRCA2<-BRCA2[,-1]
BRCA2<-BRCA2[,-1]
BRCA2<-BRCA2[,-1]
BRCA2<-BRCA2[,-1]

BRCA2<-as.data.frame(t(BRCA2))
colnames(BRCA2)<-rownames(BRCA_metadata)
pbmc <- CreateSeuratObject(counts = BRCA2, project = "pbmc3k", min.cells = 3, min.features = 200,meta.data = BRCA_metadata)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 50)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:20)
DimPlot(pbmc)
FeaturePlot(pbmc,features = "FOXP3")
VlnPlot(pbmc,features = "FOXP3",sort = TRUE,pt.size = 0)
Treg<-subset(pbmc,idents = "13")
Idents(Treg)<-Treg@meta.data$tissue
VlnPlot(Treg,features = "JMJD1C",sort = TRUE,pt.size = 0)
Treg_T_P<-subset(Treg,idents = c("TUMOR","BLOOD"))
VlnPlot(Treg_T_P,features = "JMJD1C",sort = TRUE,pt.size = 0)
library(ggpubr)
VlnPlot(Treg_T_P,features = "JMJD1C",sort = TRUE,pt.size = 0)+stat_compare_means()
