# Download public available data ( DOI：10.1101/gr.267013.120 )

library(Seurat)
library(patchwork)
library(knitr)
library(rmarkdown)
library(gt)
library(ggrepel)
library(org.Mm.eg.db)
library(clusterProfiler)
library(pathview)
library(AnnotationDbi)
library(GO.db)
library(DOSE)
library(tidyverse)
library(ComplexHeatmap)
library(future)
library(CellChat)
library(msigdbr)
library(ggsignif)
library(ggridges)

# make seurat object
# 0h
expression_matrix <- ReadMtx(mtx = "./GSM4572241_Adult.matrix.mtx.gz",cells = "./GSM4572241_Adult.barcodes.tsv.gz",features = "./GSM4572241_Adult.features.tsv.gz")
hep_0h <- CreateSeuratObject(counts = expression_matrix,project = "0h", min.cells = 3, min.features = 200)

# 24h
expression_matrix <- ReadMtx(mtx = "./GSM4572242_24.matrix.mtx.gz",cells = "./GSM4572242_24.barcodes.tsv.gz",features = "./GSM4572242_24.features.tsv.gz")
hep_24h <- CreateSeuratObject(counts = expression_matrix,project = "24h", min.cells = 3, min.features = 200)

# 48h
expression_matrix <- ReadMtx(mtx = "./GSM4572243_48.matrix.mtx.gz",cells = "./GSM4572243_48.barcodes.tsv.gz",features = "./GSM4572243_48.features.tsv.gz")
hep_48h <- CreateSeuratObject(counts = expression_matrix,project = "48h", min.cells = 3, min.features = 200)

# 96h
expression_matrix <- ReadMtx(mtx = "./GSM4572244_96.matrix.mtx.gz",cells = "./GSM4572244_96.barcodes.tsv.gz",features = "./GSM4572244_96.features.tsv.gz")
hep_96h <- CreateSeuratObject(counts = expression_matrix,project = "96h", min.cells = 3, min.features = 200)

saveRDS(hep_0h,"./make_obj/hep_0h_ver1")
saveRDS(hep_24h,"./make_obj/hep_24h_ver1")
saveRDS(hep_48h,"./make_obj/hep_48h_ver1")
saveRDS(hep_96h,"./make_obj/hep_96h_ver1")

remove(list = ls())

# Memory can be safely released at this point.


# merge object
hep_0h <- readRDS("./make_obj/hep_0h_ver1")
hep_24h <- readRDS("./make_obj/hep_24h_ver1")
hep_48h <- readRDS("./make_obj/hep_48h_ver1")
hep_96h <- readRDS("./make_obj/hep_96h_ver1")
hep_0h$percent_mt <- PercentageFeatureSet(hep_0h, pattern = "^mt-")
hep_24h$percent_mt <- PercentageFeatureSet(hep_24h, pattern = "^mt-")
hep_48h$percent_mt <- PercentageFeatureSet(hep_48h, pattern = "^mt-")
hep_96h$percent_mt <- PercentageFeatureSet(hep_96h, pattern = "^mt-")

hep_0h <- subset(hep_0h,subset = nCount_RNA < 2500 & nFeature_RNA < 2000 & percent_mt<10)
hep_24h <- subset(hep_24h,subset = nCount_RNA < 2500 & nFeature_RNA < 2000 & percent_mt<10)
hep_48h <- subset(hep_48h,subset = nCount_RNA < 2500 & nFeature_RNA < 2000 & percent_mt<10)
hep_96h <- subset(hep_96h,subset = nCount_RNA < 2500 & nFeature_RNA < 2000 & percent_mt<10)

saveRDS(hep_0h,"./make_obj/hep_0h_ver2")
saveRDS(hep_24h,"./make_obj/hep_24h_ver2")
saveRDS(hep_48h,"./make_obj/hep_48h_ver2")
saveRDS(hep_96h,"./make_obj/hep_96h_ver2")

# Memory can be safely released at this point.
remove(list = ls())
hep_0h <- readRDS("./make_obj/hep_0h_ver2")
hep_24h <- readRDS("./make_obj/hep_24h_ver2")
hep_48h <- readRDS("./make_obj/hep_48h_ver2")
hep_96h <- readRDS("./make_obj/hep_96h_ver2")

merge.list <- list(hep_0h,hep_24h,hep_48h,hep_96h)
merge.list <- lapply(X = merge.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = merge.list)
merge.list <- lapply(X = merge.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

obj.anchors <- FindIntegrationAnchors(object.list = merge.list, anchor.features = features, reduction = "rpca")

hep <- IntegrateData(anchorset = obj.anchors,normalization.method = "LogNormalize")
saveRDS(hep,"./make_obj/hepatectomy_merge_ver1")

# Memory can be safely released at this point.

# QC
remove(list = ls())
hep <- readRDS("./make_obj/hepatectomy_merge_ver1")
DefaultAssay(hep) <- "RNA"
hep <- subset(hep,subset = nCount_RNA < 2000 & nFeature_RNA < 1000 & percent_mt<7.5)


# Create three objects, each containing 5,000 randomly sampled cells from each time point.
cell_0h <- rownames(hep@meta.data[hep@meta.data$orig.ident=="0h",])
cell_24h <- rownames(hep@meta.data[hep@meta.data$orig.ident=="24h",])
cell_48h <- rownames(hep@meta.data[hep@meta.data$orig.ident=="48h",])
cell_96h <- rownames(hep@meta.data[hep@meta.data$orig.ident=="96h",])
cell_list <- list(cell_0h,cell_24h,cell_48h,cell_96h)
names(cell_list) <- c("cell_0h","cell_24h","cell_48h","cell_96h")

cell_sample <- function(x){
  sample_0h <- sample(x[[1]],size = 5000)
  sample_24h <- sample(x[[2]],size = 5000)
  sample_48h <- sample(x[[3]],size = 5000)
  sample_96h <- sample(x[[4]],size = 5000)
  reduction_cell <- c(sample_0h,sample_24h,sample_48h,sample_96h)
  return(reduction_cell)
}

for (i in 1:3) {
  set.seed(i*10)
  cell <- cell_sample(cell_list)
  a <- subset(hep,cells = cell)
  assign(paste0("hep_redu_",i),a)
}

saveRDS(hep_redu_1,"./make_obj/hepatectomy_redu1_ver1")
saveRDS(hep_redu_2,"./make_obj/hepatectomy_redu2_ver1")
saveRDS(hep_redu_3,"./make_obj/hepatectomy_redu3_ver1")


# Clustering
remove(list = ls())

 # redu1 object
hep1 <- readRDS("./make_obj/hepatectomy_redu1_ver1")
hep1 <- FindVariableFeatures(hep1, selection.method = "vst", nfeatures = 2000)
gene_all <- rownames(hep1)
hep1 <- ScaleData(hep1, features = gene_all)
hep1 <- RunPCA(hep1, features = VariableFeatures(object = hep1))
hep1 <- FindNeighbors(object = hep1,dims = 1:15)
hep1 <- FindClusters(object = hep1,resolution = 1.2)
hep1 <- RunUMAP(hep1,dims = 1:15)

 # redu2 object
hep2 <- readRDS("./make_obj/hepatectomy_redu2_ver1")
hep2 <- FindVariableFeatures(hep2, selection.method = "vst", nfeatures = 2000)
gene_all <- rownames(hep2)
hep2 <- ScaleData(hep2, features = gene_all)
hep2 <- RunPCA(hep2, features = VariableFeatures(object = hep2))
hep2 <- FindNeighbors(object = hep2,dims = 1:15)
hep2 <- FindClusters(object = hep2,resolution = 1.2)
hep2 <- RunUMAP(hep2,dims = 1:15)

 # redu3 object
hep3 <- readRDS("./make_obj/hepatectomy_redu3_ver1")
hep3 <- FindVariableFeatures(hep3, selection.method = "vst", nfeatures = 2000)
gene_all <- rownames(hep3)
hep3 <- ScaleData(hep3, features = gene_all)
hep3 <- RunPCA(hep3, features = VariableFeatures(object = hep3))
hep3 <- FindNeighbors(object = hep3,dims = 1:15)
hep3 <- FindClusters(object = hep3,resolution = 1.2)
hep3 <- RunUMAP(hep3,dims = 1:15)

saveRDS(hep1,"./make_obj/hepatectomy_redu1_ver2")
saveRDS(hep2,"./make_obj/hepatectomy_redu2_ver2")
saveRDS(hep3,"./make_obj/hepatectomy_redu3_ver2")