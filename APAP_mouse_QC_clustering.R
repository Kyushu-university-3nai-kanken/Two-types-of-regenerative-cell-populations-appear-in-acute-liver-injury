# Download public available data from https://doi.org/10.1016/j.stem.2022.04.008

library(tidyverse)
library(Seurat)
library(patchwork)
library(readr)
library(knitr)
library(rmarkdown)
library(gt)
library(hdf5r)
library(ggrepel)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationDbi)

# 48h 72h sample
rawdata <- Read10X_h5(filename = "hep_48h_72h_raw_feature_bc_matrix.h5")
rawdata.48.72 <- rawdata$`Gene Expression`
rawsample <- rawdata$`Antibody Capture`
hep48_72 <- CreateSeuratObject(counts = rawdata.48.72,min.cells = 3,min.features = 100)
hep48_72 <- NormalizeData(object = hep48_72,normalization.method ="LogNormalize",scale.factor = 10000)
joint.bcd <- intersect(rownames(hep48_72@meta.data),colnames(rawsample))
rawsample <- rawsample[,joint.bcd]
hep48_72[["HTO"]] <- CreateAssayObject(counts = rawsample)
hep48_72 <- NormalizeData(object = hep48_72,assay = "HTO",normalization.method = "CLR")


# 0h sample
zt6a_raw <- read.table(file = "ZT06A.txt.gz",header = T,sep = "\t")
zt6b_raw<- read.table(file = "ZT06B.txt.gz",header = T,sep = "\t")
 #remove deuplicate gene
zt6a_raw <- zt6a_raw[-9718,]
zt6b_raw <- zt6b_raw[-9718,]
rownames(zt6a_raw) <- zt6a_raw$gene_names
zt6a_raw$gene_names <- NULL
rownames(zt6a_raw) <- str_to_title(rownames(zt6a_raw))
rownames(zt6b_raw) <- zt6b_raw$gene_names
zt6b_raw$gene_names <- NULL
rownames(zt6b_raw) <- str_to_title(rownames(zt6b_raw))

zt6a <- CreateSeuratObject(counts = zt6a_raw,project = "zt6a",min.cells = 3,min.features = 100)
zt6b <- CreateSeuratObject(counts = zt6b_raw,project = "zt6b",min.cells = 3,min.features = 100)
hep0 <- merge(x = zt6a,y = zt6b,project="hep0")


# merge sample
m1_48h <- subset(hep48_72,HTO_maxID=="m1-48h")
m2_48h <- subset(hep48_72,HTO_maxID=="m2-48h")
m1_72h <- subset(hep48_72,HTO_maxID=="m1-72h")
m2_72h <- subset(hep48_72,HTO_maxID=="m2-72h")

m1_48h_count <- m1_48h@assays$RNA@counts
m2_48h_count <- m2_48h@assays$RNA@counts
m1_72h_count <- m1_72h@assays$RNA@counts
m2_72h_count <- m2_72h@assays$RNA@counts

rownames(m1_48h_count) <- str_to_title(rownames(m1_48h_count))
rownames(m1_72h_count) <- str_to_title(rownames(m1_72h_count))
rownames(m2_48h_count) <- str_to_title(rownames(m2_48h_count))
rownames(m2_72h_count) <- str_to_title(rownames(m2_72h_count))

m1_48h <- CreateSeuratObject(counts = m1_48h_count,project = "m1_48h")
m2_48h <- CreateSeuratObject(counts = m2_48h_count,project = "m2_48h")
m1_72h <- CreateSeuratObject(counts = m1_72h_count,project = "m1_72h")
m2_72h <- CreateSeuratObject(counts = m2_72h_count,project = "m2_72h")

hep48_72 <- merge(x = m1_48h,y = c(m2_48h,m1_72h,m2_72h),
                  add.cell.ids = c("m1_48h","m2_48h","m1_72h","m2_72h"),project="hep48_72")

remove(m1_48h,m1_48h_count,m1_72h,m1_72h_count,m2_48h,m2_48h_count,m2_72h,m2_72h_count)
remove(zt6a,zt6a_raw,zt6b,zt6b_raw)

merge.list <- list(hep0,hep48_72)
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


hep0_48_72h <- IntegrateData(anchorset = obj.anchors,normalization.method = "LogNormalize")

hep_m1 <- subset(hep0_48_72h,orig.ident=="m1_48h" |orig.ident=="m1_72h" |
                   orig.ident=="zt6a" |orig.ident=="zt6b")

saveRDS(hep_m1,file = "hep_m1_ver1")

# Memory can be safely released at this point.


# Quality Control
hep_all <- readRDS("hep_m1_ver1")
gene_all <- rownames(hep_all@assays$RNA@data)
DefaultAssay(hep_all) <- "RNA"
hep_all@meta.data$percent_mt<- as.vector(PercentageFeatureSet(hep_all,pattern = c("^Mt-"))[,1])
 # please check sample quality
 # The following are our quality control standards.
hep_all <- subset(hep_all,nCount_RNA>=1000 & nCount_RNA<=3500 & percent_mt<=20)
VlnPlot(hep_all,features = c("nCount_RNA","nFeature_RNA","percent_mt"),pt.size = 0)

# Clustering
hep_all <- FindVariableFeatures(hep_all, selection.method = "vst", nfeatures = 2000)
gene_all <- rownames(hep_all)
hep_all <- ScaleData(hep_all, features = gene_all)
hep_all <- RunPCA(hep_all, features = VariableFeatures(object = hep_all))

hep_all@meta.data$Time <- "0h"
m <- hep_all@meta.data$orig.ident=="m1_48h"
hep_all@meta.data$Time[m] <- "48h"
m <- hep_all@meta.data$orig.ident=="m1_72h"
hep_all@meta.data$Time[m] <- "72h"

hep_all <- FindNeighbors(object = hep_all,dims = 1:20)
hep_all <- FindClusters(object = hep_all,resolution = 1.2)
hep_all <- RunUMAP(object = hep_all,dims = 1:20)

saveRDS(hep_all,"hep_m1_ver3")

# Memory can be safely released at this point.










