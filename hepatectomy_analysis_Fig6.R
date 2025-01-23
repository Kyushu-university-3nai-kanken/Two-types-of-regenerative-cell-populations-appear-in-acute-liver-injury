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



# readfile
hep1 <- hep1 <- readRDS("./make_obj/hepatectomy_redu1_ver2")
colnames(hep1@meta.data)[1] <- "Time"
 # Please load the objects created during the APAP mouse analysis
hep_origi <- readRDS("./make_obj/hep_m1_ver3")


# Fig.S6A
DimPlot(hep1,reduction = "umap",group.by = "seurat_clusters")+
  DimPlot(hep1,reduction = "umap",group.by = "Time")


# Calculate the Cluster 7 score.
gene_cls7_df <- FindMarkers(hep_origi,ident.1 = c(7,13))
cls7_df <- gene_cls7_df[gene_cls7_df$p_val_adj<1,]
gene_cls7 <- rownames(cls7_df[cls7_df$avg_log2FC>1,])
gene_cls7 <- gene_cls7[gene_cls7 %in% rownames(hep1@assays$RNA@data)]
gene_cls7 <- list(gene_cls7)
hep1 <- AddModuleScore(hep1,features = gene_cls7,name = "cls7_marker")


# Fig.S6B
colnames(hep1@meta.data)[7] <- "cls7_marker"
VlnPlot(hep1,features = "cls7_marker",pt.size = 0)
FeaturePlot(hep1,features = "cls7_marker")+
  theme(axis.text = element_text(size = 18))+
  ggtitle("cluster7 score")+
  scale_color_gradient2(low = "lightgrey",mid = "#F6F6F6",high = "blue",midpoint = 0)


# Fig.S6C
FeaturePlot(hep1,features = c("Malat1","Neat1","Trp53inp1","Jun","Egfr"))