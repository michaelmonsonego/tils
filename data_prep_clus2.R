

## 1. Selection of cluster 0
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(sleepwalk)
library(gridExtra)
library(future)
library(reshape2)
library(reticulate)
library(ggplot2)
library(plyr)
library(stringr)
library(RColorBrewer)
library(grid)
library(cowplot)
library(ggsignif)
library(DT)
library(Cairo)
library(snakecase)
library(glue)

setwd("D:/Michael/git_check/tils")
# setwd("C:/Users/michael monsonego/Documents/tils") #M# update between computers(git)

x1=4.6
T_cells = readRDS("objects/tils_all_.45_integrgate_annotated_merge_prolif.rds")
DimPlot(T_cells, reduction = "tsne",repel = T, label = TRUE, label.size = 5)
DefaultAssay(T_cells) <- "RNA" #M# important before subseting! otherwise doesnt find cd8
T_cells_2 <- subset(T_cells, idents = c("2" = "2_cytotoxic_CD8")) #M# also subset out CD8 cells
table(T_cells_2$orig.ident)
table(T_cells_2$seurat_clusters)

#M# save cd4 to move to cluster 0 cd4 analysis
# cd4 <- subset(T_cells_2, subset = CD4 > 0.5)
# table(cd4$orig.ident)
# table(cd4$seurat_clusters)
# saveRDS(cd4, "objects/cd4_from_clus_1.rds")

#M# clean cd4 from this analysis
T_cells_2 <- subset(T_cells_2, subset = CD4 < 0.05)###########
table(T_cells_2$orig.ident)
table(T_cells_2$seurat_clusters)
DefaultAssay(T_cells_2) <- "RNA"
#M# merge cd8 cells from cluster 0 to this analysis
#cd8_from_0 <- readRDS("objects/cd8_from_clus_0.rds")
#table(cd8_from_0$orig.ident)
#table(cd8_from_0$seurat_clusters)
#DefaultAssay(cd8_from_0) <- "RNA"

# T_cells_2 = merge(T_cells_2, cd8_from_0) #M# this didnt work. 

T_cells_2 = NormalizeData(T_cells_2, assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000) 
T_cells_2 = FindVariableFeatures(T_cells_2, assay = "RNA", selection.method = "vst", nfeatures = 2000) 

#cd8_from_0 = NormalizeData(cd8_from_0, assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000)  
#cd8_from_0 = FindVariableFeatures(cd8_from_0, assay = "RNA", selection.method = "vst", nfeatures = 2000) 

#features <- SelectIntegrationFeatures(object.list =list(T_cells_2, cd8_from_0))

#anchors <- FindIntegrationAnchors(object.list = list(T_cells_2, cd8_from_0), anchor.features = features)
#T_cells_2 <- IntegrateData(anchorset = anchors)

#T_cells_2$Sample = factor(T_cells_2$Sample, levels = c("N1","N2","N3","N4","R1","R2","R3","R4"))
#T_cells_2$orig.ident= factor(T_cells_2$orig.ident, levels=c("N1","N2","N3","N4","R1","R2","R3","R4"))

#table(T_cells_2$orig.ident)
#table(T_cells_2$seurat_clusters)
#table(T_cells_2$Sample)
#table(T_cells_2$Treatment)


DefaultAssay(T_cells_2) <- "integrated" #M# is this needed here?

T_cells_2 = ScaleData(T_cells_2, verbose = FALSE)
T_cells_2 = RunPCA(T_cells_2, npcs = 30, verbose = FALSE)
DimPlot(T_cells_2, reduction = "pca", group.by="Treatment", cols = c("#FDCDAC", "#B3E2CD"), pt.size=1.1)
ElbowPlot(T_cells_2, ndims = 30) + theme_classic()
ggsave(file = "figures/data_prep_1/elbow_plot.png", dpi=300, width=10, height=10)

T_cells_2 <- FindNeighbors(T_cells_2, reduction = "pca", dims = 1:14) #M# choose ?
T_cells_2 = RunTSNE(T_cells_2, dims = 1:14)
#M# T_cells_2 = RunUMAP(T_cells_2, dims = 1:18)

res_seq <- c(.05,.3, .35, .4, .45, .5)
for(res in res_seq){
  T_cells_2.res_test <- FindClusters(T_cells_2, resolution = res)
  
  res_tSNE  <- DimPlot(T_cells_2.res_test, reduction = "tsne",
                       repel = T, label = TRUE, label.size = 5) +
    theme(legend.position = "none") + 
    plot_annotation(title = paste("Res of", res))
  assign(paste0("tSNE_",res), res_tSNE)
}

tSNE_ls <- list(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5)
all_tSNE <- plot_grid(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5, ncol = 3) 
ggsave(all_tSNE,filename = 'figures/data_prep_2/tils_0_14_dim_res_test.png', dpi=300, height=10, width=16)

T_cells_2 <- FindClusters(T_cells_2, resolution = 0.5)
# saveRDS(T_cells_2, file = "objects/tils_2.rds")
# T_cells_2 <- readRDS("objects/tils_2.rds")




# Lineage identification
### a. T cells
Ts <- c('Cd4', 'Cd8a', 'Cd3d','Cd3e','Cd3g',"PTPRC")
Ts <- toupper(Ts)
FeaturePlot(T_cells_2, features = Ts, order=TRUE,pt.size=1, reduction="tsne", ncol=3)
ggsave(file = "figures/data_prep_2/features_T.png", dpi=300, width=30, height=20)



T_cells_2.markers = FindAllMarkers(T_cells_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")

Top50Markers_T =
  T_cells_2.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC)
write_csv(Top50Markers_T, "Excels/cluster2_analysis_DE genes.csv")

DE_genes <- FindMarkers(T_cells_2, ident.1 = "Responder", ident.2 = "Non_Responder", group.by = "Treatment", assay = "RNA")
DE_genes <- rownames_to_column(DE_genes, var = "gene")
write_csv(DE_genes, "excels/DE_genes_res_nonRes_clus2.csv") 
Top100Markers <- DE_genes %>% 
  top_n(n = 100, wt = avg_log2FC) %>% 
  as.data.frame
write_csv(Top100Markers, "excels/Top100Markers_per_treatment_tils_2.csv") 




T_cells_2.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_T
DoHeatmap(T_cells_2, features = top10_T$gene) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))
ggsave(file = "figures/data_prep_2/T cells markers - heatmap.png", dpi=300, width=12, height=20)









