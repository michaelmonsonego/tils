

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
T_cells = readRDS("objects/tils_all_.45_integrgate_annotated.rds")
DimPlot(T_cells, reduction = "tsne",repel = T, label = TRUE, label.size = 5)
DefaultAssay(T_cells) <- "RNA" #M# important before subseting! otherwise doesnt find cd8

#M# all cells with CD4 > 0.5 and no CD8
cd8_high_cells <- subset(T_cells, subset = CD4 < 0.005 & (CD8A > 0.5 | CD8B < 0.5))
expression_matrix <- GetAssayData(cd8_high_cells, assay = "RNA", slot = "data")
average_cd4_expression_all_cells <- mean(expression_matrix["CD8A", ])
# The average expression across all cells that express cd8A > 0.5
print(average_cd4_expression_all_cells)

table(cd8_high_cells$orig.ident)
table(cd8_high_cells$seurat_clusters)

cd8s <- cd8_high_cells
cd8s = NormalizeData(cd8s, assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000) 
cd8s = FindVariableFeatures(cd8s, assay = "RNA", selection.method = "vst", nfeatures = 2000) 

DefaultAssay(cd8s) <- "integrated"

cd8s = ScaleData(cd8s, verbose = FALSE)
cd8s = RunPCA(cd8s, npcs = 30, verbose = FALSE)
DimPlot(cd8s, reduction = "pca", group.by="Treatment", cols = c("#FDCDAC", "#B3E2CD"), pt.size=1.1)
ElbowPlot(cd8s, ndims = 30) + theme_classic()
ggsave(file = "figures/data_prep_cd8/elbow_plot.png", dpi=300, width=10, height=10)

cd8s <- FindNeighbors(cd8s, reduction = "pca", dims = 1:14) #M# choose ?
cd8s = RunTSNE(cd8s, dims = 1:14)
#M# cd8s = RunUMAP(cd8s, dims = 1:18)

res_seq <- c(.05,.3, .35, .4, .45, .5)
for(res in res_seq){
  cd8s.res_test <- FindClusters(cd8s, resolution = res)
  
  res_tSNE  <- DimPlot(cd8s.res_test, reduction = "tsne",
                       repel = T, label = TRUE, label.size = 5) +
    theme(legend.position = "none") + 
    plot_annotation(title = paste("Res of", res))
  assign(paste0("tSNE_",res), res_tSNE)
}

tSNE_ls <- list(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5)
all_tSNE <- plot_grid(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5, ncol = 3) 
ggsave(all_tSNE,filename = 'figures/data_prep_cd8/tils_0_14_dim_res_test.png', dpi=300, height=10, width=16)

cd8s <- FindClusters(cd8s, resolution = 0.3)
saveRDS(cd8s, file = "objects/cd8_.3.rds") #M# change if adjusted
cd8s <- readRDS("objects/cd8_.3.rds")




# Lineage identification
# T cells
Ts <- c('Cd4', 'Cd8a', 'Cd8b', 'Cd3d','Cd3e','Cd3g',"FOXP3")
Ts <- toupper(Ts)
FeaturePlot(cd8s, features = Ts, order=TRUE,pt.size=1, reduction="tsne", ncol=3)
ggsave(file = "figures/data_prep_cd8/features_T.png", dpi=300, width=30, height=20)


cd8s.markers = FindAllMarkers(cd8s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")

Top50Markers_T =
  cd8s.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC)
write_csv(Top50Markers_T, "Excels/cd8_analysis_DE genes.csv")

cd8s.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_T
DoHeatmap(cd8s, features = top10_T$gene) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))
ggsave(file = "figures/data_prep_cd8/T cells markers - heatmap.png", dpi=300, width=12, height=20)








