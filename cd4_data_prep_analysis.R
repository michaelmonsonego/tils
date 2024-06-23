

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
# cd4s <- subset(T_cells, idents = c("0" = "0_activated_memory_CD4")) #M# also subset out CD8 cells
# table(cd4s$orig.ident)
# table(cd4s$seurat_clusters)
#M# save cd8 to move to cluster 1 cd8 analysis
# cd8 <- subset(cd4s, subset = CD8A > 0.05 & CD8B > 0.05)
# table(cd8$orig.ident)
# table(cd8$seurat_clusters)
# saveRDS(cd8, "objects/cd8_from_clus_0.rds")

#M# clean cd8 from this analysis
#cd4s <- subset(cd4s, subset = CD8A < 0.005 & CD8B < 0.005)
# table(cd4s$orig.ident)
# table(cd4s$seurat_clusters)
# DefaultAssay(cd4s)

#M# check how many cells from cd4s express cd4 > 0.5
# check1 <- subset(cd4s, subset = CD4 > 0.5)
# table(check1$orig.ident)
# table(check1$seurat_clusters)



#M# all cells with CD4 > 0.5 and no CD8
cd4_high_cells <- subset(T_cells, subset = CD4 > 0.5 & CD8A < 0.005 & CD8B < 0.005)
expression_matrix <- GetAssayData(cd4_high_cells, assay = "RNA", slot = "data")
average_cd4_expression_all_cells <- mean(expression_matrix["CD4", ])
# The average expression across all cells that express cd4 > 0.5
print(average_cd4_expression_all_cells)

table(cd4_high_cells$orig.ident)
table(cd4_high_cells$seurat_clusters)


#M# all cells with CD4 > 1 
# cd4_high_cells <- subset(T_cells, subset = CD4 > 1 & CD8A < 0.005 & CD8B < 0.005)
# expression_matrix <- GetAssayData(cd4_high_cells, assay = "RNA", slot = "data")
# average_cd4_expression_all_cells <- mean(expression_matrix["CD4", ])
# The average expression across all cells that express cd4 > 1
# print(average_cd4_expression_all_cells)

# table(cd4_high_cells$orig.ident)
# table(cd4_high_cells$seurat_clusters)



#M# this is for taking all cluster 0 cells and all cd4 high expressing cells
# unique_cells <- unique(c(Cells(cd4s), Cells(cd4_high_cells)))
# cd4s <- subset(T_cells, cells = unique_cells)
# table(cd4s$orig.ident)
# table(cd4s$seurat_clusters)
# DefaultAssay(cd4s)

#M# addition of cd4 cells from other object  
#cd4_from_1 <- readRDS("objects/cd4_from_clus_1.rds")
#table(cd4_from_1$orig.ident)
#table(cd4_from_1$seurat_clusters)
#DefaultAssay(cd4_from_1) <- "RNA"
# cd4s = merge(cd4s, cd4_from_1) #M# this didnt work. 



#M# 23.6 after conversation with Ayelet we chose this option : all cd4 high cells, without taking all cells from cluster 0 first
cd4s <- cd4_high_cells
cd4s = NormalizeData(cd4s, assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000) 
cd4s = FindVariableFeatures(cd4s, assay = "RNA", selection.method = "vst", nfeatures = 2000) 

#M# for addition of cells from other object
#cd4_from_1 = NormalizeData(cd4_from_1, assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000)  
#cd4_from_1 = FindVariableFeatures(cd4_from_1, assay = "RNA", selection.method = "vst", nfeatures = 2000) 

#features <- SelectIntegrationFeatures(object.list =list(cd4s, cd4_from_1))

#anchors <- FindIntegrationAnchors(object.list = list(cd4s, cd4_from_1), anchor.features = features)
#cd4s <- IntegrateData(anchorset = anchors)

#cd4s$Sample = factor(cd4s$Sample, levels = c("N1","N2","N3","N4","R1","R2","R3","R4"))
#cd4s$orig.ident= factor(cd4s$orig.ident, levels=c("N1","N2","N3","N4","R1","R2","R3","R4"))

#table(cd4s$orig.ident)
#table(cd4s$seurat_clusters)
#table(cd4s$Sample)
#table(cd4s$Treatment)










DefaultAssay(cd4s) <- "integrated"

cd4s = ScaleData(cd4s, verbose = FALSE)
cd4s = RunPCA(cd4s, npcs = 30, verbose = FALSE)
DimPlot(cd4s, reduction = "pca", group.by="Treatment", cols = c("#FDCDAC", "#B3E2CD"), pt.size=1.1)
ElbowPlot(cd4s, ndims = 30) + theme_classic()
ggsave(file = "figures/data_prep_cd4/elbow_plot.png", dpi=300, width=10, height=10)

cd4s <- FindNeighbors(cd4s, reduction = "pca", dims = 1:14) #M# choose ?
cd4s = RunTSNE(cd4s, dims = 1:14)
#M# cd4s = RunUMAP(cd4s, dims = 1:18)

res_seq <- c(.05,.3, .35, .4, .45, .5)
for(res in res_seq){
  cd4s.res_test <- FindClusters(cd4s, resolution = res)
  
  res_tSNE  <- DimPlot(cd4s.res_test, reduction = "tsne",
                       repel = T, label = TRUE, label.size = 5) +
    theme(legend.position = "none") + 
    plot_annotation(title = paste("Res of", res))
  assign(paste0("tSNE_",res), res_tSNE)
}

tSNE_ls <- list(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5)
all_tSNE <- plot_grid(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5, ncol = 3) 
ggsave(all_tSNE,filename = 'figures/data_prep_cd4/tils_0_14_dim_res_test.png', dpi=300, height=10, width=16)

cd4s <- FindClusters(cd4s, resolution = 0.3)
saveRDS(cd4s, file = "objects/cd4_.3.rds") #M# change if adjusted
cd4s <- readRDS("objects/cd4_.3.rds")




# Lineage identification
# T cells
Ts <- c('Cd4', 'Cd8a', 'Cd8b', 'Cd3d','Cd3e','Cd3g',"FOXP3")
Ts <- toupper(Ts)
FeaturePlot(cd4s, features = Ts, order=TRUE,pt.size=1, reduction="tsne", ncol=3)
ggsave(file = "figures/data_prep_0/features_T.png", dpi=300, width=30, height=20)



cd4s.markers = FindAllMarkers(cd4s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")

Top50Markers_T =
  cd4s.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC)
write_csv(Top50Markers_T, "Excels/cd4_analysis_DE genes.csv")

cd4s.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_T
DoHeatmap(cd4s, features = top10_T$gene) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))
ggsave(file = "figures/data_prep_cd4/cd4_heatmap.png", dpi=300, width=12, height=20)








