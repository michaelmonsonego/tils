
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
library(ggrepel)
library(fgsea)


setwd("D:/Michael/git_check/tils")

T_cells = readRDS("objects/tils_all_.45_integrgate_annotated (2).rds")
table(T_cells$orig.ident)
table(T_cells$seurat_clusters)
levels(Idents(T_cells))

DimPlot(T_cells, reduction = "tsne",repel = T, label = TRUE, label.size = 5)
DefaultAssay(T_cells) <- "RNA" #M# important before subseting! otherwise doesnt find cd8
# pipeline-----------
prolif <- subset(T_cells, idents = c("5" = "5_CD8_cell_cycle",
                                     "7" = "7_CD8_proliferation_like",
                                     "8" = "8_CD8_proliferation")) #M# also subset out CD8 cells
table(prolif$orig.ident)
table(prolif$seurat_clusters)


prolif = NormalizeData(prolif, assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000) 
prolif = FindVariableFeatures(prolif, assay = "RNA", selection.method = "vst", nfeatures = 2000) 

DefaultAssay(prolif) <- "integrated"

prolif = ScaleData(prolif, verbose = FALSE)
prolif = RunPCA(prolif, npcs = 30, verbose = FALSE)
DimPlot(prolif, reduction = "pca", group.by="Treatment", cols = c("#FDCDAC", "#B3E2CD"), pt.size=1.1)
ElbowPlot(prolif, ndims = 30) + theme_classic()
ggsave(file = "figures/prolif/elbow_plot.png", dpi=300, width=10, height=10)

prolif <- FindNeighbors(prolif, reduction = "pca", dims = 1:10) #M# choose ?
prolif = RunTSNE(prolif, dims = 1:10)
#M# prolif = RunUMAP(prolif, dims = 1:10)

res_seq <- c(.05,.3, .35, .4, .45, .5)
for(res in res_seq){
  prolif.res_test <- FindClusters(prolif, resolution = res)
  
  res_tSNE  <- DimPlot(prolif.res_test, reduction = "tsne",
                       repel = T, label = TRUE, label.size = 5) +
    theme(legend.position = "none") + 
    plot_annotation(title = paste("Res of", res))
  assign(paste0("tSNE_",res), res_tSNE)
}

tSNE_ls <- list(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5)
all_tSNE <- plot_grid(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5, ncol = 3) 
ggsave(all_tSNE,filename = 'figures/prolif/tils_0_14_dim_res_test.png', dpi=300, height=10, width=16)

prolif <- FindClusters(prolif, resolution = 0.3)
# saveRDS(prolif, file = "objects/prolif.rds") #M# change if adjusted
prolif <- readRDS("objects/prolif.rds")



DefaultAssay(prolif) <- "RNA"
Ts <- c('Cd4', 'Cd8a', 'Cd8b', 'Cd3d','Cd3e','Cd3g',"FOXP3")
Ts <- toupper(Ts)
FeaturePlot(prolif, features = Ts, order=TRUE,pt.size=1, reduction="tsne", ncol=3)
ggsave(file = "figures/prolif/features_T.png", dpi=300, width=30, height=20)

prolif.markers = FindAllMarkers(prolif, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")

Top50Markers_T =
  prolif.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC)
write_csv(Top50Markers_T, "Excels/prolif_analysis_DE genes.csv")

#M# find DE genes between responders to non-responders
DE_genes <- FindMarkers(prolif, ident.1 = "Responder", ident.2 = "Non_Responder", group.by = "Treatment", assay = "RNA")
DE_genes <- rownames_to_column(DE_genes, var = "gene")
write_csv(DE_genes, "excels/DE_genes_res_nonRes_prolif.csv") 
Top100Markers <- DE_genes %>% 
  top_n(n = 100, wt = avg_log2FC) %>% 
  as.data.frame
write_csv(Top100Markers, "excels/Top100Markers_per_treatment_prolif.csv") 









































