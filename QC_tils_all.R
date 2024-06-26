library(Seurat)
library(plyr)
library(tidyverse)  
library(patchwork)
library(stringr)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(cowplot)
library(future)
library(ggsignif)
library(DT)
library(Cairo)
library(snakecase)
library(scCustomize)

8482 -> saved.seed
set.seed(saved.seed)

theme_set(theme_bw(base_size = 14))
memory.limit(size=10^9)

#M# human rna data integration - sent from keren
setwd("D:/Michael/git_check/tils")

# set up data ----
#Sample_Names_vec <- list.files(path ='/Users/kerenreshef/Desktop/TAU/PhD/Madi_lab/Pemphigus/scRNA-seq/exp.2/samples_2exp',full.names = F)

#seurat_object_list <- list()

# Create Seurat objects in a loop
#for (sample in Sample_Names_vec) {
 # seurat_object_list[[sample]] <- CreateSeuratObject(counts =  Read10X(data.dir = paste0('/Users/kerenreshef/Desktop/TAU/PhD/Madi_lab/Pemphigus/scRNA-seq/exp.2/samples_2exp/',sample,"/filtered_feature_bc_matrix/")),project = sample)
#}

# Load the datasets
seurat_object_list <- list()

N1_data = Read10X(data.dir = "G:/michael/tils/Tils/N1_exp72/outs/filtered_feature_bc_matrix")
N1 = CreateSeuratObject(counts = N1_data, project = "N1")
N1$Sample = "N1"
N1$Treatment = "Non_Responder"
seurat_object_list[["N1"]] <- N1

N2_data = Read10X(data.dir = "G:/michael/tils/Tils/N2_exp72/outs/filtered_feature_bc_matrix")
N2 = CreateSeuratObject(counts = N2_data, project = "N2")
N2$Sample = "N2"
N2$Treatment = "Non_Responder"
seurat_object_list[["N2"]] <- N2

N3_data = Read10X(data.dir = "G:/michael/tils/Tils/N3_exp72/outs/filtered_feature_bc_matrix")
N3 = CreateSeuratObject(counts = N3_data, project = "N3")
N3$Sample = "N3"
N3$Treatment = "Non_Responder"
seurat_object_list[["N3"]] <- N3

N4_data = Read10X(data.dir = "G:/michael/tils/Tils/N4_exp72/outs/filtered_feature_bc_matrix")
N4 = CreateSeuratObject(counts = N4_data, project = "N4")
N4$Sample = "N4"
N4$Treatment = "Non_Responder"
seurat_object_list[["N4"]] <- N4

R1_data = Read10X(data.dir = "G:/michael/tils/Tils/R1_exp72/outs/filtered_feature_bc_matrix")
R1 = CreateSeuratObject(counts = R1_data, project = "R1")
R1$Sample = "R1"
R1$Treatment = "Responder"
seurat_object_list[["R1"]] <- R1

R2_data = Read10X(data.dir = "G:/michael/tils/Tils/R2_exp72/outs/filtered_feature_bc_matrix")
R2 = CreateSeuratObject(counts = R2_data, project = "R2")
R2$Sample = "R2"
R2$Treatment = "Responder"
seurat_object_list[["R2"]] <- R2

R3_data = Read10X(data.dir = "G:/michael/tils/Tils/R3_exp72/outs/filtered_feature_bc_matrix")
R3 = CreateSeuratObject(counts = R3_data, project = "R3")
R3$Sample = "R3"
R3$Treatment = "Responder"
seurat_object_list[["R3"]] <- R3

R4_data = Read10X(data.dir = "G:/michael/tils/Tils/R4_exp72/outs/filtered_feature_bc_matrix")
R4 = CreateSeuratObject(counts = R4_data, project = "R4")
R4$Sample = "R4"
R4$Treatment = "Responder"
seurat_object_list[["R4"]] <- R4

#downnsample
seurat_object_sub_list <- list()
Sample_Names_vec <- c("N1","N2","N3","N4","R1","R2","R3","R4")

for (sample in Sample_Names_vec) {
  seurat_object_sub_list[[sample]] <- seurat_object_list[[sample]][, sample(colnames(seurat_object_list[[sample]]), size =6900, replace=F)]
}

#mito and ribo
for (sample in Sample_Names_vec){
  seurat_object_sub_list[[sample]][["percent.mt"]] <- PercentageFeatureSet(seurat_object_sub_list[[sample]], pattern = "^MT-")
  seurat_object_sub_list[[sample]][["percent.ribo"]] <- PercentageFeatureSet(seurat_object_sub_list[[sample]], pattern =  "^RP[SL]")
}
# Visualize QC metrics as a violin plot
for (sample in Sample_Names_vec){
  VlnPlot(seurat_object_sub_list[[sample]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo") ,ncol = 4)
  ggsave(filename = paste0(sample,'_QC_Violin_Before.png'), path = 'figures/integration/before', dpi=300, height=7, width=12, device = 'png')
}

#subset
for (sample in Sample_Names_vec){ #m# todo : choose values : make sure with friends
  seurat_object_sub_list[[sample]] <- subset(seurat_object_sub_list[[sample]],
                            nFeature_RNA > 200 &
                              nFeature_RNA < 5500 &
                              nCount_RNA > 200 &
                              nCount_RNA < 40000 &
                              percent.mt < 15 &
                              percent.ribo < 45
)
}

# Visualize QC metrics as a violin plot
for (sample in Sample_Names_vec){
  VlnPlot(seurat_object_sub_list[[sample]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo") ,ncol = 4)
  ggsave(filename = paste0(sample,'_QC_Violin_After.png'), path = 'figures/integration/after', dpi=300, height=7, width=12, device = 'png')
}

#norm
for (sample in Sample_Names_vec){
  seurat_object_sub_list[[sample]] <- NormalizeData(seurat_object_sub_list[[sample]], normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_object_sub_list[[sample]] <- FindVariableFeatures(seurat_object_sub_list[[sample]], selection.method = "vst", nfeatures = 3000)
}  

#pre-integrate
features <- SelectIntegrationFeatures(object.list =seurat_object_sub_list)

#integration
immune.anchors <- FindIntegrationAnchors(object.list = seurat_object_sub_list, anchor.features = features)
#======from this part the analysis was done on the server===== 
# saveRDS(immune.anchors, file = "objects/immune.anchors.tils_all.rds") 
immune.anchors <- readRDS("objects/immune.anchors.tils_all.rds")
immune.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
ep <- ElbowPlot(immune.combined, ndims = 30) 
ggsave(ep, filename = 'figures/integration/Elbow.integrate.jpg', dpi=300, height=7, width=12)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:15) #M# 15 works for me
immune.combined = RunTSNE(immune.combined, dims = 1:15)
merged = RunUMAP(immune.combined, dims = 1:15)

res_seq <- c(.25,.3, .35, .4, .45, .5)
for(res in res_seq){
  immune.combined.res_test <- FindClusters(immune.combined, resolution = res)
  
  res_tSNE  <- DimPlot(immune.combined.res_test, reduction = "tsne",
                       repel = T, label = TRUE, label.size = 5) +
    theme(legend.position = "none") + 
    plot_annotation(title = paste("Res of", res))
  assign(paste0("tSNE_",res), res_tSNE)
}

tSNE_ls <- list(tSNE_0.25,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5)
all_tSNE <- plot_grid(tSNE_0.25,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5, ncol = 3) 
ggsave(all_tSNE,filename = 'figures/integration/tils_15_dim_res_test.png', dpi=300, height=10, width=16, device = 'png') #M# todo : change name

#M# same for umap (if needed) : 
for(res in res_seq){
  immune.combined.res_test <- FindClusters(immune.combined, resolution = res)
  
  res_umap  <- DimPlot(immune.combined.res_test, reduction = "umap",
                       repel = T, label = TRUE, label.size = 5) +
    theme(legend.position = "none") + 
    plot_annotation(title = paste("Res of", res))
  assign(paste0("umap_",res), res_umap)
}
umsp_ls <- list(umap_0.25,umap_0.3,umap_0.35,umap_0.4,umap_0.45,umap_0.5)
all_umap <- plot_grid(umap_0.25,umap_0.3,umap_0.35,umap_0.4,umap_0.45,umap_0.5, ncol = 3) 
ggsave(all_umap,filename = 'figures/integration/tils_8_dim_res_test_umap.png', dpi=300, height=10, width=16, device = 'png') #M# todo : change name





immune.combined <- FindClusters(immune.combined, resolution = .45) #M# choose res here : important for rest of analysis
#M# saveRDS(immune.combined, file = "objects/tils_all_.45_integrgate.rds")
immune.combined <- readRDS("objects/tils_all_.45_integrgate.rds")
table(immune.combined$orig.ident)

# immune.combined <- JoinLayers(immune.combined, assay = "RNA") #M# this did not work : no layers in my seurat
allmarkers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
Top50Markers <- allmarkers %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame %>% 
  arrange(cluster, -avg_log2FC)

write_csv(Top50Markers, "excels/Top50Markers_perClust.tils.45integrgate.csv") 





