

library(Seurat)
library(tidyverse)
library(sleepwalk)
library(gridExtra)
library(plyr)
library(writexl)
library(patchwork)
library(RColorBrewer)
library(scales)
library(matlab)
# library(fgsea)
library(stringr)
library(reticulate)
library(dplyr)
library(future)
library(reshape2)
library(ggplot2)
library(grid)
library(cowplot)
library(ggsignif)
library(DT)
library(Cairo)
library(snakecase)
library(glue)


setwd("D:/Michael/git_check/tils")
# setwd("C:/Users/michael monsonego/Documents/tils") #M# update between computers(git)

use_python((r"(C:\Users\user\AppData\Local\Programs\Python\Python311)"), required = TRUE)
source("//asafmadilabfs/asafmadilab$/michael/Madi lab/Signature_propo/signature_utils.R")
source('D:/Michael/git_check/3-groups/sig_genes.R')

# Signature score function
SignatureScore <- function(object, name){
  merged_responder <- subset(object, subset = Treatment == "Responder")
  merged_NON_responder <- subset(object, subset = Treatment == "Non_Responder")
  
  
  pval <- unlist(wilcox.test(unlist(merged_NON_responder[[paste0(name)]]), 
                             unlist(merged_responder[[paste0(name)]]), 
                             alternative = "two.sided"))[2]
  
  object$normalized = (object[[name]]-min(object[[name]]))/(max(object[[name]])-min(object[[name]]))
  
  # print(pval)
  # y.max <-  1.14*max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  y.max <- max(object[['normalized']])
  annot.max.y <- y.max*0.98
  
  label <- as.character(symnum(as.numeric(pval), corr = FALSE, na = FALSE,
                               cutpoints = c(0,0.0001, 0.001, 0.01, 0.05,1),
                               symbols = c('****',"***", "**", "*", "ns")))
  # print(label)
  print(paste0("Pvalue = ",pval, " == ", label))
  
  p <- VlnPlot(object, features='normalized', group.by = "Treatment",
               y.max = y.max, pt.size = 0,
               cols= c("#1f77b4","#ff7f0e")) +  theme_classic(base_size = 14) +
    theme(text = element_text(size=18, colour = "black")) + RotatedAxis() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    labs(title = "", y = name,  x="") + theme(legend.position="right") +
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") +
    annotate("text", x = 1.5, y=annot.max.y, label = label, size = 6) +
    scale_y_continuous(labels = comma)
  return (p)
  
}

#M# clus 1 sigs
clus1_with_sigs = readRDS("objects/clus1_with_sigs.rds")
clus_2_with_sigs <- readRDS("objects/tils_2_with_sigs.rds")

#M# sig up genes from paper
persistance_up <- c('KLRB1','ZNF683','ITGB1','C1orf162','IL7R','LIME1', "S1PR1", "TIMP1", "C10orf54", "TBXAS1", "KLF2", "LTB", "UBXN11", "CD40LG", "AMICA1", "FAM65B", "VCL", "RASA3", "SCML4", "MYC", "P2RY8")
persistance_up_sig <- make_signature(clus1_with_sigs, persistance_up,idents = c("Responder","Non_Responder"), format_flag = FALSE)
persistance_up_sig_obj <- persistance_up_sig[[1]]
persistance_up_sig[[3]]
ggsave("figures/clus1_with_sigs/co_stim_sig_box.png", width = 20, height = 20, units = "cm")
persistance_up_sig[[2]]
ggsave("figures/clus1_with_sigs/co_stim_sig_tsne.png", width = 20, height = 20, units = "cm")


#M# test dor addmodulescore() function from vis.rmd
# violin plot per treatment
VlnPlot(clus1_with_sigs, features = "Tregs1", group.by = "Treatment", pt.size = 0)+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+
  geom_boxplot(alpha=0.3,show.legend = FALSE)
ggsave("figures/clus1_with_sigs/violin_sig_co_stimulatory_by_treatment.png", width = 20, height = 20, units = "cm")

# violin plot per treatment per cluster
gglist <-  list()
name <- "il21"

DefaultAssay(clus1_with_sigs) <- "RNA" 

for(clus in levels(clus1_with_sigs$seurat_clusters)){
  print(clus)
  obj <- subset(clus1_with_sigs, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'il21') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus1/il2_sig_per_cluster_treatment.png", dpi=300, width=16, height=10)

DefaultAssay(clus1_with_sigs) <- "integrated"

# violin plot per treatment per cluster
gglist <-  list()
name <- "il2_r1"

DefaultAssay(clus1_with_sigs) <- "RNA" 

for(clus in levels(clus1_with_sigs$seurat_clusters)){
  print(clus)
  obj <- subset(clus1_with_sigs, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'il2_r1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus1/il2r_sig_per_cluster_treatment.png", dpi=300, width=16, height=10)


# ____

# violin plot per treatment per cluster
gglist <-  list()
name <- "persistance_up1"

DefaultAssay(clus1_with_sigs) <- "RNA" 

for(clus in levels(clus1_with_sigs$seurat_clusters)){
  print(clus)
  obj <- subset(clus1_with_sigs, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'persistance_up1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus1/persistance_up_sig_per_cluster_treatment.png", dpi=300, width=16, height=10)


# ____

# violin plot per treatment per cluster
gglist <-  list()
name <- "CD8_memory1"

DefaultAssay(clus1_with_sigs) <- "RNA" 

for(clus in levels(clus1_with_sigs$seurat_clusters)){
  print(clus)
  obj <- subset(clus1_with_sigs, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'CD8_memory1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus1/CD8_memory_sig_per_cluster_treatment.png", dpi=300, width=16, height=10)


# ___

# violin plot per treatment per cluster
gglist <-  list()
name <- "effector_memory1"

DefaultAssay(clus1_with_sigs) <- "RNA" 

for(clus in levels(clus1_with_sigs$seurat_clusters)){
  print(clus)
  obj <- subset(clus1_with_sigs, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'effector_memory1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus1/effector_memory_sig_per_cluster_treatment.png", dpi=300, width=16, height=10)



DefaultAssay(clus1_with_sigs) <- "integrated"


# ____

#M# same for cluster 2 : 

# violin plot per treatment per cluster
gglist <-  list()
name <- "il21"

DefaultAssay(clus_2_with_sigs) <- "RNA" 

for(clus in levels(clus_2_with_sigs$seurat_clusters)){
  print(clus)
  obj <- subset(clus_2_with_sigs, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'il21') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus2/il2_sig_per_cluster_treatment.png", dpi=300, width=16, height=10)

DefaultAssay(clus_2_with_sigs) <- "integrated"

# violin plot per treatment per cluster
gglist <-  list()
name <- "il2_r1"

DefaultAssay(clus_2_with_sigs) <- "RNA" 

for(clus in levels(clus_2_with_sigs$seurat_clusters)){
  print(clus)
  obj <- subset(clus_2_with_sigs, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'il2_r1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus2/il2r_sig_per_cluster_treatment.png", dpi=300, width=16, height=10)


# ____

# violin plot per treatment per cluster
gglist <-  list()
name <- "persistance_up1"

DefaultAssay(clus_2_with_sigs) <- "RNA" 

for(clus in levels(clus_2_with_sigs$seurat_clusters)){
  print(clus)
  obj <- subset(clus_2_with_sigs, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'persistance_up1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus2/persistance_up_sig_per_cluster_treatment.png", dpi=300, width=16, height=10)


# ____

# violin plot per treatment per cluster
gglist <-  list()
name <- "CD8_memory1"

DefaultAssay(clus_2_with_sigs) <- "RNA" 

for(clus in levels(clus_2_with_sigs$seurat_clusters)){
  print(clus)
  obj <- subset(clus_2_with_sigs, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'CD8_memory1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus2/CD8_memory_sig_per_cluster_treatment.png", dpi=300, width=16, height=10)


# ___

# violin plot per treatment per cluster
gglist <-  list()
name <- "effector_memory1"

DefaultAssay(clus_2_with_sigs) <- "RNA" 

for(clus in levels(clus_2_with_sigs$seurat_clusters)){
  print(clus)
  obj <- subset(clus_2_with_sigs, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'effector_memory1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus2/effector_memory_sig_per_cluster_treatment.png", dpi=300, width=16, height=10)
























































