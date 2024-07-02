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


#M# clus 1 sigs
clus1 = readRDS("objects/tils_1_.5.rds")

#M# sig up genes from paper
persistance_up <- c('KLRB1','ZNF683','ITGB1','C1orf162','IL7R','LIME1', "S1PR1", "TIMP1", "C10orf54", "TBXAS1", "KLF2", "LTB", "UBXN11", "CD40LG", "AMICA1", "FAM65B", "VCL", "RASA3", "SCML4", "MYC", "P2RY8")
persistance_up_sig <- make_signature(clus1, persistance_up,idents = c("Responder","Non_Responder"), format_flag = FALSE)
persistance_up_sig_obj <- persistance_up_sig[[1]]
persistance_up_sig[[3]]
ggsave("figures/clus1/co_stim_sig_box.png", width = 20, height = 20, units = "cm")
persistance_up_sig[[2]]
ggsave("figures/clus1/co_stim_sig_tsne.png", width = 20, height = 20, units = "cm")

# violin plot per treatment
VlnPlot(persistance_up_sig_obj, features = "SigUint", group.by = "Treatment", pt.size = 0)+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+
  geom_boxplot(alpha=0.3,show.legend = FALSE)
ggsave("figures/clus1/violin_sig_co_stimulatory_by_treatment.png", width = 20, height = 20, units = "cm")

# violin plot per treatment per cluster
gglist <-  list()
name <- "persistance up regulated"
for(clus in levels(persistance_up_sig_obj$seurat_clusters)){
  print(clus)
  obj <- subset(persistance_up_sig_obj, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'SigUint') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus1/persistance_up_sig_per_cluster_treatment.png", dpi=300, width=16, height=10)




























































