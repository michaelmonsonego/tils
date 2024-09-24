

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
x1=4.6


# stem like analysis : CD39 & CD69 negative population ------------

VlnPlot(
    T_cells, 
    features = c('ENTPD1', 'CD69'), 
    assay = "RNA", 
    stack = TRUE, 
    flip = TRUE, 
    split.by = "Treatment"
  ) + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 24, face = "italic"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.ticks.y = element_line(size = 0.5),
      strip.text.y = element_text(angle = 0, size = 16, face = "bold")
    ) +
    geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/Tcells/stem_like.png", dpi=300, width=18, height=8)

# important!
DefaultAssay(T_cells) <- "RNA"
dim(T_cells)
# subset for double negative
T_cells_stem_like <- subset(T_cells, subset = ENTPD1 == 0 & CD69 == 0 & CD4 < 0.0000000005)

dim(T_cells_stem_like) # 2367 cells in the double negative population
table(T_cells_stem_like$Treatment)
DimPlot(T_cells_stem_like, reduction = "tsne")

responders <- subset(T_cells_stem_like, subset = Treatment == "Responder")
non_responders <- subset(T_cells_stem_like, subset = Treatment == "Non_Responder")
library(cowplot)
p1 <- DimPlot(responders, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) +
  ggtitle("Responders")
p2 <- DimPlot(non_responders, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) +
  ggtitle("Non Responders")
plot_grid(p1, p2, ncol = 2)
ggsave(file = "figures/Tcells/stem_like_tsne_by_treatment.png", dpi=300, width=18, height=8)


# bar plots 
frque <- table(Idents(T_cells_stem_like), T_cells_stem_like$Treatment)
frque <- as.data.frame(frque)
frque <- frque %>% dplyr::rename(Cluster = Var1, Treatment = Var2)
frque$Treatment <- factor(frque$Treatment)

frque <- ddply(frque, .(Treatment),  transform, percentperstatus   = Freq/sum(Freq))
frque <- ddply(frque, .(Cluster), transform, percentpercluster  = Freq/sum(Freq))
frque <- ddply(frque, .(Cluster), transform, normalizedppstatus = percentperstatus/sum(percentperstatus)*100)
frque <- ddply(frque, .(Treatment),  transform, normalizedpcluster = percentpercluster/sum(percentpercluster))

mfrque <- table(Idents(T_cells_stem_like), T_cells_stem_like$orig.ident)
mfrque <- as.data.frame(mfrque)
mfrque <- mfrque %>% dplyr::rename(Cluster = Var1, patient = Var2)

mfrque <- ddply(mfrque, .(patient),  transform, percentperstatus   = Freq/sum(Freq))
mfrque <- ddply(mfrque, .(Cluster), transform, percentpercluster  = Freq/sum(Freq))
mfrque <- ddply(mfrque, .(Cluster), transform, normalizedppstatus = percentperstatus/sum(percentperstatus)*100)
mfrque <- ddply(mfrque, .(patient),  transform, normalizedpcluster = percentpercluster/sum(percentpercluster))


p3 <- ggplot(mfrque, aes(x = Cluster, y = normalizedppstatus, fill = patient)) +
  theme_minimal() +
  geom_bar(stat = "identity", width = 0.7,position="dodge") +
  scale_fill_manual(values= c("#CFBAE1", "#A8A5CA", "#8291B3", "#5C7D9D", "#FF6F61", "#D78459", "#AF9A52", "#88B04B"))+
  ggtitle("Normalized Cluster Counts") +
  ylab("Normalized Proportions in %") +
  theme(
    axis.text.x  = element_text(angle = 67, hjust = 1, size = 11, color = "black"),
    axis.text.y  = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13, color = "black"),
    title = element_text(size = 13),
    legend.justification = "top",
    panel.border = element_blank()) +
  geom_hline(yintercept = 20, linetype='dotted')


p1 <- ggplot(frque, aes(x = Cluster, y = normalizedppstatus, fill = Treatment)) +
  theme_minimal() +
  geom_bar(stat = "identity", width = 0.7,position="dodge") +
  scale_fill_manual(values= c("#ff7f0e", "#1f77b4"))+
  ggtitle("Normalized Cluster Counts") +
  ylab("Normalized Proportions in %") +
  theme(
    axis.text.x  = element_text(angle = 67, hjust = 1, size = 11, color = "black"),
    axis.text.y  = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13, color = "black"),
    title = element_text(size = 13),
    legend.justification = "top",
    panel.border = element_blank()) +
  geom_hline(yintercept = 50, linetype='dotted') 

p2<-ggplot(frque, aes(x= Cluster, y= Freq, fill = Treatment)) +
  theme_minimal() +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values= c("#ff7f0e", "#1f77b4"))+
  ggtitle("Cluster Counts") + ylab("Cell counts\n") +
  theme(
    axis.text.x  = element_text(angle = 67, hjust = 1, size = 11, color = "black"),
    axis.text.y  = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13, color = "black"),
    title = element_text(size = 13),
    legend.justification = "top",
    panel.border = element_blank())

p2+p1+plot_layout(guides = "collect")+plot_annotation(tag_levels = 'A') 
ggsave(file = "figures/T_cells_stem_like/unnormalized_and_normalized_contribution_per_patient.png", dpi=300, width=14, height=10)

p3
ggsave(file = "figures/T_cells_stem_like/normalized_cluster_counts.png", dpi=300, width=14, height=10)









# cluster 2 DE genes visualisation vln -----------------
DefaultAssay(clus2) <- "RNA"

allmarkers <- read.delim("excels/allmarkers.csv", sep = ",", header = T, quote="\"", check.names=FALSE)
clus2.markers = allmarkers %>%
  filter(cluster=="2_cytotoxic_CD8") %>% # 19.9.24
  top_n(n = 80, wt = avg_log2FC) %>%
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC) %>% 
  pull(gene)

#M# in clus2
i <- 1
while(i<length(clus2.markers)){
  VlnPlot(
    clus2, 
    features = c(clus2.markers[i], clus2.markers[i+1], clus2.markers[i+2], clus2.markers[i+3], clus2.markers[i+4], clus2.markers[i+5], clus2.markers[i+6], clus2.markers[i+7]), 
    assay = "RNA", 
    stack = TRUE, 
    flip = TRUE, 
    split.by = "Treatment"
  ) + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 24, face = "italic"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.ticks.y = element_line(size = 0.5),
      strip.text.y = element_text(angle = 0, size = 16, face = "bold")
    ) +
    geom_boxplot(alpha = 0.3, show.legend = FALSE)
  ggsave(file = paste0("figures/clus2/de_gene_search_by_cluster/","interesting_DE_run_", i,  ".png"), dpi=300, width=16, height=12)
  i <- i+8
}
i <- 1

#M# in all cells
i <- 1
while(i<length(clus2.markers)){
  VlnPlot(
    T_cells, 
    features = c(clus2.markers[i], clus2.markers[i+1], clus2.markers[i+2], clus2.markers[i+3], clus2.markers[i+4], clus2.markers[i+5], clus2.markers[i+6], clus2.markers[i+7]), 
    assay = "RNA", 
    stack = TRUE, 
    flip = TRUE, 
    split.by = "Treatment"
  ) + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 24, face = "italic"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.ticks.y = element_line(size = 0.5),
      strip.text.y = element_text(angle = 0, size = 16, face = "bold")
    ) +
    geom_boxplot(alpha = 0.3, show.legend = FALSE)
  ggsave(file = paste0("figures/clus2/de_gene_search_by_cluster_allCellsClusters/","all_cells_interesting_DE_run_", i,  ".png"), dpi=300, width=16, height=12)
  i <- i+8
}
i <- 1

# save tsne with cluster number as label for all cells -------------------------
DimPlot(T_cells, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5, label.size = 16)
ggsave(file = "figures/Tcells/tsne_number_labels.png", dpi=300, width=15, height=10)

# info for rons tool ---------------------


cell_annotations <- data.frame(Cell_ID = names(Idents(T_cells)), 
                               Identity = Idents(T_cells), row.names = NULL)
write.csv(cell_annotations, file = "excels/T_cells_idents_annotations.csv", row.names = FALSE)


rna_assay <- GetAssayData(T_cells, assay = "RNA", slot = "counts")
rna_assay_df <- as.data.frame(as.matrix(rna_assay))
write.csv(rna_assay_df, file = "excels/T_cells_RNA_assay_counts.csv")



# merge proliferation clusters(5,7,8) in all cells--------------------------------------------

new.cluster.ids <- c("0" = "0_activated_memory_CD4",
                     "1" = "1_cytotoxic_CD8",
                     "2" = "2_cytotoxic_CD8",
                     "3" = "3_CD8_DNA_replication",
                     "4" = "4_CD8 heat_shock_proteins",
                     "5" = "5_CD8_proliferation",
                     "6" = "6_CD8_(CD4)_exhausted",
                     "7" = "5_CD8_proliferation",
                     "8" = "5_CD8_proliferation")
names(new.cluster.ids) <- levels(T_cells)
T_cells <- RenameIdents(T_cells, new.cluster.ids)
DimPlot(T_cells, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

#M# saveRDS(T_cells, file = "objects/tils_all_.45_integrgate_annotated_merge_prolif.rds")


# info of proliferation cluster in all cells ----------------------------------
cluster_5_cells <- subset(T_cells, idents = c("5_CD8_proliferation"))

Res_cell_counts <- table(cluster_5_cells$Treatment)["Responder"]
Non_res_cell_counts <- table(cluster_5_cells$Treatment)["Non_Responder"]
total_cells_in_cluster <- Res_cell_counts + Non_res_cell_counts

res_percentage <- (Res_cell_counts / total_cells_in_cluster) * 100
non_res_percentage <-  (Non_res_cell_counts / total_cells_in_cluster) * 100



cluster_5_stats <- data.frame(
  Treatment = names(cell_counts),
  Cluster_5_Cell_Counts = as.vector(cell_counts),
  Total_Cell_Counts = as.vector(total_cells),
  Percentage_in_Cluster_5 = as.vector(percentages)
)

print(cluster_5_stats)


# signature score function ------------
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
               cols= c("#ff7f0e", "#1f77b4")) +  theme_classic(base_size = 14) +
    theme(text = element_text(size=18, colour = "black")) + RotatedAxis() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    labs(title = "", y = name,  x="") + theme(legend.position="right") +
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") +
    annotate("text", x = 1.5, y=annot.max.y, label = label, size = 6) +
    scale_y_continuous(labels = comma)
  return (p)
  
}
# cytokine tool ------------------------------------------------------------
cyto_df <- read_csv("excels/df_irea_clus_0.csv")
cy <- cyto_df %>% select(c("Enrichment Score", "Cytokine")) %>% rename('Cluster_0' = 'Enrichment Score')

cyto_df_clus_1 <- read_csv("excels/df_irea_clus_1.csv") %>% select(c("Enrichment Score", "Cytokine")) %>% rename('Cluster_1' = 'Enrichment Score')

cyto_df_clus_2 <- read_csv("excels/df_irea_clus_2.csv") %>% select(c("Enrichment Score", "Cytokine")) %>% rename('Cluster_2' = 'Enrichment Score')

merged_data <- cy %>%
  full_join(cyto_df_clus_1, by = "Cytokine") %>%
  full_join(cyto_df_clus_2, by = "Cytokine") %>%  
  column_to_rownames(var = "Cytokine")

long_data <- merged_data %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = starts_with("Cluster"), names_to = "Cluster", values_to = "Score")

# write_csv(long_data, "excels/long_cyto_enrich.csv") 
# merged_data_with_genes <- merged_data %>%
 #  rownames_to_column(var = "Gene")
# write_csv(merged_data_with_genes, "excels/cyto_by_clusteR_enrich.csv") 

long_data <- read_csv("excels/long_cyto_enrich.csv")
merged_data <- read_csv("excels/cyto_by_clusteR_enrich.csv")

clus2_top_10 <- merged_data %>% 
  top_n(n=10, wt=Cluster_2) %>% 
  pull(Gene)

clus_2_lowest_10 <- merged_data %>% 
  top_n(n=-10, wt=Cluster_2) %>% 
  pull(Gene)

selected_genes <- c(clus_2_lowest_10, clus2_top_10)

filtered_data <- long_data %>%
  filter(Gene %in% selected_genes)

#M# control gene order in plot : have to put all genes in list
gene_order <- c(clus2_top_10, clus_2_lowest_10)
filtered_data$Gene <- factor(filtered_data$Gene, levels = gene_order)

# Create the bar plot
ggplot(filtered_data, aes(x = Gene, y = Score, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Scores for Selected Genes Across Clusters",
       x = "Gene",
       y = "Score",
       fill = "Cluster") + 
  coord_flip() + # Flip the coordinates
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file = "figures/cytokine tool/clus_2.png", dpi=300, width=5, height=5)

#M# cluster 0 10 highest and lowest scored cytokines
clus0_top_10 <- merged_data %>% 
  top_n(n=10, wt=Cluster_0) %>% 
  pull(Gene)

clus_0_lowest_10 <- merged_data %>% 
  top_n(n=-10, wt=Cluster_0) %>% 
  pull(Gene)

selected_genes <- c(clus_0_lowest_10, clus0_top_10)

filtered_data <- long_data %>%
  filter(Gene %in% selected_genes)

#M# control gene order in plot : have to put all genes in list
gene_order <- c(clus0_top_10, clus_0_lowest_10)
filtered_data$Gene <- factor(filtered_data$Gene, levels = gene_order)

# Create the bar plot
ggplot(filtered_data, aes(x = Gene, y = Score, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Scores for Selected Genes Across Clusters",
       x = "Gene",
       y = "Score",
       fill = "Cluster") + 
  coord_flip() + # Flip the coordinates
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file = "figures/cytokine tool/clus_0.png", dpi=300, width=5, height=5)

#M# cluster 1 10 highest and lowest scored cytokines
clus_1_top_10 <- merged_data %>% 
  top_n(n=10, wt=Cluster_1) %>% 
  pull(Gene)

clus_1_lowest_10 <- merged_data %>% 
  top_n(n=-10, wt=Cluster_1) %>% 
  pull(Gene)

selected_genes <- c(clus_1_lowest_10, clus_1_top_10)

filtered_data <- long_data %>%
  filter(Gene %in% selected_genes)

#M# control gene order in plot : have to put all genes in list
gene_order <- c(clus_1_top_10, clus_1_lowest_10)
filtered_data$Gene <- factor(filtered_data$Gene, levels = gene_order)

# Create the bar plot
ggplot(filtered_data, aes(x = Gene, y = Score, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Scores for Selected Genes Across Clusters",
       x = "Gene",
       y = "Score",
       fill = "Cluster") + 
  coord_flip() + # Flip the coordinates
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file = "figures/cytokine tool/clus_1.png", dpi=300, width=5, height=5)


# conversion --------------------------------------------------------------
#M# convert genes from human(uppercase) to mouse format for cytokyne tool

# Define the function
capitalize_genes <- function(gene_string) {
  gene_list <- strsplit(gene_string, "\n")[[1]]
  gene_list <- trimws(gene_list)
  capitalized_genes <- sapply(gene_list, function(gene) {
    paste0(toupper(substr(gene, 1, 1)), tolower(substr(gene, 2, nchar(gene))))
  })
  capitalized_gene_string <- paste(capitalized_genes, collapse = "\n")
  return(capitalized_gene_string)
}

genes <- "MSRB2
CFAP36
SVIP
MALSU1
PFDN4
ARCN1
ZNHIT3
GLRX5
RPS26
MRPS33
RPLP0
RPLP1
MED28
RTRAF
PSMD4
EEF1B2
RPS10
NDUFB4
PPA1
RPL23
RPL36A
GAS5
PRDX1
SYF2
ZNF706
RPS17
EEF1G
RPL37A
ADA
ATP5MPL
RPL13A
RPL37
RPS4X
LAMTOR5
LDHB
RPL4
UQCRB
RPS20
GNG10
ABRACL
RPS18
MRPL20
COX14
RPS6
RSL24D1
RPS27L
NOSIP
RPS8
RPL31
NDUFS5
"

capitalized_genes <- capitalize_genes(genes)
cat(capitalized_genes)


# make list of genes into vector ------------------------------------------
gene_list <- "TYROBP
FCER1G
NCR2
KIT
TRDC
RHOBTB3
ADAM12
CXXC5
RASSF8
KRT81
SYK
FAM49A
ZBTB16
TLE1
AREG
TRIO
KRT86
TMIGD2
AL136456.1
GSN
PLD1
TNFRSF18
NCR1
APOL4
SH2D1B
KLRC1
IKZF2
PDE6G
DOCK5
FUCA1
FES
CST3
LINC00299
NCAM1
TSPAN4
IL4I1
FSD1
RNF165
ITM2C
CHKA
AHDC1
SLFN12
PLEK
CEBPD
HIP1R
ZNF169
SKAP2
DPF3
ITGAX
GALE
"
gene_list <- strsplit(gene_list, "\n")[[1]]
gene_list <- trimws(gene_list)
gene_list[2]
a <- gene_list
b <- gene_list
a %in% b

# figures for presentation ------------------------------------------------
#M# all cells 
# exhaustion
VlnPlot(
  T_cells, 
  features = c("HAVCR2", "LAG3", "TIGIT", "TOX"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
    ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(filename = "vln_activation_exhastion_cluster_by_treatment_1.png" , path = "figures/Tcells/", dpi=300, width=12, height=10)


VlnPlot(
  T_cells,
  features = c("HAVCR2","LAG3","TIGIT","TOX"),
  assay = "RNA",
  stack=TRUE,
  flip= TRUE,
  group.by = "Treatment",
  fill.by ="ident"
  )+ 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/Tcells/vln_activation_exhastion_by_treatment_1.png", dpi=300, width=6, height=10)



# effector markers
VlnPlot(
  T_cells, 
  features = c("GZMB","PRF1","IFNG","TNF", "ICOS"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/Tcells/Vln_T_effector_cluster_by_treatment_1.png", dpi=300, width=12, height=10, limitsize=FALSE)


VlnPlot(
  T_cells,
  features = c("GZMB","PRF1","IFNG","TNF", "ICOS"),
  assay = "RNA",
  stack=TRUE,
  flip= TRUE,
  group.by = "Treatment",
  fill.by ="ident"
)+ 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/Tcells/Vln_T_effector_by_treatment_1.png", dpi=300, width=6, height=10, limitsize=FALSE)


# naive markers
VlnPlot(
  T_cells, 
  features = c("CD69","ENTPD1"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/Tcells/Vln_T_naive_cluster_by_treatment_1.png", dpi=300, width=12, height=10, limitsize=FALSE)


VlnPlot(
  T_cells,
  features = c("CD69","ENTPD1"),
  assay = "RNA",
  stack=TRUE,
  flip= TRUE,
  group.by = "Treatment",
  fill.by ="ident"
)+ 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/Tcells/Vln_T_naive_by_treatment_1.png", dpi=300, width=6, height=10, limitsize=FALSE)


#_____________________________________________________________________________________________________________________
#M# cluster 0
#M# subset for cluster 0 only for vln plots : to avoid non-informative cluster plotting
cd4_cluster0 <- subset(cd4_cells, idents = "0_Tefm")
DefaultAssay(cd4_cluster0) <- "RNA"
# exhaustion
VlnPlot(
  cd4_cluster0, 
  features = c("HAVCR2", "LAG3", "TIGIT", "TOX"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(filename = "vln_activation_exhastion_cluster_by_treatment_1.png" , path = "figures/cd4_cluster0/", dpi=300, width=8, height=10)


VlnPlot(
  cd4_cluster0,
  features = c("HAVCR2","LAG3","TIGIT","TOX"),
  assay = "RNA",
  stack=TRUE,
  flip= TRUE,
  group.by = "Treatment",
  fill.by ="ident"
)+ 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/cd4_cluster0/vln_activation_exhastion_by_treatment_1.png", dpi=300, width=6, height=10)



# effector markers
VlnPlot(
  cd4_cluster0, 
  features = c("GZMB","PRF1","IFNG","TNF", "ICOS"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/cd4_cluster0/Vln_T_effector_cluster_by_treatment_1.png", dpi=300, width=8, height=10, limitsize=FALSE)


VlnPlot(
  cd4_cluster0,
  features = c("GZMB","PRF1","IFNG","TNF", "ICOS"),
  assay = "RNA",
  stack=TRUE,
  flip= TRUE,
  group.by = "Treatment",
  fill.by ="ident"
)+ 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/cd4_cluster0/Vln_T_effector_by_treatment_1.png", dpi=300, width=6, height=10, limitsize=FALSE)


# naive markers
VlnPlot(
  cd4_cluster0, 
  features = c("CD69","ENTPD1"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/cd4_cells/cluster0_only/Vln_T_naive_cluster_by_treatment_1.png", dpi=300, width=8, height=10, limitsize=FALSE)


VlnPlot(
  cd4_cluster0,
  features = c("CD69","ENTPD1"),
  assay = "RNA",
  stack=TRUE,
  flip= TRUE,
  group.by = "Treatment",
  fill.by ="ident"
)+ 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/cd4_cluster0/Vln_T_naive_by_treatment_1.png", dpi=300, width=6, height=10, limitsize=FALSE)

#_____________________________________________________________________________________________________________________

#M# cluster 1
# exhaustion
VlnPlot(
  clus1, 
  features = c("HAVCR2", "LAG3", "TIGIT", "TOX"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(filename = "vln_activation_exhastion_cluster_by_treatment_1.png" , path = "figures/clus1/", dpi=300, width=12, height=10)


VlnPlot(
  clus1,
  features = c("HAVCR2","LAG3","TIGIT","TOX"),
  assay = "RNA",
  stack=TRUE,
  flip= TRUE,
  group.by = "Treatment",
  fill.by ="ident"
)+ 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus1/vln_activation_exhastion_by_treatment_1.png", dpi=300, width=6, height=10)



# effector markers
VlnPlot(
  clus1, 
  features = c("GZMB","PRF1","IFNG","TNF", "ICOS"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus1/Vln_T_effector_cluster_by_treatment_1.png", dpi=300, width=12, height=10, limitsize=FALSE)


VlnPlot(
  clus1,
  features = c("GZMB","PRF1","IFNG","TNF", "ICOS"),
  assay = "RNA",
  stack=TRUE,
  flip= TRUE,
  group.by = "Treatment",
  fill.by ="ident"
)+ 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus1/Vln_T_effector_by_treatment_1.png", dpi=300, width=6, height=10, limitsize=FALSE)


# naive markers
VlnPlot(
  clus1, 
  features = c("CD69","ENTPD1"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus1/Vln_T_naive_cluster_by_treatment_1.png", dpi=300, width=12, height=10, limitsize=FALSE)


VlnPlot(
  clus1,
  features = c("CD69","ENTPD1"),
  assay = "RNA",
  stack=TRUE,
  flip= TRUE,
  group.by = "Treatment",
  fill.by ="ident"
)+ 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus1/Vln_T_naive_by_treatment_1.png", dpi=300, width=6, height=10, limitsize=FALSE)


#_____________________________________________________________________________________________________________________

#M# cluster 2
# exhaustion
VlnPlot(
  clus2, 
  features = c("HAVCR2", "TIGIT"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(filename = "vln_activation_exhastion_cluster_by_treatment_1.png" , path = "figures/clus2/", dpi=300, width=12, height=10)


VlnPlot(
  clus2,
  features = c("HAVCR2", "TIGIT"),
  assay = "RNA",
  stack=TRUE,
  flip= TRUE,
  group.by = "Treatment",
  fill.by ="ident"
)+ 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus2/vln_activation_exhastion_by_treatment_1.png", dpi=300, width=6, height=10)



# effector markers
VlnPlot(
  clus2, 
  features = c("GZMB","PRF1","IFNG","TNF", "ICOS"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus2/Vln_T_effector_cluster_by_treatment_1.png", dpi=300, width=12, height=10, limitsize=FALSE)


VlnPlot(
  clus2,
  features = c("GZMB","PRF1","IFNG","TNF", "ICOS"),
  assay = "RNA",
  stack=TRUE,
  flip= TRUE,
  group.by = "Treatment",
  fill.by ="ident"
)+ 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus2/Vln_T_effector_by_treatment_1.png", dpi=300, width=6, height=10, limitsize=FALSE)

# naive markers
VlnPlot(
  clus2, 
  features = c("CD69","ENTPD1"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus2/Vln_T_naive_cluster_by_treatment_1.png", dpi=300, width=12, height=10, limitsize=FALSE)


VlnPlot(
  clus2,
  features = c("CD69","ENTPD1"),
  assay = "RNA",
  stack=TRUE,
  flip= TRUE,
  group.by = "Treatment",
  fill.by ="ident"
)+ 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus2/Vln_T_naive_by_treatment_1.png", dpi=300, width=6, height=10, limitsize=FALSE)


VlnPlot(
  clus2, 
  features = c("ISG20", "TGFBR3"), 
  assay = "RNA", 
  stack = TRUE, 
  flip = TRUE, 
  split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus2/de_gene_search/interesting_DE_genes_for_presentation.png", dpi=300, width=16, height=12)

# make tsne's by responders/non responders ---------------
# all cells
responders <- subset(T_cells, subset = Treatment == "Responder")
non_responders <- subset(T_cells, subset = Treatment == "Non_Responder")
library(cowplot)
p1 <- DimPlot(responders, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) +
  ggtitle("Responders")
p2 <- DimPlot(non_responders, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) +
  ggtitle("Non Responders")
plot_grid(p1, p2, ncol = 2)
ggsave(file = "figures/Tcells/tsne_by_treatment_together.png", dpi=300, width=18, height=8)

# cluster 0 
responders <- subset(cd4_cells, subset = Treatment == "Responder")
non_responders <- subset(cd4_cells, subset = Treatment == "Non_Responder")
p1 <- DimPlot(responders, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) +
  ggtitle("Responders")
p2 <- DimPlot(non_responders, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) +
  ggtitle("Non Responders")
plot_grid(p1, p2, ncol = 2)
ggsave(file = "figures/cd4_cells/tsne_by_treatment_together.png", dpi=300, width=18, height=8)

# cluster 1
responders <- subset(clus1, subset = Treatment == "Responder")
non_responders <- subset(clus1, subset = Treatment == "Non_Responder")
p1 <- DimPlot(responders, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) +
  ggtitle("Responders")
p2 <- DimPlot(non_responders, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) +
  ggtitle("Non Responders")
plot_grid(p1, p2, ncol = 2)
ggsave(file = "figures/clus1/tsne_by_treatment_together.png", dpi=300, width=14, height=5)

# cluster 2
responders <- subset(clus2, subset = Treatment == "Responder")
non_responders <- subset(clus2, subset = Treatment == "Non_Responder")
p1 <- DimPlot(responders, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) +
  ggtitle("Responders")
p2 <- DimPlot(non_responders, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) +
  ggtitle("Non Responders")
plot_grid(p1, p2, ncol = 2)
ggsave(file = "figures/clus2/tsne_by_treatment_together.png", dpi=300, width=14, height=5)




































# calculate cd8/cd4 ratio in responders and non ---------

# add cell type for ratio calculation
T_cells$cell_type <- ifelse(T_cells@assays$RNA@data["CD8A", ] > 0 | T_cells@assays$RNA@data["CD8B", ] > 0, 
                            "CD8", 
                            ifelse(T_cells@assays$RNA@data["CD4", ] > 0, 
                                   "CD4", 
                                   "Other"))
table(T_cells$cell_type)

responders <- subset(T_cells, subset = Treatment == "Responder")
non_responders <- subset(T_cells, subset = Treatment == "Non_Responder")
cd8_count_responders <- sum(responders$cell_type == "CD8")
cd4_count_responders <- sum(responders$cell_type == "CD4")
cd8_count_non_responders <- sum(non_responders$cell_type == "CD8")
cd4_count_non_responders <- sum(non_responders$cell_type == "CD4")
cd8_cd4_ratio_responders <- cd8_count_responders / cd4_count_responders
cd8_cd4_ratio_non_responders <- cd8_count_non_responders / cd4_count_non_responders
cat("CD8/CD4 ratio in responders:", cd8_cd4_ratio_responders, "\n")
cat("CD8/CD4 ratio in non-responders:", cd8_cd4_ratio_non_responders, "\n")

DefaultAssay(T_cells) <- "RNA"
#M# calculate ratio per patient
df <- T_cells@meta.data %>%
  as.data.frame() %>% 
  group_by(Sample, cell_type) %>%
  summarize(Freq = n())
#M# why doesnt this work?! 11.9.24


df <-  T_cells@meta.data %>%
  as.data.frame() %>% 
  group_by(Sample, cell_type) %>%
  tally(name = "Freq") %>% 
  as.data.frame()



calculate_ratio <- function(df, sample_id) {
  cd4_count <- df[df$Sample == sample_id & df$cell_type == "CD4", "Freq"]
  cd8_count <- df[df$Sample == sample_id & df$cell_type == "CD8", "Freq"]
  if (cd4_count > 0) {
    ratio <- cd8_count / cd4_count
  } else {
    ratio <- NA
  }
  return(ratio = ratio)
}
patients <- unique(T_cells$Sample)
treatments <- sapply(patients, function(p) unique(subset(T_cells, Sample == p)$Treatment))
ratios_list <- list()
for (patient in patients) {
  ratios_list[[patient]] <- calculate_ratio(df, patient)  
}
ratios_df <- do.call(rbind, lapply(ratios_list, as.data.frame))

ratios_df <- data.frame(Sample = patients, Treatment = treatments, ratios_df)
ratios_df$ratio <- as.numeric(ratios_df$ratio)
print(ratios_df)

colnames(ratios_df)[colnames(ratios_df) == "X..i.."] <- "ratio"

summary_stats_responders <- ratios_df %>%
  filter(Treatment == "Responder") %>%
  summarize(
    Mean_Ratio = mean(ratio, na.rm = TRUE),
    SD_Ratio = sd(ratio, na.rm = TRUE)
  )
print(summary_stats_responders)
summary_stats_non_responders <- ratios_df %>%
  filter(Treatment == "Non_Responder") %>%
  summarize(
    Mean_Ratio = mean(ratio, na.rm = TRUE),
    SD_Ratio = sd(ratio, na.rm = TRUE)
  )
print(summary_stats_non_responders)
summary_stats <- data.frame(
  Treatment = c("Responder", "Non_Responder"),
  Mean_Ratio = c(summary_stats_responders$Mean_Ratio, summary_stats_non_responders$Mean_Ratio),
  SD_Ratio = c(summary_stats_responders$SD_Ratio, summary_stats_non_responders$SD_Ratio)
)
ggplot(summary_stats, aes(x = Treatment, y = Mean_Ratio, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) +
  geom_errorbar(aes(ymin = Mean_Ratio - SD_Ratio, ymax = Mean_Ratio + SD_Ratio),
                width = 0.2, position = position_dodge(0.6)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 40),
    axis.text.y = element_text(size = 12)
  )
ggsave(file = "figures/Tcells/cd8_cd4_ratio_2.png", dpi=300, width=4, height=6)


















# cd40lg zoom in --------
# all cells
VlnPlot(T_cells, features = c("CD40LG"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment",
        pt.size = 0
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/Tcells/cd40lg.png", dpi=300, width=10, height=6)


# cd4 cluster
VlnPlot(cd4_cells, features = c("CD40LG"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment",
        pt.size = 0
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/cd4_cells/cd40lg_cluster_0.png", dpi=300, width=10, height=6)

# clus 1
VlnPlot(clus1, features = c("CD40LG"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment",
        pt.size = 0
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus1/cd40lg_cluster_1.png", dpi=300, width=10, height=6)

# clus 2
VlnPlot(clus2, features = c("CD40LG"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment",
        pt.size = 0
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus2/cd40lg_cluster_2.png", dpi=300, width=10, height=6)

# 41bb zoom in----------
VlnPlot(cd4_cells, features = c("TNFRSF9"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/cd4_cells/41bb_cluster_0.png", dpi=300, width=10, height=6)


VlnPlot(clus1, features = c("TNFRSF9"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus1/41bb_cluster_1.png", dpi=300, width=10, height=6)


VlnPlot(clus2, features = c("TNFRSF9"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus2/41bb_cluster_2.png", dpi=300, width=10, height=6)

# IL21 check-------------

VlnPlot(cd4_cells, features = c("IL21"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/cd4_cells/IL21_cluster_0.png", dpi=300, width=10, height=6)

VlnPlot(clus1, features = c("IL21"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus1/IL21_cluster_1.png", dpi=300, width=10, height=6)

VlnPlot(clus2, features = c("IL21"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus2/IL21_cluster_2.png", dpi=300, width=10, height=6)

# IL1 check----------

VlnPlot(cd4_cells, features = c("IL1A"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/cd4_cells/IL1_cluster_0.png", dpi=300, width=10, height=6)

VlnPlot(clus1, features = c("IL1A"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus1/IL1_cluster_1.png", dpi=300, width=10, height=6)

VlnPlot(clus2, features = c("IL1A"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus2/IL1_cluster_2.png", dpi=300, width=10, height=6)
# checking ------------

VlnPlot(T_cells, features = c("KLF2", "SLAMF6", "CD27", "IL4"),
        assay = "RNA", 
        stack = TRUE, 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/Tcells/check_paper_DN.png", dpi=300, width=10, height=6)

VlnPlot(clus1, features = c("TGFBR3", "HIVEP1"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus1/1.png", dpi=300, width=10, height=6)

VlnPlot(clus2, features = c("TGFBR3", "HIVEP1"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus2/2.png", dpi=300, width=10, height=6)

















# immune dictionary rds ---------

immune_dict = readRDS("D:/user/Downloads/ref_data_T_cell_CD8.RDS")

immune_dict$
unique(immune_dict$celltype)

DimPlot(immune_dict , reduction = "tsne", label= T, pt.size=0.5, label.size = 10)



# filter genes from excel--------------
#M# cd4
df = read.csv("D:/Michael/git_check/tils/excels/DE_genes_res_nonRes_CD4.csv")
filtered_df <- df[df$gene %in% gene_list, ]
filtered_df <- filtered_df[filtered_df$p_val_adj <= 0.0011, ]
write_csv(filtered_df, "excels/filtered_DE_genes_res_nonRes_CD4.csv")

#M# cd8 late (clus 1)
df = read.csv("D:/Michael/git_check/tils/excels/DE_genes_res_nonRes_Clus1.csv")
filtered_df <- df[df$gene %in% gene_list, ]
filtered_df <- filtered_df[filtered_df$p_val_adj <= 0.0011, ]
write_csv(filtered_df, "excels/filtered_DE_genes_res_nonRes_Clus1.csv")

#M# cd8 early (clus 2)
df = read.csv("D:/Michael/git_check/tils/excels/DE_genes_res_nonRes_Clus2.csv")
filtered_df <- df[df$gene %in% gene_list, ]
filtered_df <- filtered_df[filtered_df$p_val_adj <= 0.0011, ]
write_csv(filtered_df, "excels/filtered_DE_genes_res_nonRes_Clus2.csv")



# genesets --------------
# install.packages("BiocManager")
# BiocManager::install("fgsea")

il2_biocarta <- gmtPathways("geneSets/BIOCARTA_IL2_PATHWAY.v2023.2.Hs.gmt")[[1]]

# search of interesting DE genes between responders & non-responders in cd4 cluster---------
cd4_cells = readRDS("objects/tils_0_.35.rds")
filtered_DE_genes_res_nonRes_CD4 <- read.csv("excels/filtered_DE_genes_res_nonRes_CD4.csv")

genes <- filtered_DE_genes_res_nonRes_CD4 %>% 
  pull(gene)

#M# one-by-one gene violin plot
for (gene in genes) {
  VlnPlot(cd4_cells, features = c(gene),
          assay = "RNA", 
          flip = TRUE, 
          split.by = "Treatment"
  ) + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 24, face = "italic"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.ticks.y = element_line(size = 0.5),
      strip.text.y = element_text(angle = 0, size = 16, face = "bold")
    ) + 
    geom_boxplot(alpha = 0.3, show.legend = FALSE)
  ggsave(file = paste0("figures/cd4_cells/de_gene_search/",gene, ".png"), dpi=300, width=10, height=6)
  
}
#M# or same but more efficient for large scale picture
i <- 1
while(i<length(genes)){
  VlnPlot(
    cd4_cells, 
    features = c(genes[i], genes[i+1], genes[i+2], genes[i+3], genes[i+4], genes[i+5], genes[i+6], genes[i+7]), 
    assay = "RNA", 
    stack = TRUE, 
    flip = TRUE, 
    split.by = "Treatment"
  ) + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 24, face = "italic"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.ticks.y = element_line(size = 0.5),
      strip.text.y = element_text(angle = 0, size = 16, face = "bold")
    ) +
    geom_boxplot(alpha = 0.3, show.legend = FALSE)
  ggsave(file = paste0("figures/cd4_cells/de_gene_search/","run_", i,  ".png"), dpi=300, width=16, height=12)
  i <- i+8
}
i <- 1




# search of interesting DE genes between responders & non-responders in cluster 1---------
clus1 = readRDS("objects/clus1_no_cd4_second_sub_14_0.5.rds")
filtered_DE_genes_res_nonRes_Clus1 <- read.csv("excels/filtered_DE_genes_res_nonRes_Clus1.csv")

genes <- filtered_DE_genes_res_nonRes_Clus1 %>% 
  pull(gene)

#M# or same but more efficient for large scale picture
i <- 1
while(i<length(genes)){
  VlnPlot(
    clus1, 
    features = c(genes[i], genes[i+1], genes[i+2], genes[i+3], genes[i+4], genes[i+5], genes[i+6], genes[i+7]), 
    assay = "RNA", 
    stack = TRUE, 
    flip = TRUE, 
    split.by = "Treatment"
  ) + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 24, face = "italic"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.ticks.y = element_line(size = 0.5),
      strip.text.y = element_text(angle = 0, size = 16, face = "bold")
    ) +
    geom_boxplot(alpha = 0.3, show.legend = FALSE)
  ggsave(file = paste0("figures/clus1/de_gene_search/","interesting_DE_genes_run_", i,  ".png"), dpi=300, width=16, height=12)
  i <- i+8
}
i <- 1


# #M# one-by-one gene violin plot
# for (gene in genes) {
#   VlnPlot(clus1, features = c(gene),
#           assay = "RNA", 
#           flip = TRUE, 
#           split.by = "Treatment"
#   ) + 
#     theme_classic() +
#     theme(
#       axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
#       axis.title.x = element_blank(),
#       axis.text.y = element_text(size = 24, face = "italic"),
#       axis.title.y = element_text(size = 20, face = "bold"),
#       axis.ticks.y = element_line(size = 0.5),
#       strip.text.y = element_text(angle = 0, size = 16, face = "bold")
#     ) + 
#     geom_boxplot(alpha = 0.3, show.legend = FALSE)
#   ggsave(file = paste0("figures/clus1/de_gene_search/",gene, ".png"), dpi=300, width=10, height=6)
#   
# }

# search of interesting DE genes between responders & non-responders in cluster 2---------
clus2 = readRDS("objects/tils_2.rds")
filtered_DE_genes_res_nonRes_Clus2 <- read.csv("excels/filtered_DE_genes_res_nonRes_Clus2.csv")

genes <- filtered_DE_genes_res_nonRes_Clus2 %>% 
  pull(gene)

#M# or same but more efficient for large scale picture
i <- 1
while(i<length(genes)){
  VlnPlot(
    clus2, 
    features = c(genes[i], genes[i+1], genes[i+2], genes[i+3], genes[i+4], genes[i+5], genes[i+6], genes[i+7]), 
    assay = "RNA", 
    stack = TRUE, 
    flip = TRUE, 
    split.by = "Treatment"
  ) + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 24, face = "italic"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.ticks.y = element_line(size = 0.5),
      strip.text.y = element_text(angle = 0, size = 16, face = "bold")
    ) +
    geom_boxplot(alpha = 0.3, show.legend = FALSE)
  ggsave(file = paste0("figures/clus2/de_gene_search/","interesting_DE_genes_run_", i,  ".png"), dpi=300, width=16, height=12)
  i <- i+8
}
i <- 1

# 
# #M# one-by-one gene violin plot
# for (gene in genes) {
#   VlnPlot(clus2, features = c(gene),
#           assay = "RNA", 
#           flip = TRUE, 
#           split.by = "Treatment"
#   ) + 
#     theme_classic() +
#     theme(
#       axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
#       axis.title.x = element_blank(),
#       axis.text.y = element_text(size = 24, face = "italic"),
#       axis.title.y = element_text(size = 20, face = "bold"),
#       axis.ticks.y = element_line(size = 0.5),
#       strip.text.y = element_text(angle = 0, size = 16, face = "bold")
#     ) + 
#     geom_boxplot(alpha = 0.3, show.legend = FALSE)
#   ggsave(file = paste0("figures/clus2/de_gene_search/",gene, ".png"), dpi=300, width=10, height=6)
#   
# }









# cluster 0 of cd4 cluster (Tefm) annotation -------------------

cd4_cells = readRDS("objects/tils_0_.35.rds")
x1=4.6

Tefm_markers <- c("GPR183", "AQP3", "TNFRSF4", "LTB", "KLF2", "RASA3")

DefaultAssay(cd4_cells) <- "RNA"

FeaturePlot(cd4_cells, features = Tefm_markers, label.size = 8, pt.size = 1, label=F, ncol = 3, order = T,reduction = "tsne") + 
  scale_color_gradient(low = "grey", high = "blue")
ggsave(file = "figures/cd4_cells/Tefm_markers.png", dpi=300, width=14, height=x1*2)


# FeaturePlot(cd4_cells, features = "RASA3", label.size = 8, pt.size = 0.2, label=F, order = T,reduction = "tsne") + 
#   scale_color_gradientn(colors = c("grey", "blue"), 
#                         limits = c(0, 3),
#                         na.value = "grey")
# ggsave(file = "figures/cd4_cells/rasa.png", dpi=300, width=4, height=3)
# 
# 
FeaturePlot(cd4_cells, features = "LTB", label.size = 8, pt.size = 1, label=F, reduction = "tsne") +
  scale_color_gradient(low = "grey", high = "blue")
ggsave(file = "figures/cd4_cells/LTB.png", dpi=300, width=4, height=3)



# TNF family research cd4 cluster -------------

tnf_family <- read_csv("excels/TNF_family_18.8.24.csv") %>% 
  pull("Gene")

cd4_cells = readRDS("objects/tils_0_.35.rds")
#M# cd4 cluster : 
i <- 1
while(i<length(tnf_family)){
  VlnPlot(
    cd4_cells, 
    features = c(tnf_family[i], tnf_family[i+1], tnf_family[i+2], tnf_family[i+3], tnf_family[i+4], tnf_family[i+5], tnf_family[i+6], tnf_family[i+7]), 
    assay = "RNA", 
    stack = TRUE, 
    flip = TRUE, 
    split.by = "Treatment"
  ) + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 24, face = "italic"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.ticks.y = element_line(size = 0.5),
      strip.text.y = element_text(angle = 0, size = 16, face = "bold")
    ) +
    geom_boxplot(alpha = 0.3, show.legend = FALSE)
  ggsave(file = paste0("figures/cd4_cells/de_gene_search/","tnf_search_run_", i,  ".png"), dpi=300, width=16, height=12)
  i <- i+8
}
i <- 1


# TNF family research cluster 1-------------
tnf_family <- read_csv("excels/TNF_family_18.8.24.csv") %>% 
  pull("Gene")

clus1 = readRDS("objects/clus1_no_cd4_14_0.5.rds")

i <- 1
while(i<length(tnf_family)){
  VlnPlot(
    clus1, 
    features = c(tnf_family[i], tnf_family[i+1], tnf_family[i+2], tnf_family[i+3], tnf_family[i+4], tnf_family[i+5], tnf_family[i+6], tnf_family[i+7]), 
    assay = "RNA", 
    stack = TRUE, 
    flip = TRUE, 
    split.by = "Treatment"
  ) + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 24, face = "italic"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.ticks.y = element_line(size = 0.5),
      strip.text.y = element_text(angle = 0, size = 16, face = "bold")
    ) +
    geom_boxplot(alpha = 0.3, show.legend = FALSE)
  ggsave(file = paste0("figures/clus1/de_gene_search/","tnf_search_run_", i,  ".png"), dpi=300, width=16, height=12)
  i <- i+8
}
i <- 1




# TNF family research cluster 2-------------
tnf_family <- read_csv("excels/TNF_family_18.8.24.csv") %>% 
  pull("Gene")

clus2 = readRDS("objects/tils_2.rds")

i <- 1
while(i<length(tnf_family)){
  VlnPlot(
    clus2, 
    features = c(tnf_family[i], tnf_family[i+1], tnf_family[i+2], tnf_family[i+3], tnf_family[i+4], tnf_family[i+5], tnf_family[i+6], tnf_family[i+7]), 
    assay = "RNA", 
    stack = TRUE, 
    flip = TRUE, 
    split.by = "Treatment"
  ) + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 24, face = "italic"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.ticks.y = element_line(size = 0.5),
      strip.text.y = element_text(angle = 0, size = 16, face = "bold")
    ) +
    geom_boxplot(alpha = 0.3, show.legend = FALSE)
  ggsave(file = paste0("figures/clus2/de_gene_search/","tnf_search_run_", i,  ".png"), dpi=300, width=16, height=12)
  i <- i+8
}
i <- 1
















# cd40 visualisation : (cd40lg receptor)--------
VlnPlot(cd4_cells, features = c("CD40"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/cd4_cells/CD40_0.png", dpi=300, width=10, height=6)

VlnPlot(clus1, features = c("CD40"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus1/CD40_1.png", dpi=300, width=10, height=6)

VlnPlot(clus2, features = c("CD40"),
        assay = "RNA", 
        flip = TRUE, 
        split.by = "Treatment"
) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) + 
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/clus2/CD40_2.png", dpi=300, width=10, height=6)














# feature plots for annotation all cells for presentation ----------

T_cells = readRDS("objects/tils_all_.45_integrgate_annotated (2).rds")
DefaultAssay(T_cells) <- "RNA"

#cluster 0
Activated_memory_CD4 <- c("IL7R", "TNFRSF4", "AQP3", "GPR183", "CD4", "LTB")
FeaturePlot(T_cells, features = Activated_memory_CD4, label.size = 8, pt.size = 1, label=F, ncol = 3, order = T,reduction = "tsne")
ggsave(file = "figures/Tcells/Activated_memory_CD4_Markers.png", dpi=300, width=20, height=x1*2)

#cluster 1
Cytotoxic_CD8 <- c("NKG7", "GZMB", "GNLY", "PRF1", "CD8A", "JUNB")
FeaturePlot(T_cells, features = Cytotoxic_CD8, label.size = 8, pt.size = 1, label=F, ncol = 3, order = T,reduction = "tsne")
ggsave(file = "figures/Tcells/Cytotoxic_CD8_Markers.png", dpi=300, width=20, height=x1*2)

#cluster 2
Cytotoxic_CD8_early <- c("CD27", "CCL5")
FeaturePlot(T_cells, features = Cytotoxic_CD8_early, label.size = 8, pt.size = 1, label=F, ncol = 2, order = T,reduction = "tsne")
ggsave(file = "figures/Tcells/Cytotoxic_CD8_early_Markers.png", dpi=300, width=10, height=x1)

#cluster 3
CD8_DNA_replication <- c("GINS2", "MCM7", "CHEK1")
FeaturePlot(T_cells, features = CD8_DNA_replication, label.size = 8, pt.size = 1, label=F, ncol = 3, order = T,reduction = "tsne")
ggsave(file = "figures/Tcells/CD8_DNA_replication_Markers.png", dpi=300, width=20, height=x1)

#cluster 4 
CD8_Heat_Shock_Proteins <- c("HSPA6", "HSPA1A", "CHORDC1", "HSPD1")
FeaturePlot(T_cells, features = CD8_Heat_Shock_Proteins, label.size = 8, pt.size = 1, label=F, ncol = 3, order = T,reduction = "tsne")
ggsave(file = "figures/Tcells/CD8_Heat_Shock_Proteins_Markers.png", dpi=300, width=20, height=x1*3)

#cluster 5 HIST1H1B, HIST1H4C
CD8_Cell_cycle <- c("HIST1H1B", "HIST1H4C")
FeaturePlot(T_cells, features = CD8_Cell_cycle, label.size = 8, pt.size = 1, label=F, ncol = 2, order = T,reduction = "tsne")
ggsave(file = "figures/Tcells/CD8_Cell_cycle_Markers.png", dpi=300, width=10, height=x1)

#cluster 6
CD8_CD4_exhausted <- c("MIAT", "SETBP1")
FeaturePlot(T_cells, features = CD8_CD4_exhausted, label.size = 8, pt.size = 1, label=F, ncol = 2, order = T,reduction = "tsne")
ggsave(file = "figures/Tcells/CD8_CD4_exhausted_Markers.png", dpi=300, width=10, height=x1)

#cluster 7
CD8_Proliferating_7 <- c("CCNB1", "CCNB2", "TUBB4B", "TOP2A", "MKI67")
FeaturePlot(T_cells, features = CD8_Proliferating_7, label.size = 8, pt.size = 1, label=F, ncol = 3, order = T,reduction = "tsne")
ggsave(file = "figures/Tcells/CD8_Proliferating_7_Markers.png", dpi=300, width=20, height=x1*2)

#cluster 8
CD8_Proliferating_cluster_8 <- c("TOP2A", "MKI67")
FeaturePlot(T_cells, features = CD8_Proliferating_cluster_8, label.size = 8, pt.size = 0.5, label=F, ncol = 2, order = T,reduction = "tsne")
ggsave(file = "figures/Tcells/CD8_Proliferating_cluster_8_Markers.png", dpi=300, width=10, height=x1)





















#proliferation signature between response groups------
prolif = readRDS("objects/prolif.rds")
T_cells = readRDS("objects/tils_all_.45_integrgate_annotated_merge_prolif.rds")

proliferation_sig <- c("MKI67", "STMN1","TOP2A","Nolc1", "Npm1","Ccne2")
proliferation_sig <- toupper(proliferation_sig)
proliferation_sig <- list(proliferation_sig)

# proliferation object (clusters 5, 7, 8)
prolif = AddModuleScore(object = prolif, features = proliferation_sig, name = "proliferation", assay = "RNA")
a <- FeaturePlot(object = prolif, features = "proliferation1", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "proliferation_sig", subtitle = "MKI67, STMN1, TOP2A, Nolc1, Npm1,Ccne2..")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/prolif/proliferation_signature.png.png", dpi=300, width=5, height=5)

gglist <-  list()
name <- "proliferation1"
DefaultAssay(prolif) <- "RNA" 
for(clus in levels(prolif$seurat_clusters)){
  print(clus)
  obj <- subset(prolif, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'proliferation1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/prolif/proliferation1_vln_signature.png", dpi=300, width=16, height=10)

# all cells
T_cells = AddModuleScore(object = T_cells, features = proliferation_sig, name = "proliferation", assay = "RNA")
a <- FeaturePlot(object = T_cells, features = "proliferation1", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "proliferation_sig", subtitle = "MKI67, STMN1, TOP2A, Nolc1, Npm1,Ccne2..")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/Tcells/proliferation_signature.png", dpi=300, width=5, height=5)

gglist <-  list()
name <- "proliferation1"
DefaultAssay(T_cells) <- "RNA" 
for(clus in levels(T_cells$seurat_clusters)){
  print(clus)
  obj <- subset(T_cells, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'proliferation1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/Tcells/proliferation1_vln_signature.png", dpi=300, width=16, height=10)

# responder vs. non responder in all cells all together
p <- SignatureScore(T_cells, 'proliferation1') + ggtitle("proliferation sig")
ggsave(file = "figures/Tcells/proliferation1_all_cells_combined_signature.png", dpi=300, width=16, height=10)


























