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

setwd("D:/Michael/git_check/tils")

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

genes <- "AC243829.4
NELL2
GZMH
CD27
GZMK
TXK
TRIM22
KLRK1
VAV3
AF165147.1
SLAMF7
PECAM1
GZMB
DDX3Y
YBX3
IKZF3
CCL4
AOAH
AL158071.1
NR4A2
CCL5
LINC02446
CRIM1
AC015849.1
ZNF331
NKG7
EPAS1
LYST
CD8B
CCL3
CBLB
SYNE1
LINC00861
SCML4
BPGM
RNF213-AS1
RASA3
TNFSF8
CD8A
GIMAP4
XCL2
MVB12B
IVNS1ABP
FYN
TIGIT
METRNL
AC011476.3
YES1
CD226
PRF1
"

capitalized_genes <- capitalize_genes(genes)
cat(capitalized_genes)


# make list of genes into vector ------------------------------------------
gene_list <- "CD3G
CLN8
CRIP1
CXCR4
GAPDH
GNLY
ISG15
LBH
LSP1
MT-CO1
MYL12A
RASGEF1B
SRGN
TNFAIP3
TSC22D3
VPS37B
YPEL5
TMSB10
AIMP1
IL2RB
PIM1
SSH2
IFRD1
MIER1
"
gene_list <- strsplit(gene_list, "\n")[[1]]
gene_list <- trimws(gene_list)
gene_list[2]


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
  features = c("IL7R","CD44","CD69","ENTPD1"), 
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
  features = c("IL7R","CD44","CD69","ENTPD1"),
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
# exhaustion
VlnPlot(
  cd4_cells, 
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
ggsave(filename = "vln_activation_exhastion_cluster_by_treatment_1.png" , path = "figures/cd4_cells/", dpi=300, width=12, height=10)


VlnPlot(
  cd4_cells,
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
ggsave(file = "figures/cd4_cells/vln_activation_exhastion_by_treatment_1.png", dpi=300, width=6, height=10)



# effector markers
VlnPlot(
  cd4_cells, 
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
ggsave(file = "figures/cd4_cells/Vln_T_effector_cluster_by_treatment_1.png", dpi=300, width=12, height=10, limitsize=FALSE)


VlnPlot(
  cd4_cells,
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
ggsave(file = "figures/cd4_cells/Vln_T_effector_by_treatment_1.png", dpi=300, width=6, height=10, limitsize=FALSE)


# naive markers
VlnPlot(
  cd4_cells, 
  features = c("IL7R","CD44","CD69","ENTPD1"), 
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
ggsave(file = "figures/cd4_cells/Vln_T_naive_cluster_by_treatment_1.png", dpi=300, width=12, height=10, limitsize=FALSE)


VlnPlot(
  cd4_cells,
  features = c("IL7R","CD44","CD69","ENTPD1"),
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
ggsave(file = "figures/cd4_cells/Vln_T_naive_by_treatment_1.png", dpi=300, width=6, height=10, limitsize=FALSE)

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
  features = c("IL7R","CD44","CD69","ENTPD1"), 
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
  features = c("IL7R","CD44","CD69","ENTPD1"),
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
  features = c("IL7R","CD44","CD69","ENTPD1"), 
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
  features = c("IL7R","CD44","CD69","ENTPD1"),
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
T_cells$cell_type <- ifelse(T_cells@assays$RNA@data["CD8A", ] > 0.2 | T_cells@assays$RNA@data["CD8B", ] > 0.2, 
                            "CD8", 
                            ifelse(T_cells@assays$RNA@data["CD4", ] > 0.2, 
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

#M# calculate ratio per patient
df_wide <- df %>%
  pivot_wider(names_from = Var2, values_from = Freq)
calculate_ratio <- function(sample_id) {
  cd4_count <- df_wide[df_wide$Var1 == sample_id, "CD4"]
  cd8_count <- df_wide[df_wide$Var1 == sample_id, "CD8"]
  if (cd4_count > 0) {
    ratio <- cd8_count / cd4_count
  } else {
    ratio <- NA
  }
  return(c(cd4_count = cd4_count, cd8_count = cd8_count, ratio = ratio))
}
patients <- unique(T_cells$Sample)
treatments <- sapply(patients, function(p) unique(subset(T_cells, Sample == p)$Treatment))
ratios_list <- lapply(patients, calculate_ratio)
ratios_df <- do.call(rbind, ratios_list)
ratios_df <- data.frame(Sample = patients, Treatment = treatments, ratios_df)
ratios_df$ratio <- as.numeric(ratios_df$ratio)
print(ratios_df)

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
ggsave(file = "figures/Tcells/cd8_cd4_ratio.png", dpi=300, width=4, height=6)


















# cd40lg zoom in --------
VlnPlot(cd4_cells, features = c("CD40LG"),
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
ggsave(file = "figures/cd4_cells/cd40lg_cluster_0.png", dpi=300, width=10, height=6)

VlnPlot(clus1, features = c("CD40LG"),
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
ggsave(file = "figures/clus1/cd40lg_cluster_1.png", dpi=300, width=10, height=6)

VlnPlot(clus2, features = c("CD40LG"),
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

VlnPlot(cd4_cells, features = c("ZNF831", "HIVEP1"),
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
ggsave(file = "figures/cd4_cells/0.png", dpi=300, width=10, height=6)

VlnPlot(clus1, features = c("ZNF831", "HIVEP1"),
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

VlnPlot(clus2, features = c("ZNF831", "HIVEP1"),
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











