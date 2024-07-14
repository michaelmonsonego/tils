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

clus2_top_10 <- merged_data %>% 
  top_n(n=10, wt=Cluster_2) %>% 
  rownames()

clus_2_lowest_10 <- merged_data %>% 
  top_n(n=-10, wt=Cluster_2) %>% 
  rownames()

selected_genes <- c(clus_2_lowest_10, clus2_top_10)

filtered_data <- long_data %>%
  filter(Gene %in% selected_genes)

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
  rownames()

clus_0_lowest_10 <- merged_data %>% 
  top_n(n=-10, wt=Cluster_0) %>% 
  rownames()

selected_genes <- c(clus_0_lowest_10, clus0_top_10)

filtered_data <- long_data %>%
  filter(Gene %in% selected_genes)

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
  rownames()

clus_1_lowest_10 <- merged_data %>% 
  top_n(n=-10, wt=Cluster_1) %>% 
  rownames()

selected_genes <- c(clus_1_lowest_10, clus_1_top_10)

filtered_data <- long_data %>%
  filter(Gene %in% selected_genes)

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

