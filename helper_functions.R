

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
gene_list <- "DDR1
CX3CL1
PLAT
LAMP1
KIT
LPL
SDC1
PEBP1
SLC2A4
VCAM1
FGG
CD24A
SCARB1
HSPA5
CD9
CD34
ADAM9
CD83
HMGB1
USP14
IL2RG
CD81
MIF
LAMP2
MCAM
HSP90AB1
M6PR
CD82
BGN
GABARAPL1
CRYAB
AIMP1
PDIA4
EPCAM
SFRP1
PLA2G1B
AMBP
SERPINE2
BACE2
PTPRU
APOH
DPP4
PLG
ATP5B
ADAM8
CLPTM1
LY6D
TRPV2
HSPA2
LPAR1
PDGFRB
LY6A
TMEM123
CD14
ENG
H13
TNFRSF1A
CAV2
F3
TGFB3
RAMP1
SERPINF2
FLOT2
TRPV4
IL2RB
CXCL12
SLC11A2
CD28
CD276
CALR
CD79B
SDC4
PVRL2
TJP1
ITGA2B
APOA4
PLA2G5
SYNJ2BP
FCGR1
GPR97
VLDLR
GHR
ADA
B4GALT1
EPHB6
NRP1
TRIP10
CAR4
TLR4
STX2
IL12RB1
RAMP2
THSD1
IL15
FCER1G
HBEGF
CD5
GPR124
ITGA7
CD97
TSPAN32
CAV3
SCNN1B
HSPA9
FURIN
TNFRSF12A
GPR162
ASTN1
NOTCH1
CXCL9
CAPN5
TIRAP
CD59A
PDGFA
LPAR3
ITGB7
CD55
CD2
TEK
FASL
SERPINA5
CHRNA1
PAM
CXCL10
CORIN
MPP3
SCNN1G
FGFBP1
HCST
PDGFC
ITGAX
TLR2
ACVR2B
CD163
CD3G
CD37
ICOSL
TNFRSF11A
CD96
ENPP1
FZD4
TNFRSF13C
PDPN
CTLA4
STAB2
SPAM1
CLEC2D
SELL
C3AR1
IL1R2
FLT3
ISG20
LGALS1
CD93
TNF
CCR1
BMPR2
CD4
PROM1
WNT6
CD7
CD274
RTN4R
CD22
GAB2
TREML2
P2RX7
CSF1R
TMX3
RC3H2
PSEN1
TNFRSF4
CD86
KLRB1B
DCBLD2
EPHA5
SELP
LPAR2
CFTR
GPR84
APP
HSPA8
TGFB1
ARNT2
KCNE1
IL17RB
IL2RA
CLEC7A
HPN
CD247
CHRNB2
KLRA5
ART2B
IL13
GDI2
SEMA4D
NTRK2
FGFR2
WNT7B
IL17RA
VEGFA
MFGE8
TIMP2
IL4RA
MRC2
GPC4
TRPC1
GPR56
ITGA4
PTPN11
ITGAV
CHRNA4
CIITA
TREM2
PTPRT
IL6ST
PECAM1
IL18RAP
KLRA2
H2-M3
CLEC5A
ANPEP
HHIP
CDON
KCNJ3
MS4A2
ITGB3
CD46
CEACAM2
GRIN2A
TRPM8
IL5RA
IL12RB2
TMC1
SLC6A2
CYSLTR2
CD1D2
NLGN1
CCR4
IL17A
GPR98
STRC
SELE
TNFSF4
BMP10
FCER1A
KLRC3
EFNA5
LDLR
BACE1
ABCA1
FGFR3
ADAM17
NR3C1
PDGFRA
EPHA4
ICOS
TGFA
L1CAM
NCAM1
ITGA3
CRHR2
TGFBR2
TNFRSF22
SEMA7A
ITGAM
CDH5
GABBR1
MSR1
WNT3A
ACVR1B
TNFRSF23
CD3E
FLT3L
FCER2A
PTPRC
SPN
PLAU
CLEC2I
C5AR1
CD200R1
SCN5A
CCR5
IL6RA
GLRA1
TACR1
CD40LG
CCR8
TNFRSF18
GP1BA
IL1RL1
ART1
IL15RA
ITGA6
GLRB
CR1L
WNT5B
ACHE
ADIPOQ
EBAG9
IRAK1BP1
TLR3
CNTNAP2
CXCR6
EPS8
CD3D
KCND2
CD84
TNFSF9
FZD5
HSPB1
FPR2
AGTRAP
TFRC
MME
FZD1
GFRA2
GYPA
IL1RN
HNRNPU
PAQR4
IDE
THY1
RALA
CD36
TNFRSF13B
MS4A1
ERP44
ITGA5
PTPRK
SLC34A1
SORT1
WNT7A
CLIC4
ADRB1
PDIA3
PGRMC1
CCR7
SEPT2
SLC46A2
JAM3
NID2
CDH13
ABCG1
S1PR1
TRPC4
AXL
BMP2
ATP6AP2
IFITM3
CD44
PVR
1600029D21RIK
PDLIM2
ICAM1
IGF2R
IFITM1
FGA
REEP2
CST8
ANGPTL3
PSTPIP1
SLIT2
BCAM
KCNMA1
BST2
EGFR
EMR4
KLRC1
SCNN1A
ACE2
FCGR4
PEAR1
CR2
CD8A
H2-K1
GREM1
ITGAL
WNT1
AQP4
KLRA8
H2-AB1
BMPR1A
CD74
H2-D1
ANXA5
SLAMF1
PTPRJ
AIPL1
TGFBR3
CEP290
KLHL20
NFAM1
HMMR
PSEN2
AMOT
IFNG
KLRB1F
PTPRR
CD72
CD209B
CD8B1
KLRA7
CD209A
ROBO4
ALCAM
SLC1A3
HSPD1
RYK
GPM6A
HSP90AA1
AGRN
LRPAP1
GPR125
BOC
ITGB1
PCSK6
ACE
ENOX2
ROBO1
CD48
MUC3
ITGB4
DSCAML1
FZD9
SHH
GHRHR
CD80
CCRL2
TNFRSF9
GPRC6A
KLRA1
FGB
ADAM10
NRXN1
REEP4
CD69
XPOT
RPS6KB1
CD99
ADCYAP1R1
PAQR3
HFE2
AQP11
HYAL4
RSPO2
CNRIP1
GUCY2G
SULF2
CD200R3
ULBP1
SCARA5
ANXA9
P2RY12
MUC16
APOE
INTU
NR4A2
CD38
ECE1
FERMT2
ATPIF1
ISLR2
CNTN2
P2RX2
GRIA1
H2-AA
CD200R4
VWF
FCGR2B
ACVRL1
SCUBE3
TLN1
TNR
SULF1
IRAK2
GRK5
WNT5A
RTN4RL1
NOTCH4
F2R
ACVR2A
PCSK9
P4HB
IL10RA
GRIN1
RGMA
5830411N06RIK
HYAL5
UMODL1
CD40
RTN4RL2
ADAM6A
ITGA1
LY6E
HEG1
FZD10
GPR116
UNC5D
CHRNA7
GPR174
WNT4
KCNH5
TNN
LAYN
DLK1
LRP1
TRPV1
NTRK1
CD226
GHSR
LBP
ITGAE
IFNGR1
CDH1
C1QBP
ADAM3
THBD
IL7R
CD53
FCGR3
ENPEP
HYAL2
CXCR4
ICAM2
OCLN
NCOR2
IL1R1
CD1D1
LY9
CD68
GPR65
MUC1
NDP
B2M
PTGER2
LY75
RAPSN
KDR
AOC3
KCNE2
IL27RA
KLRB1C
PDCD1
EPHB4
SCUBE1
IL4
VPREB1
CASR
LAG3
CXCR3
CD70
CD244
CX3CR1
PLA2R1
FGF22
KLRB1A
IL6
GABRR1
KLRC2
PDGFB
MRC1
IL21R
KISS1R
KLRK1
ITGA2
CD33
LY6F
CD19
ITGB2
FSHR
FUT4
TDGF1
FOLR1
LRP6
SFRP4
EMR1
CTSL
STX4A
HAVCR2
AMELX
FOLR2
FGF8
NOTCH2
CD6
SLC6A1
ADAMTS7
CD27
TNFRSF14
PLAUR
TREML1
GPIHBP1
GPR160
TMEM102
SLAMF7
ERP29
AAMP
NLGN2
PTPN3
BTLA
BSG
PPFIA2
PKD2L1
VTCN1
ROBO2
KLRE1
CD52
MSLN
KLRD1
FAS
THBS1
ARSA
ADAM2
ENTPD6
CD164
CD320
ENTPD5
CD248
CAP1
AMN
SIVA1
TRAF3
MS4A6B
CD47
ABCB1A
IGLL1
CD160
IDO1
PROCR
CD2AP
IL18R1
TNFRSF8
IL8RB
CXCR5
TNFRSF10B
ABCG2
ICAM4
ENTPD1
GYPC
CD99L2
PRR3
SLC3A2
CAST
NT5E
PGP
CCND2
CD3EAP
SCO1
DARC
SLC44A1
PTGFRN
NRP
LSM1
CTSD
PRNP
GSS
PTPRCAP
CD2BP2
PVRL3
CD200
CD302
IFNAR1
TLR1
CD5L
CCR6
PEMT
LAP3
PLXNC1
SCARB2
IL13RA1
SP1
CD151
TNFSF10
IGSF8
TIGIT
LILRB4
Gene
PDCD1LG2
CTLA2B
CTLA2A
IL12A
SPP1
TNFSF18
OSM
LIF
BMP3
WNT2
SLURP1
PRL7D1
GDF1
CRLF1
IL17F
GDF15
IL17C
GDF7
GDF6
GDF5
TSLP
GDF3
IFNL2
IFNL3
IL23A
CCL1
GREM2
IL1A
BMP15
IL1B
IFNK
IL27
BMP7
GDF11
EBI3
GPI1
IL7
TNFSF15
LEFTY2
LEFTY1
TNFSF8
CTF1
CTF2
CSF1
IL18
BMP6
IL17B
KITL
LTB
IL33
LTA
MSTN
BMP1
EDN1
NODAL
IL9
RPL13A
IL25
CMTM5
TNFSF14
CLCF1
IL31
TNFSF13B
CER1
FAM3B
BMP8B
BMP8A
IL1F8
IL1F9
IL1F6
IL1F5
IL1F10
CCL17
NAMPT
IL12B
IL22
ILTIFB
IL11
GRN
IFNB1
GM12597
IFNA14
IFNA9
IFNA12
IFNA13
GM13280
IFNA2
IFNAB
GM13271
GM13283
GM13290
TNFSF11
GM13289
GM13272
IFNZ
GM13276
GM13277
GM13278
GM13275
GM13279
GM13285
GM13287
GM13288
IFNA7
IFNA11
IFNA6
IFNA5
IFNA4
IFNA1
IFNE
CMTM2A
CMTM6
CMTM7
CMTM2B
CMTM8
CMTM3
CMTM4
C1QTNF4
SCGB3A1
IL16
IL17D
SCG2
GDF10
GDF2
PGLYRP1
CCL20
INHBA
IL34
AREG
TNFSF12
BC096441
TNFSF13
GDF9
IL5
THPO
CSF2
IL3
BMP4
CCL24
IL2
IL21
FGF2
CSF3
IL24
IL20
IL19
IL10
BMP5
CCL21C
CCL27A
CXCL13
GM21541
GM13304
GM13306
CCL28
CCL21B
GM10591
GM2564
CCL27B
CCL19
CCL21A
XCL1
CXCL16
CCL2
CCL7
CCL11
CCL12
CCL8
CCL5
CCL9
CCL6
CCL3
CCL4
CXCL14
CCL25
CCL22
CKLF
GRAMD2
CXCL5
PPBP
PF4
CXCL3
CXCL15
CXCL1
CXCL2
CXCL11
CCL26
TRPM4
ARHGEF5
RETNLG
RPS19
FLT1
MYO9B
CALCA
PTPRO
RAC1
CXADR
PRKCA
SYK
SLC37A4
AMICA1
NCKAP1L
TGFB2
EDN3
EDN2
S100A8
S100A9
CSF3R
CXCR2
ITGA9
PDE4B
PDE4D
CORO1A
LYST
SBDS
CCR2
GAS6
HRH1
NUP85
EDNRB
ROCK1
MSN
EZR
OLR1
FERMT3
TNIP1
GCNT1
PODXL2
LEP
SELPLG
GOLPH3
CHST4
STK10
FN1
IGHG2C
PLA2G2A
REG4
F11
PTGDS
KLKB1
OLFM4
BGLAP2
SPOCK3
C8G
SERPINC1
OLFM1
CTRB1
OGN
C1QTNF7
GPLD1
EGF
IL18BP
UCMA
CFI
MMP1A
MMP8
APOC2
IAPP
PTPRG
C8B
GIF
COL6A2
C8A
SERPINE3
MYOC
ADAMTS20
F7
FGF10
CTS7
SERPINB10
CTSB
WNT9A
NEPN
POMC
APOD
PRL3D1
CEL
COL25A1
PRL3B1
BCHE
CNTF
CEACAM10
SMPD1
HPX
WFDC1
TFF2
FBLN1
SERPINI2
TFF1
SEZ6
FBLN5
ADCYAP1
F5
CNP
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

VlnPlot(cd4_cells, features = c("TGFBR3", "HIVEP1"),
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
clus1 = readRDS("objects/clus1_no_cd4_14_0.5.rds")
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
  ggsave(file = paste0("figures/clus1/de_gene_search/","run_", i,  ".png"), dpi=300, width=16, height=12)
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
  ggsave(file = paste0("figures/clus2/de_gene_search/","run_", i,  ".png"), dpi=300, width=16, height=12)
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
