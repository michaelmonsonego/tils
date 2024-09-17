

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
library(scales)


setwd("D:/Michael/git_check/tils")
# setwd("C:/Users/michael monsonego/Documents/tils") #M# update between computers(git)

# use_python((r"(C:\Users\user\AppData\Local\Programs\Python\Python311)"), required = TRUE)
# source("//asafmadilabfs/asafmadilab$/michael/Madi lab/Signature_propo/signature_utils.R")
# source('D:/Michael/git_check/3-groups/sig_genes.R')



# TGFb up signature -------------------------------

up_tgfb <- "ACVR1	APC	ARID4B	BCAR3	BMP2	BMPR1A	BMPR2	CDH1	CDK9	CDKN1C	CTNNB1	ENG	FKBP1A	FNTA	FURIN	HDAC1	HIPK2	ID1	ID2	ID3	IFNGR2	JUNB	KLF10	LEFTY2	LTBP2	MAP3K7	NCOR2	NOG	PMEPA1	PPM1A	PPP1CA	PPP1R15A	RAB31	RHOA	SERPINE1	SKI	SKIL	SLC20A1	SMAD1	SMAD3	SMAD6	SMAD7	SMURF1	SMURF2	SPTBN1	TGFB1	TGFBR1	TGIF1	THBS1	TJP1	TRIM33	UBE2D3	WWTR1	XIAP
"

up_tgfb <- strsplit(up_tgfb, "\t")[[1]]
up_tgfb <- toupper(up_tgfb)
up_tgfb <- trimws(up_tgfb)

up_tgfb <- list(up_tgfb)
clus2 = AddModuleScore(object = clus2, features = up_tgfb, name = "up_tgfb", assay = "RNA")
a <- FeaturePlot(object = clus2, features = "up_tgfb1", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "up_tgfb", subtitle = "up_tgfb")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/clus2/clus_2_up_tgfb.png", dpi=300, width=5, height=5)


gglist <-  list()
name <- "up_tgfb1"
DefaultAssay(clus2) <- "RNA" 
for(clus in levels(clus2$seurat_clusters)){
  print(clus)
  obj <- subset(clus2, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'up_tgfb1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus2/clus_2_up_tgfb_vln.png", dpi=300, width=16, height=10)

#M# vlnplot for all clusters combined 
SignatureScore(clus2, 'up_tgfb1')
ggsave(file = "figures/clus2/clus2_up_tgfb_vln_all_clusters.png", dpi=300, width=8, height=6)


# clus 1
clus1 = AddModuleScore(object = clus1, features = up_tgfb, name = "up_tgfb", assay = "RNA")
a <- FeaturePlot(object = clus1, features = "up_tgfb1", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "up_tgfb", subtitle = "up_tgfb")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/clus1/clus_1_up_tgfb.png", dpi=300, width=5, height=5)


gglist <-  list()
name <- "up_tgfb1"
DefaultAssay(clus1) <- "RNA" 
for(clus in levels(clus1$seurat_clusters)){
  print(clus)
  obj <- subset(clus1, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'up_tgfb1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus1/clus_1_up_tgfb_vln.png", dpi=300, width=16, height=10)

# cluster 0
cd4_cells = AddModuleScore(object = cd4_cells, features = up_tgfb, name = "up_tgfb", assay = "RNA")
a <- FeaturePlot(object = cd4_cells, features = "up_tgfb1", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "up_tgfb", subtitle = "up_tgfb")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/cd4_cells/clus_0_up_tgfb.png", dpi=300, width=5, height=5)


gglist <-  list()
name <- "up_tgfb1"
DefaultAssay(cd4_cells) <- "RNA" 
for(clus in levels(cd4_cells$seurat_clusters)){
  print(clus)
  obj <- subset(cd4_cells, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'up_tgfb1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/cd4_cells/clus_0_up_tgfb_vln.png", dpi=300, width=16, height=10)


# IFNg down signature cluster 2-------------------------------
down_ifng <- "Akr1c18
Apex1
BC094916
Bpgm
Capg
Cbx3-ps6
Ccr4
Cnga1
Dgcr8
Dhx9
Fabp5
Fam101b
Fgf2
Ftsj3
Gar1
Gbp2
Gm13237
Gm13577
Gm14085
Gm3636
Gm3851
Gpatch4
Grwd1
Hmgb1
Hmgn5
Hnrnpd
Id3
Ifi47
Igfbp7
Iigp1
Il7r
Itgav
Klrd1
Mat2a
Mlf1
Mrpl13
Mxd4
Nop58
Nrn1
Oxct1
Pdzk1ip1
Pinx1
Plekhb2
Pno1
Polr1b
Prdx1
Qtrtd1
Rab11fip4
Rcsd1
RP23-77B7.1
Rrp15
Slc9b2
Socs2
Sytl2
Thy1
Tiam1
Timm10
Tma16
Tmem38b
Tns1
Tsr2
Tyms
Utp20
Vav2
Xcl1
"

down_ifng <- strsplit(down_ifng, "\n")[[1]]
down_ifng <- toupper(down_ifng)
down_ifng <- trimws(down_ifng)

down_ifng <- list(down_ifng)
clus2 = AddModuleScore(object = clus2, features = down_ifng, name = "down_ifng", assay = "RNA")
a <- FeaturePlot(object = clus2, features = "down_ifng1", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "down_ifng", subtitle = "down_ifng")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/clus2/clus_2_down_ifng.png", dpi=300, width=5, height=5)


gglist <-  list()
name <- "down_ifng1"
DefaultAssay(clus2) <- "RNA" 
for(clus in levels(clus2$seurat_clusters)){
  print(clus)
  obj <- subset(clus2, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'down_ifng1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus2/clus_2_down_ifng_vln.png", dpi=300, width=16, height=10)

#M# vlnplot for all clusters combined 
SignatureScore(clus2, 'down_ifng1')
ggsave(file = "figures/clus2/clus_2_down_ifng_vln_all_clusters.png", dpi=300, width=8, height=6)

# IFNg hallmark signature -------------------------------

up_ifng <- "ADAR
APOL6
ARID5B
ARL4A
AUTS2
B2M
BANK1
BATF2
BPGM
BST2
BTG1
C1R
C1S
CASP1
CASP3
CASP4
CASP7
CASP8
CCL2
CCL5
CCL7
CD274
CD38
CD40
CD69
CD74
CD86
CDKN1A
CFB
CFH
CIITA
CMKLR1
CMPK2
CMTR1
CSF2RB
CXCL10
CXCL11
CXCL9
DDX58
DDX60
DHX58
EIF2AK2
EIF4E3
EPSTI1
FAS
FCGR1A
FGL2
FPR1
GBP4
GBP6
GCH1
GPR18
GZMA
HELZ2
HERC6
HIF1A
HLA-A
HLA-B
HLA-DMA
HLA-DQA1
HLA-DRB1
HLA-G
ICAM1
IDO1
IFI27
IFI30
IFI35
IFI44
IFI44L
IFIH1
IFIT1
IFIT2
IFIT3
IFITM2
IFITM3
IFNAR2
IL10RA
IL15
IL15RA
IL18BP
IL2RB
IL4R
IL6
IL7
IRF1
IRF2
IRF4
IRF5
IRF7
IRF8
IRF9
ISG15
ISG20
ISOC1
ITGB7
JAK2
KLRK1
LAP3
LATS2
LCP2
LGALS3BP
LY6E
LYSMD2
MARCHF1
METTL7B
MT2A
MTHFD2
MVP
MX1
MX2
MYD88
NAMPT
NCOA3
NFKB1
NFKBIA
NLRC5
NMI
NOD1
NUP93
OAS2
OAS3
OASL
OGFR
P2RY14
PARP12
PARP14
PDE4B
PELI1
PFKP
PIM1
PLA2G4A
PLSCR1
PML
PNP
PNPT1
PSMA2
PSMA3
PSMB10
PSMB2
PSMB8
PSMB9
PSME1
PSME2
PTGS2
PTPN1
PTPN2
PTPN6
RAPGEF6
RBCK1
RIPK1
RIPK2
RNF213
RNF31
RSAD2
RTP4
SAMD9L
SAMHD1
SECTM1
SELP
SERPING1
SLAMF7
SLC25A28
SOCS1
SOCS3
SOD2
SP110
SPPL2A
SRI
SSPN
ST3GAL5
ST8SIA4
STAT1
STAT2
STAT3
STAT4
TAP1
TAPBP
TDRD7
TNFAIP2
TNFAIP3
TNFAIP6
TNFSF10
TOR1B
TRAFD1
TRIM14
TRIM21
TRIM25
TRIM26
TXNIP
UBE2L6
UPP1
USP18
VAMP5
VAMP8
VCAM1
WARS1
XAF1
XCL1
ZBP1
ZNFX1
"

up_ifng <- strsplit(up_ifng, "\n")[[1]]
up_ifng <- toupper(up_ifng)
up_ifng <- trimws(up_ifng)

up_ifng <- list(up_ifng)

# clus 2
clus2 = AddModuleScore(object = clus2, features = up_ifng, name = "up_ifng", assay = "RNA")
a <- FeaturePlot(object = clus2, features = "up_ifng1", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "up_ifng", subtitle = "up_ifng")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/clus2/clus_2_up_ifng.png", dpi=300, width=5, height=5)


gglist <-  list()
name <- "up_ifng1"
DefaultAssay(clus2) <- "RNA" 
for(clus in levels(clus2$seurat_clusters)){
  print(clus)
  obj <- subset(clus2, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'up_ifng1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus2/clus_2_up_ifng_vln.png", dpi=300, width=16, height=10)

#M# vlnplot for all clusters combined 
SignatureScore(clus2, 'up_ifng1')
ggsave(file = "figures/clus2/clus_2_up_ifng_vln_all_clusters.png", dpi=300, width=8, height=6)


# clus 1
clus1 = AddModuleScore(object = clus1, features = up_ifng, name = "up_ifng", assay = "RNA")
a <- FeaturePlot(object = clus1, features = "up_ifng1", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "up_ifng", subtitle = "up_ifng")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/clus1/clus_1_up_ifng.png", dpi=300, width=5, height=5)


gglist <-  list()
name <- "up_ifng1"
DefaultAssay(clus1) <- "RNA" 
for(clus in levels(clus1$seurat_clusters)){
  print(clus)
  obj <- subset(clus1, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'up_ifng1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/clus1/clus_1_up_ifng_vln.png", dpi=300, width=16, height=10)


# clus 0
cd4_cells = AddModuleScore(object = cd4_cells, features = up_ifng, name = "up_ifng", assay = "RNA")
a <- FeaturePlot(object = cd4_cells, features = "up_ifng1", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "up_ifng", subtitle = "up_ifng")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/cd4_cells/clus_0_up_ifng.png", dpi=300, width=5, height=5)


gglist <-  list()
name <- "up_ifng1"
DefaultAssay(cd4_cells) <- "RNA" 
for(clus in levels(cd4_cells$seurat_clusters)){
  print(clus)
  obj <- subset(cd4_cells, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'up_ifng1') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/cd4_cells/clus_0_up_ifng_vln.png", dpi=300, width=16, height=10)



# telomere signature ------------------------------------------------------
clus1 = readRDS("objects/clus1_no_cd4_second_sub_14_0.5.rds")
clus2 = readRDS("objects/tils_2.rds")

#M# downregulated in cd8 telomere signature-------------
down_cd8 <- "ADTRP
AIM1
BIRC3
CCNG1
DUSP4
EEF1A1
EIF1
FTH1
GIMAP5
HINT1
KLF6
LINC00152
PDCD4
PIM3
PPP6C
RASA3
RPL12
RPL13
RPL29
RPL39
RPS18
RPS8
SESN3
SMS
TMEM66
ZFP36
"

down_cd8 <- strsplit(down_cd8, "\n")[[1]]
down_cd8 <- trimws(down_cd8)

stem_cd8 <- list(down_cd8)
clus2 = AddModuleScore(object = clus2, features = stem_cd8, name = "stem_cd8", assay = "RNA")
a <- FeaturePlot(object = clus2, features = "stem_cd81", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "stem_cd8", subtitle = "stem_cd8")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/signatures/clus_2_down_stem_cd8.png", dpi=300, width=5, height=5)

stem_cd8 <- list(down_cd8)
clus1 = AddModuleScore(object = clus1, features = stem_cd8, name = "stem_cd8", assay = "RNA")
a <- FeaturePlot(object = clus1, features = "stem_cd81", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "stem_cd8", subtitle = "stem_cd8")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/signatures/clus_1_down_stem_cd8.png", dpi=300, width=5, height=5)


gglist <-  list()
name <- "stem_cd81"
DefaultAssay(clus1) <- "RNA" 
for(clus in levels(clus1$seurat_clusters)){
  print(clus)
  obj <- subset(clus1, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'stem_cd81') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/signatures/stem_cd8_clus_1_down_vln.png", dpi=300, width=16, height=10)

#M# vlnplot for all clusters combined 
SignatureScore(clus1, 'stem_cd81')
ggsave(file = "figures/clus1/clus_1_down_telomere_vln_all_clusters.png", dpi=300, width=8, height=6)

gglist <-  list()
name <- "stem_cd81"
DefaultAssay(clus2) <- "RNA" 
for(clus in levels(clus2$seurat_clusters)){
  print(clus)
  obj <- subset(clus2, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'stem_cd81') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/signatures/stem_cd8_clus_2_down_vln.png", dpi=300, width=16, height=10)

#M# vlnplot for all clusters combined 
SignatureScore(clus2, 'stem_cd81')
ggsave(file = "figures/clus2/clus2_down_telomere_vln_all_clusters.png", dpi=300, width=8, height=6)




#M# upregulated in cd8 telomere signature---------
# DNAJA1 #M# memory cd8 DE gene for telomere legth

up_cd8 <- "CORO1A
DBI
DNAJA1
FGFBP2
GADD45B
HSPH1
IFITM2
ITGB1BP1
ITGB2
NDUFA12
NDUFB7
NR4A2
PAIP2
PTP4A1
RANBP2
RPS26
RSRC2
SKAP1
SRP14
SSBP4
TBC1D10C
TMBIM1
TMEM59
TMSB4X
TRA2B
VAMP2
ZNF331
"
up_cd8 <- strsplit(up_cd8, "\n")[[1]]
up_cd8 <- trimws(up_cd8)

stem_cd8 <- list(up_cd8)
clus2 = AddModuleScore(object = clus2, features = stem_cd8, name = "stem_cd8", assay = "RNA")
a <- FeaturePlot(object = clus2, features = "stem_cd81", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "stem_cd8", subtitle = "stem_cd8")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/signatures/clus_2_up_stem_cd8.png", dpi=300, width=5, height=5)

stem_cd8 <- list(up_cd8)
clus1 = AddModuleScore(object = clus1, features = stem_cd8, name = "stem_cd8", assay = "RNA")
a <- FeaturePlot(object = clus1, features = "stem_cd81", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "stem_cd8", subtitle = "stem_cd8")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/signatures/clus_1_up_stem_cd8.png", dpi=300, width=5, height=5)

gglist <-  list()
name <- "stem_cd81"
DefaultAssay(clus2) <- "RNA" 
for(clus in levels(clus2$seurat_clusters)){
  print(clus)
  obj <- subset(clus2, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'stem_cd81') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/signatures/stem_cd8_clus_2_up_vln.png", dpi=300, width=16, height=10)

#M# vlnplot for all clusters combined 
SignatureScore(clus2, 'stem_cd81')
ggsave(file = "figures/clus2/clus2_up_telomere_vln_all_clusters.png", dpi=300, width=8, height=6)


gglist <-  list()
name <- "stem_cd81"
DefaultAssay(clus1) <- "RNA" 
for(clus in levels(clus1$seurat_clusters)){
  print(clus)
  obj <- subset(clus1, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'stem_cd81') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/signatures/stem_cd8_clus_1_up_vln.png", dpi=300, width=16, height=10)

#M# vlnplot for all clusters combined 
SignatureScore(clus1, 'stem_cd81')
ggsave(file = "figures/clus1/clus1_up_telomere_vln_all_clusters.png", dpi=300, width=8, height=6)

#M# cd4 upregulated genes telomere signature ---------
cd4_cells = readRDS("objects/tils_0_.35.rds")

up_cd4_gene_list <- "CD3G
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
IFRD1
MIER1
"
up_cd4_gene_list <- strsplit(up_cd4_gene_list, "\n")[[1]]
up_cd4_gene_list <- trimws(up_cd4_gene_list)


stem_cd4 <- list(up_cd4_gene_list)
cd4_cells = AddModuleScore(object = cd4_cells, features = stem_cd4, name = "stem_cd4", assay = "RNA")
a <- FeaturePlot(object = cd4_cells, features = "stem_cd41", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "stem_cd4", subtitle = "stem_cd4")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/signatures/cd4_stem_up.png", dpi=300, width=5, height=5)


gglist <-  list()
name <- "stem_cd41"
DefaultAssay(cd4_cells) <- "RNA" 
for(clus in levels(cd4_cells$seurat_clusters)){
  print(clus)
  obj <- subset(cd4_cells, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'stem_cd41') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/signatures/stem_cd4_up_vln.png", dpi=300, width=16, height=10)


#down cd4 telomere signature -----------
down_cd4 <- "ARID5A
BATF
BAZ1A
CCR7
CMSS1
DDX21
GIMAP4
GIMAP7
IL7R
MTHFD2
MYC
NAP1L1
PABPC1
PIM2
PSMC3
RPL36
RPS29
SELL
SOCS3
AIMP1
IL2RB
PIM1
SSH2
"
down_cd4 <- strsplit(down_cd4, "\n")[[1]]
down_cd4 <- trimws(down_cd4)


stem_cd4 <- list(down_cd4)
cd4_cells = AddModuleScore(object = cd4_cells, features = stem_cd4, name = "stem_cd4", assay = "RNA")
a <- FeaturePlot(object = cd4_cells, features = "stem_cd41", reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"))+labs(title = "stem_cd4", subtitle = "stem_cd4")+ theme(plot.subtitle = element_text(hjust = 0.5))
a
ggsave(file = "figures/signatures/cd4_stem_down.png", dpi=300, width=5, height=5)

gglist <-  list()
name <- "stem_cd41"
DefaultAssay(cd4_cells) <- "RNA" 
for(clus in levels(cd4_cells$seurat_clusters)){
  print(clus)
  obj <- subset(cd4_cells, subset = seurat_clusters == clus)
  p <- SignatureScore(obj, 'stem_cd41') + ggtitle(obj@active.ident)
  gglist[[(as.numeric(clus)+1)]] <- p
}
cowplot::plot_grid(plotlist = gglist, ncol = 4, nrow = 2) + 
  ggtitle(name)
ggsave(file = "figures/signatures/stem_cd4_down_vln.png", dpi=300, width=16, height=10)


# signature score function ------------------------------------------------
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
               cols= c("#ff7f0e", "#1f77b4")) +  theme_classic(base_size = 14) +
    theme(text = element_text(size=18, colour = "black")) + RotatedAxis() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    labs(title = "", y = name,  x="") + theme(legend.position="right") +
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") +
    annotate("text", x = 1.5, y=annot.max.y, label = label, size = 6) +
    scale_y_continuous(labels = comma)
  return (p)
  
}


# signatures(effector, il2..) -----------------------------------------------------------------
#M# clus 1 sigs
clus1_with_sigs = readRDS("objects/clus1_with_sigs.rds")
clus_2_with_sigs <- readRDS("objects/tils_2_with_sigs.rds")

#M# sig up genes from paper
# persistance_up <- c('KLRB1','ZNF683','ITGB1','C1orf162','IL7R','LIME1', "S1PR1", "TIMP1", "C10orf54", "TBXAS1", "KLF2", "LTB", "UBXN11", "CD40LG", "AMICA1", "FAM65B", "VCL", "RASA3", "SCML4", "MYC", "P2RY8")
# persistance_up_sig <- make_signature(clus1_with_sigs, persistance_up,idents = c("Responder","Non_Responder"), format_flag = FALSE)
# persistance_up_sig_obj <- persistance_up_sig[[1]]
# persistance_up_sig[[3]]
# ggsave("figures/clus1_with_sigs/co_stim_sig_box.png", width = 20, height = 20, units = "cm")
# persistance_up_sig[[2]]
# ggsave("figures/clus1_with_sigs/co_stim_sig_tsne.png", width = 20, height = 20, units = "cm")


#M# test dor addmodulescore() function from vis.rmd
# violin plot per treatment
# VlnPlot(clus1_with_sigs, features = "Tregs1", group.by = "Treatment", pt.size = 0)+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+
#   geom_boxplot(alpha=0.3,show.legend = FALSE)
# ggsave("figures/clus1_with_sigs/violin_sig_co_stimulatory_by_treatment.png", width = 20, height = 20, units = "cm")

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

SignatureScore(clus1, 'effector_memory1')
ggsave(file = "figures/clus1/clus1_effector_memory_vln_all_clusters.png", dpi=300, width=8, height=6)


DefaultAssay(clus1_with_sigs) <- "integrated"


# ____

#M# same for cluster 2 : 
clus_2_with_sigs <- readRDS("objects/tils_2_with_sigs.rds")

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

#M# vlnplot for all clusters combined 
SignatureScore(clus2, 'il21')
ggsave(file = "figures/clus2/clus_2_il21_vln_all_clusters.png", dpi=300, width=8, height=6)



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

#M# vlnplot for all clusters combined 
SignatureScore(clus2, 'persistance_up1')
ggsave(file = "figures/clus2/clus_2_persistance_up_vln_all_clusters.png", dpi=300, width=8, height=6)


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

SignatureScore(clus2, 'effector_memory1')
ggsave(file = "figures/clus2/clus_2_effector_memory_vln_all_clusters.png", dpi=300, width=8, height=6)


# another section ---------------------------------------------------------
























































