## This script performs the analyses for ssGSEA using RNAseq data from the 
## MONSTAR2 cohort. ssGSEA NES scores for hallmark pathways were generated for 
## all mss crc patients and comparisons between the early onset versus late onset
## cancers were made.

library(plyr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(tidyverse)
library(plotrix)
library(GSVA)
library(fgsea)
library(purrr)
library(heatmaply)
library(progeny)
library(factoextra)
library(janitor)
library(ggrepel)

setwd("/fsx/home/vlee/Projects/MONSTAR2/EO_vs_LO_CRC/")
hallmark <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt")
key <- read.csv("/fsx/home/vlee/Projects/MONSTAR2/Analysis_dataset/WTSExpression_SUBJID_TSB_cancertype_oncotree_key.csv")
eo_list <- unlist(read.csv("eocrcdenovomss_list.csv"))
lo_list <- unlist(read.csv("locrcdenovomss_list.csv"))
key <- key%>%
  mutate(Onset = if_else(Tumor_Sample_Barcode %in% eo_list, "EO", 
                         if_else(Tumor_Sample_Barcode%in% lo_list, "LO", "")))%>%
  mutate(NewID = gsub("-", ".", Tumor_Sample_Barcode))


eo_wts_files <- key%>%
  filter(Tumor_Sample_Barcode %in%eo_list)%>%
  filter(Version %in% c("a", "b"))%>%
  select(filename)%>%
  unlist()
eo_list2 <- key%>%
  filter(Tumor_Sample_Barcode %in%eo_list)%>%
  filter(Version %in% c("a", "b"))%>%
  select(Tumor_Sample_Barcode)%>%
  unlist()
lo_wts_files <- key%>%
  filter(Tumor_Sample_Barcode%in% lo_list)%>%
  filter(Version %in% c("a", "b"))%>%
  select(filename)%>%
  unlist()
lo_list2 <- key%>%
  filter(Tumor_Sample_Barcode%in% lo_list)%>%
  filter(Version %in% c("a", "b"))%>%
  select(Tumor_Sample_Barcode)%>%
  unlist()
mss_crc <- c(eo_wts_files, lo_wts_files)
eolo_list <- c(eo_list2, lo_list2)

genelist <- unlist(read.csv("/fsx/home/vlee/Projects/MONSTAR2/Analysis_dataset/WTS_common_genes.csv"))

## Combine MSS_CRC files
crc_datamatrix <- matrix(genelist)
colnames(crc_datamatrix) <- "X1"
setwd("/fsx/home/vlee/Projects/MONSTAR2/MONSTAR_genome_data_20240617/")
for (e in mss_crc){
  sample <- read_delim(e, comment = "#", col_names = FALSE, show_col_types = FALSE)%>%
    filter(X1 %in% genelist)%>%
    mutate(log2tpm = log2(X2+1))%>%
    select(c(X1, log2tpm))
  colnames(sample) <- c("X1", substr(e, 1, 11))
  crc_datamatrix <- merge(crc_datamatrix, sample, by = "X1", all.x=TRUE)
  
}


datamatrix <- crc_datamatrix[,-1]
datamatrix <- as.matrix(datamatrix)
rownames(datamatrix) <- genelist


annot_onset <- key%>%
  filter(filename %in% mss_crc)%>%
  select(Tumor_Sample_Barcode, Onset)


anotcol_onset = list(Onset = c( EO=  "firebrick",
                               LO = "dodgerblue"))

## GSVA ssGSEA method, hallmark and progeny_geneset

data_param <-ssgseaParam(datamatrix, hallmark) 
set.seed(178)
hallmark_res <- gsva(data_param)

sample_order <- match(colnames(hallmark_res), annot_onset$Tumor_Sample_Barcode)
annot_onset1 <- annot_onset[sample_order,]

annot_onset2 <- as.data.frame(annot_onset1[,-1])
rownames(annot_onset2) <- annot_onset1$Tumor_Sample_Barcode
colnames(annot_onset2) <- "Onset"

pheatmap(t(hallmark_res),
         annotation_row = annot_onset2,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_colors = anotcol_onset)

## Check ssGSEA result on PCA

res.pca <- prcomp(t(hallmark_res), scale. = TRUE )
fviz_eig(res.pca)
pca <- as.data.frame(res.pca$x)%>%
  select(c(PC1, PC2))
pca$Tumor_Sample_Barcode <- rownames(pca)
pca <- merge(pca, key, by ="Tumor_Sample_Barcode", all.x = TRUE)
nes_pca <- ggplot(pca, aes(x = PC1, y = PC2, colour =Version))+geom_point()+
  theme_light()+
  ggtitle("EOvsLO PCA plot of all ssGSEA NES")+
  theme(axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 22))
nes_pca


library(reshape2)

hallmark_res_df <- melt(t(hallmark_res))
colnames(hallmark_res_df) <- c("Tumor_Sample_Barcode", "Pathway", "NES")
pathway <- unique(hallmark_res_df$Pathway)
##annotate cancer type 

##generate ttest results for EO vs LO
##annotate onset
hallmark_res_df <- hallmark_res_df%>%
  mutate(onset = if_else(Tumor_Sample_Barcode %in% eo_list, "EO", "LO"))
hallmark_res_df$onset <- as.factor(hallmark_res_df$onset)


ttest_result_eolo <- c()

for (p in pathway){
  df <- hallmark_res_df%>%
    filter(Pathway ==p)
  test <- t.test(df$NES~df$onset)
  
  ttest_result_eolo <- rbind(ttest_result_eolo, cbind(p, test$statistic, test$p.value, test$estimate[1], test$estimate[2]))
  
}
ttest_result_eolo <- as.data.frame(ttest_result_eolo)
colnames(ttest_result_eolo) <- c("Pathway", "statistic", "pval", "mean_EO", "mean_LO")
ttest_result_eolo$statistic <- as.numeric(ttest_result_eolo$statistic)
ttest_result_eolo$pval <- as.numeric(ttest_result_eolo$pval)
ttest_result_eolo$mean_EO <- as.numeric(ttest_result_eolo$mean_EO)
ttest_result_eolo$mean_LO <- as.numeric(ttest_result_eolo$mean_LO)
eolo_lab <- c("HALLMARK_E2F_TARGETS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
              "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_BILE_ACID_METABOLISM")
ttest_result_eolo <- ttest_result_eolo%>%
  mutate(NES_delta = mean_EO-mean_LO)%>%
  mutate(`-log10pval` = -log10(pval))%>%
  mutate(NES_Score= if_else(NES_delta>0.05 &pval<0.05, "UP",
                            if_else(NES_delta< -0.05&pval<0.05, "DOWN", "N.S")))%>%
  mutate(label = if_else(Pathway %in% eolo_lab, Pathway, ""))

## Plot results on a volcano plot

eolo_volcplot <- ggplot(ttest_result_eolo, aes(x = NES_delta, y = `-log10pval`, col = NES_Score, label= label ))+
  geom_point(size = 4)+
  geom_text_repel()+
  scale_color_manual(values = c("black", "lightgrey", "dodgerblue3"))+
  geom_hline(yintercept = -log10(0.05), col = "red")+
  geom_vline(xintercept = c(-0.1, 0.1), col= "red")+
  theme_light()+
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13))+
  xlim(c(-0.02, 0.02))+
  ylim(c(0,1.3))

eolo_volcplot

