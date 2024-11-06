##This script performs the analyses for ssGSEA using RNAseq data from the 
##MONSTAR2 cohort. ssGSEA NES scores for hallmark pathways were generated for 
##breast HR+HER2- and all colorectal adenocarcinoma patients.

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

## Get HR and HER2 status of BC patients
pf6 <- read.csv("/fsx/home/vlee/Projects/MONSTAR2/MONSTAR2_clinical_June2024/PF6.csv")
pf7 <- read.csv("/fsx/home/vlee/Projects/MONSTAR2/MONSTAR2_clinical_June2024/PF7.csv")
hrpos <- pf7%>%
  filter(PFITEM1==1|PFITEM2==1)%>%
  select(SUBJID)%>%
  unlist()
her2neg <- pf6%>%
  filter(PFITEM1==0|PFITEM2==0|PFITEM2==1|PFITEM3==2)%>%
  select(SUBJID)%>%
  unlist()

hrposher2neg <- intersect(hrpos, her2neg)

cancers <- c("Breast cancer", "Colorectal cancer")
bc <- key%>%
  filter(CANCERTYPE =="Breast cancer")%>%
  filter(SUBJID %in% hrposher2neg)%>%
  filter(Version %in% c("a", "b"))%>%
  select(filename)%>%
  unlist()
crc <- key%>%
  filter(CANCERTYPE =="Colorectal cancer")%>%
  filter(Version %in% c("a", "b"))%>%
  select(filename)%>%
  unlist()
bio_control <- c(bc, crc)

genelist <- unlist(read.csv("/fsx/home/vlee/Projects/MONSTAR2/Analysis_dataset/WTS_common_genes.csv"))

## Combine BC and CRC files
biocontrol_datamatrix <- matrix(genelist)
colnames(biocontrol_datamatrix) <- "X1"
setwd("/fsx/home/vlee/Projects/MONSTAR2/MONSTAR_genome_data_20240617/")
for (e in bio_control){
  sample <- read_delim(e, comment = "#", col_names = FALSE, show_col_types = FALSE)%>%
    filter(X1 %in% genelist)%>%
    mutate(log2tpm = log2(X2+1))%>%
    select(c(X1, log2tpm))
  colnames(sample) <- c("X1", substr(e, 1, 11))
  biocontrol_datamatrix <- merge(biocontrol_datamatrix, sample, by = "X1", all.x=TRUE)
  
}

datamatrix <- biocontrol_datamatrix[,-1]
datamatrix <- as.matrix(datamatrix)
rownames(datamatrix) <- genelist


annot_bio <- key%>%
  filter(filename %in% bio_control)%>%
  select(Tumor_Sample_Barcode, CANCERTYPE)

anotcol_bio = list(Cancer = c(`Breast cancer` =  "hotpink2",
                                   `Colorectal cancer` = "green3"))

## GSVA ssGSEA method, hallmark and progeny_geneset

data_param <-ssgseaParam(datamatrix, hallmark) 
set.seed(178)
hallmark_res <- gsva(data_param)
sample_order <- match(colnames(hallmark_res), annot_bio$Tumor_Sample_Barcode)

annot_bio1 <- annot_bio[sample_order,]
rownames(annot_bio1) <- annot_bio1$Tumor_Sample_Barcode

annot_bio1 <- as.data.frame(annot_bio1[,-1])
colnames(annot_bio1) <- "Cancer"
set.seed(278)
pheatmap(t(hallmark_res),
         annotation_row = annot_bio1,
         show_rownames = FALSE, 
         show_colnames = TRUE,
         annotation_colors = anotcol_bio)


## Check ssGSEA result

res.pca <- prcomp(t(hallmark_res), scale. = TRUE )
fviz_eig(res.pca)
pca <- as.data.frame(res.pca$x)%>%
  select(c(PC1, PC2))
pca$Tumor_Sample_Barcode <- rownames(pca)
pca <- merge(pca, key, by ="Tumor_Sample_Barcode", all.x = TRUE)
nes_pca <- ggplot(pca, aes(x = PC1, y = PC2, colour =Version))+geom_point()+
  theme_light()+
  ggtitle("BCvsCRC PCA plot of all ssGSEA NES")+
  theme(axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 22))
nes_pca


## Perform t test for BC vs CRC

library(reshape2)

hallmark_res_df <- melt(t(hallmark_res))
colnames(hallmark_res_df) <- c("Tumor_Sample_Barcode", "Pathway", "NES")
pathway <- unique(hallmark_res_df$Pathway)
##annotate cancer type 

hallmark_res_df <- merge(hallmark_res_df, annot_bio, by ="Tumor_Sample_Barcode", all.x = TRUE )
#hallmark_res_df$CANCERTYPE <- as.factor(hallmark_res_df$CANCERTYPE)


##generate ttest results for BC vs CRC
ttest_result_bc_crc <- c()

for (p in pathway){
  df <- hallmark_res_df%>%
    filter(Pathway ==p)
  test <- t.test(df$NES~df$CANCERTYPE)
  
  ttest_result_bc_crc <- rbind(ttest_result_bc_crc, cbind(p, test$statistic, test$p.value, test$estimate[1], test$estimate[2]))
  
}
ttest_result_bc_crc <- as.data.frame(ttest_result_bc_crc)
colnames(ttest_result_bc_crc) <- c("Pathway", "statistic", "pval", "mean_BC", "mean_CRC")
ttest_result_bc_crc$statistic <- as.numeric(ttest_result_bc_crc$statistic)
ttest_result_bc_crc$pval <- as.numeric(ttest_result_bc_crc$pval)
ttest_result_bc_crc$mean_BC <- as.numeric(ttest_result_bc_crc$mean_BC)
ttest_result_bc_crc$mean_CRC <- as.numeric(ttest_result_bc_crc$mean_CRC)

bccrc_boxplot_pathways <- c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_WNT_BETA_CATENIN_SIGNALING")

ttest_result_bc_crc <- ttest_result_bc_crc%>%
  mutate(NES_delta =mean_BC-mean_CRC)%>%
  mutate(log2FC = log2(mean_BC/mean_CRC))%>%
  mutate(`-log10pval` = -log10(pval))%>%
  mutate(NES_Score= if_else(NES_delta>0.05 &pval<0.05, "UP",
                            if_else(NES_delta< -0.05&pval<0.05, "DOWN", "N.S")))%>%
  mutate(label = if_else(Pathway %in% bccrc_boxplot_pathways, Pathway, ""))


bc_crc_dotplot <- ggplot(ttest_result_bc_crc, aes(x = mean_BC, y = mean_CRC, label = label))+
  geom_point(aes(size = `-log10pval`, colour = log2FC ))+
  geom_text_repel()+
  scale_color_gradient(low = "green3", high = "hotpink1")+
  scale_size(range = c(0.1, 7))+
  theme_light()+
  geom_abline(intercept = 0, slope = 1, color = "black")+
  xlab("Mean NES score, breast cancers")+
  ylab("Mean NES score, colorectal cancers")+
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13))

bc_crc_dotplot

bc_crc_volcplot <- ggplot(ttest_result_bc_crc, aes(x = NES_delta, y = `-log10pval`, col = NES_Score, label= label ))+
  geom_point(size = 4)+
  geom_text_repel(size = 5)+
  scale_color_manual(values = c("green3", "lightgrey", "hotpink1"))+
  geom_hline(yintercept = -log10(0.05), col = "red")+
  geom_vline(xintercept = c(-0.05, 0.05), col= "red")+
  theme_light()+
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13))

bc_crc_volcplot


bc_crc_boxplotdf <- hallmark_res_df%>%
  filter(Pathway %in% bccrc_boxplot_pathways)


bc_crc_boxplot <- ggplot(bc_crc_boxplotdf, aes(x = Pathway, y = NES, fill = CANCERTYPE))+
  geom_boxplot()+
  theme_light()+
  scale_fill_manual(values = c("hotpink1", "green3"))+
  theme(axis.title = element_text(size = 18), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13))
bc_crc_boxplot
