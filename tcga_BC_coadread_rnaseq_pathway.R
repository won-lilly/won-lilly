## This script performs the ssGSEA analysis using TCGA BRCA and COADREAD RNAseq data
##  ssGSEA NES scores for hallmark pathways were generated for 
## BRCA and COADREAD patients, and comparisons were made between the 2 cancer types.

library(plyr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(tidyverse)
library(plotrix)
library(progeny)
library(SummarizedExperiment)
library(progeny)
library(reshape2)
library(ggrepel)


setwd("/fsx/home/vlee/Projects/MONSTAR2/EO_vs_LO_CRC/")
genelist <- unlist(read.csv("/fsx/home/vlee/Projects/MONSTAR2/Analysis_dataset/WTS_common_genes.csv"))

map <- read.delim("/fsx/home/xrao/resource/gene_annotation.txt")
map <- map%>%
  separate(ensembl, into = c("ensembl", "throw"), sep = "[.]")
  
map <- unique(map[,c("ensembl", "gene")])

coadread <- c("COAD", "READ")
tcga_clinical <- read.delim("coadread_tcga_pan_can_atlas_2018_clinical_data.tsv")

##read in TCGA data 
tcga_coad <- readRDS("/fsx/home/pjonsson/Data/TCGA/recount-rse-Jan2024.rds")$COAD
tcga_read <- readRDS("/fsx/home/pjonsson/Data/TCGA/recount-rse-Jan2024.rds")$READ

late <- c("STAGE III", "STAGE IIIA", "STAGE IIIB", "STAGE IIIC", "STAGE IV", "STAGE IVA", "STAGE IVB")
#Get late stage mss EO and LO sample IDs
eo_sampleID <- tcga_clinical%>%
  filter(Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code %in% late)%>%
  filter(MSIsensor.Score<3)%>%
  filter(Diagnosis.Age<50)%>%
  select(Sample.ID)%>%
  unlist()

lo_sampleID <- tcga_clinical%>%
  filter(Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code %in% late)%>%
  filter(MSIsensor.Score<3)%>%
  filter(Diagnosis.Age >=50)%>%
  select(Sample.ID)%>%
  unlist()
## get key from object
coad_key <- colData(tcga_coad)%>%
  as_tibble(rownames = "sample")%>%
  select(sample, tcga.tcga_barcode)%>%
  mutate(sampleID = substr(tcga.tcga_barcode, 1, 15))

read_key <- colData(tcga_read)%>%
  as_tibble(rownames = "sample")%>%
  select(sample, tcga.tcga_barcode)%>%
  mutate(sampleID = substr(tcga.tcga_barcode, 1, 15))

coadread_key <- rbind(coad_key, read_key)

##Generate list of samples from the key
eo_rnaseq_list <- coadread_key%>%
  filter(sampleID %in% eo_sampleID)%>%
  select(sample)%>%
  unlist()

lo_rnaseq_list <- coadread_key%>%
  filter(sampleID %in% lo_sampleID)%>%
  select(sample)%>%
  unlist()
all_tcgacoadread <- c(eo_rnaseq_list, lo_rnaseq_list)

##Create data matric with samples as columns and genes as rows
data_coad <- tcga_coad@assays@data$TPM
data_coad <- as.data.frame(data_coad)
data_coad$ensembl <- rownames(data_coad)
data_coad <- data_coad%>%
  separate(ensembl, into = c("ensembl", "throw"), sep = "[.]")

datacoad1 <- merge(data_coad, map, by = "ensembl", all.x = TRUE)
datacoad1<- datacoad1%>%
  filter(gene %in% genelist)
datamatrix_coad <- datacoad1%>%
  select(-c(ensembl, throw, gene))
datamatrix_coad <- as.matrix(datamatrix_coad)
rownames(datamatrix_coad) <- datacoad1$gene
datamatrix_coad <- log2(datamatrix_coad+1)

data_read <- tcga_read@assays@data$TPM
data_read <- as.data.frame(data_read)
data_read$ensembl <- rownames(data_read)
data_read <- data_read%>%
  separate(ensembl, into = c("ensembl", "throw"), sep = "[.]")

dataread1 <- merge(data_read, map, by = "ensembl", all.x = TRUE)
dataread1<- dataread1%>%
  filter(gene %in% genelist)
datamatrix_read <- dataread1%>%
  select(-c(ensembl, throw, gene))
datamatrix_read <- as.matrix(datamatrix_read)
rownames(datamatrix_read) <- dataread1$gene
datamatrix_read <- log2(datamatrix_read+1)


###Get breast cancer data
tcga_brca <- readRDS("/fsx/home/pjonsson/Data/TCGA/recount-rse-Jan2024.rds")$BRCA
breast_clinical <- read.delim("brca_tcga_pub2015_clinical_data.tsv")

breast_late <- c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV")
breast_list <- breast_clinical%>%
  filter(Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code %in% breast_late)%>%
  select(Sample.ID)%>%
  unique()%>%
  unlist()

##Prepare BRCA into a data matrix with samples as columns and genes as rows
data_brca <- tcga_brca@assays@data$TPM
data_brca <- as.data.frame(data_brca)
data_brca$ensembl <- rownames(data_brca)
data_brca <- data_brca%>%
  separate(ensembl, into = c("ensembl", "throw"), sep = "[.]")

data_brca1 <- merge(data_brca, map, by = "ensembl", all.x = TRUE)
data_brca1<- data_brca1%>%
  filter(gene %in% genelist)
datamatrix_brca <- data_brca1%>%
  select(-c(ensembl, throw, gene))
datamatrix_brca <- as.matrix(datamatrix_brca)
rownames(datamatrix_brca) <- data_brca1$gene
datamatrix_brca <- log2(datamatrix_brca+1)

##Get BRCA sample key and list
brca_key <- colData(tcga_brca)%>%
  as_tibble(rownames = "sample")%>%
  select(sample, tcga.tcga_barcode)%>%
  mutate(sampleID = substr(tcga.tcga_barcode, 1, 15))%>%
  filter(sampleID %in% breast_list)

brca_sample_list <- brca_key%>%
  select(sample)%>%
  unlist()

##Combine data matrices
## 2 comparisons were made:
## 1. BRCA versus COADREAD (this was performed as a sanity check for the ssGSEA analysis strategy)
## 2. COADREAD Early onset versus Late onset

## Extract list of common genes across the 3 matrices
commongenes <- intersect(rownames(datamatrix_coad), rownames(datamatrix_read))
commongenes <- intersect(commongenes, rownames(datamatrix_brca))
datamatrix_coad1 <- datamatrix_coad[commongenes,]
datamatrix_read1 <- datamatrix_read[commongenes,]
datamatrix_brca1 <- datamatrix_brca[commongenes, brca_sample_list]

## Combine COAD and READ for EO vs LO analyses
tcga_coadread <- cbind(datamatrix_coad1, datamatrix_read1)
tcga_coadread1 <- tcga_coadread[,all_tcgacoadread]
tcga_coadreadbrca <- cbind(tcga_coadread1, datamatrix_brca1)

## Prepare clinical data frames for using in heatmap
coadread_key1 <- coadread_key%>%
  filter(sample %in%all_tcgacoadread)%>%
  mutate(Onset = if_else(sample %in% eo_rnaseq_list, "EO", "LO"))%>%
  select(sample, Onset)


coadread_key_cancertype <- coadread_key1%>%
  mutate(Cancertype = "COADREAD")%>%
  select(-Onset)
brca_key_cancertype <- brca_key%>%
  mutate(Cancertype = "BRCA")%>%
  select(sample, Cancertype)
cancertypes <- rbind(coadread_key_cancertype, brca_key_cancertype)

##Generate ssGSEA NES result for EO vs LO

library(GSVA)
library(fgsea)
hallmark <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt")

data_param <-ssgseaParam(tcga_coadread1, hallmark) 
set.seed(178)
hallmark_res <- gsva(data_param)

## Match clinical row order to ssGSEA result df
order <- match(colnames(hallmark_res), coadread_key1$sample)
coadread_key2 <- coadread_key1[order, ]

##Create heatmaps Early vs Late
anno_eolo <- coadread_key2[,-1]
rownames(anno_eolo) <- coadread_key2$sample
anotcol_onset = list(Onset = c( EO=  "firebrick",
                                LO = "dodgerblue"))
pheatmap(t(hallmark_res),
         annotation_row = anno_eolo,
         show_rownames = FALSE, 
         annotation_colors = anotcol_onset)




##Create heatmaps Early vs Late
anno_eolo <- coadread_key2[,-1]
rownames(anno_eolo) <- coadread_key2$sample
anotcol_onset = list(Onset = c( EO=  "firebrick",
                                LO = "dodgerblue"))
pheatmap(t(hallmark_res),
         annotation_row = anno_eolo,
         show_rownames = FALSE, 
         annotation_colors = anotcol_onset)

##############################  T test ######################

##generate ttest results for EO vs LO
hallmark_res_df <- melt(t(hallmark_res))
colnames(hallmark_res_df) <- c("sample", "Pathway", "NES")
hallmark_res_df <- hallmark_res_df%>%
  mutate(onset = if_else(sample %in% eo_rnaseq_list, "EO", "LO"))
ttest_result_eolo <- c()

pathway <- unique(hallmark_res_df$Pathway)

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

eolo_lab <- c("HALLMARK_P53_PATHWAY", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
              "HALLMARK_IL6_JAK_STAT3_SIGNALING")
ttest_result_eolo <- ttest_result_eolo%>%
  mutate(NESdelta =mean_EO-mean_LO)%>%
  mutate(`-log10pval` = -log10(pval))%>%
  mutate(NES_Score= if_else(NESdelta>0.1 &pval<0.1, "UP",
                            if_else(NESdelta< -0.1&pval<0.1, "DOWN", "N.S")))%>%
  mutate(label = if_else(Pathway%in% eolo_lab, Pathway, ""))


write.csv(ttest_result_eolo, "TCGA_COADREAD_nes_eo_lo_ttest.csv", row.names = TRUE)

eolo_volcplot <- ggplot(ttest_result_eolo, aes(x = NESdelta, y = `-log10pval`, col = NES_Score, label= label))+
  geom_point(size = 4)+
  geom_text_repel()+
  scale_color_manual(values = c("black", "lightgrey", "hotpink1"))+
  geom_hline(yintercept = -log10(0.05), col = "red")+
  geom_vline(xintercept = c(-0.1, 0.1), col= "red")+
  theme_light()+
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13))+
  xlim(c(-0.02, 0.02))+
  ylim(c(0,1.1))

eolo_volcplot


## Repeated the same steps to compare BC versus COAREAD

data_param <-ssgseaParam(tcga_coadreadbrca, hallmark) 
set.seed(178)
hallmark_res <- gsva(data_param)

## Match row order to ssGSEA result df
order <- match(colnames(hallmark_res), cancertypes$sample)
cancertypes2 <-cancertypes[order,] 

##Plot heatmap for results
anno_bio <- cancertypes2[,-1]
anotcol_bio = list(Cancertype = c(BRCA =  "hotpink2",
                                  COADREAD = "green3"))
pheatmap(t(hallmark_res),
         annotation_row = anno_bio,
         show_rownames = FALSE, 
         annotation_colors = anotcol_bio)



#############################################################################
##generate ttest results for BRCA vs COADREAD

hallmark_res_df <- melt(t(hallmark_res))
colnames(hallmark_res_df) <- c("sample", "Pathway", "NES")
hallmark_res_df <- hallmark_res_df%>%
  mutate(onset = if_else(sample %in% all_tcgacoadread, "COADREAD", "BRCA"))

ttest_bc_crc <- c()
pathway <- unique(hallmark_res_df$Pathway)

for (p in pathway){
  df <- hallmark_res_df%>%
    filter(Pathway ==p)
  test <- t.test(df$NES~df$onset)
  
  ttest_bc_crc <- rbind(ttest_bc_crc, cbind(p, test$statistic, test$p.value, test$estimate[1], test$estimate[2]))
  
}
ttest_bc_crc <- as.data.frame(ttest_bc_crc)
colnames(ttest_bc_crc) <- c("Pathway", "statistic", "pval", "mean_BRCA", "mean_COADREAD")
ttest_bc_crc$statistic <- as.numeric(ttest_bc_crc$statistic)
ttest_bc_crc$pval <- as.numeric(ttest_bc_crc$pval)
ttest_bc_crc$mean_BRCA <- as.numeric(ttest_bc_crc$mean_BRCA)
ttest_bc_crc$mean_COADREAD <- as.numeric(ttest_bc_crc$mean_COADREAD)

pathlabs <- c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_WNT_BETA_CATENIN_SIGNALING")
ttest_bc_crc <- ttest_bc_crc%>%
  mutate(NESdelta =mean_BRCA-mean_COADREAD )%>%
  mutate(`-log10pval` = -log10(pval))%>%
  mutate(NES_status= if_else(NESdelta>0.05 &pval<0.05, "UP",
                            if_else(NESdelta< -0.05&pval<0.05, "DOWN", "N.S")))%>%
  mutate(label = if_else(Pathway %in% pathlabs, Pathway, ""))


bc_crc_volcplot <- ggplot(ttest_bc_crc, aes(x = NESdelta, y = `-log10pval`, col = NES_status, label= label ))+
  geom_point(size = 4)+
  geom_text_repel(size = 7)+
  scale_color_manual(values = c("green3", "lightgrey", "hotpink1"))+
  geom_hline(yintercept = -log10(0.05), col = "red")+
  geom_vline(xintercept = c(-0.05, 0.05), col= "red")+
  theme_light()+
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13))

bc_crc_volcplot
