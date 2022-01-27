###############
### Download GSE DATA
GSE113182 <- read.table("GSE113182_raw_count.txt", header = T, sep='\t') #readcount data
head(GSE113182)
dim(GSE113182)

### Create metadata for GSE data
ID <- c("CMP_ID236","CMP_ID235", "CMP_ID185", "GMP_ID236","GMP_ID235","GMP_ID185",
	"HSC_ID236","HSC_ID235","HSC_ID185","MEP_ID236","MEP_ID235","MEP_ID185",
	"CMP_ID132","CMP_ID234","GMP_ID132","GMP_ID234","HSC_ID132", "HSC_ID234","MEP_ID132","MEP_ID234")

index <- c(rep(c("common_myeloid_progenitors","granulocyte-macrophage_progenitors",
	"hematopoietic_stem_cells","megakaryocyte-erythrocyte_progenitors"),c(3,3,3,3)), rep(c(
	"common_myeloid_progenitors","granulocyte-macrophage_progenitors","hematopoietic_stem_cells",
	"megakaryocyte-erythrocyte_progenitors"),c(2,2,2,2)))

metadata <- data.frame(index)
rownames(metadata) <- ID
###############


###############
### Pre-Filtering

## protein-coding-gene filtering
library(dplyr)

# Upload gmt file(whole gene info)
ann = rtracklayer::import("/spstorage/DB/NGS_reference/GRCh38_v34/gencode.v34.primary_assembly.annotation.gtf") %>% as.data.frame()
pcg = ann %>% filter(gene_type == "protein_coding") %>% distinct(gene_id)
pcg$gene_id <- gsub("\\.[0-9]+","",pcg$gene_id)

# Applying to GSEdata
pcglist <- ifelse(GSE113182$Feature %in% pcg$gene_id, "True", "False")
GSE113182_pcg <- GSE113182[which(pcglist == 'True'),]
head(GSE113182_pcg)
dim(GSE113182_pcg)

## gene-expression-average
## Finding data with an average expression of 5 or higher
# HSC
GSE113182_pcg %>% select(contains("HSC")) %>% head
GSE113182_HSC <- GSE113182_pcg %>% select(Feature.ID, contains("HSC")) %>%
				mutate(rowmean = (HSC_ID236 + HSC_ID235 + HSC_ID185 + HSC_ID132 + HSC_ID234)/5) %>%
				filter(rowmean >= 5) %>% select(Feature.ID:HSC_ID234)
# CMP
GSE113182_pcg %>% select(contains("CMP")) %>% head
GSE113182_CMP <- GSE113182_pcg %>% select(Feature.ID, contains("CMP")) %>%
				mutate(rowmean = (CMP_ID236 + CMP_ID235 + CMP_ID185 + CMP_ID132 + CMP_ID234)/5) %>%
	 			filter(rowmean >= 5) %>% select(Feature.ID:CMP_ID234)
# GMP
GSE113182_pcg %>% select(contains("GMP")) %>% head
GSE113182_GMP <- GSE113182_pcg %>% select(Feature.ID, contains("GMP")) %>%
				mutate(rowmean = (GMP_ID236 + GMP_ID235 + GMP_ID185 + GMP_ID132 + GMP_ID234)/5) %>%
				filter(rowmean >= 5) %>% select(Feature.ID:GMP_ID234)
# MEP
GSE113182_pcg %>% select(contains("MEP")) %>% head
GSE113182_MEP <- GSE113182_pcg %>% select(Feature.ID, contains("MEP")) %>%
				mutate(rowmean = (MEP_ID236 + MEP_ID235 + MEP_ID185 + MEP_ID132 + MEP_ID234)/5) %>%
				filter(rowmean >= 5) %>% select(Feature.ID:MEP_ID234)
###############


###############
### Check DEG(differentially expressed gene)

## 1. Grouping by stage of differentiation
## HSC -> CMP
HSC_CMP <- GSE113182_HSC %>% inner_join(GSE113182_CMP, by = "Feature.ID")
## CMP -> GMP
CMP_GMP <- GSE113182_CMP %>% inner_join(GSE113182_GMP, by = "Feature.ID")
## CMP -> MEP
CMP_MEP <- GSE113182_CMP %>% inner_join(GSE113182_MEP, by = "Feature.ID")


## 2. Create metadata in different stages
## HSC->CMP: metadata_1
ID_1 <- c("HSC_ID236","HSC_ID235","HSC_ID185","HSC_ID132", "HSC_ID234",
	"CMP_ID236","CMP_ID235", "CMP_ID185","CMP_ID132","CMP_ID234")
index_1 <- rep(c("hematopoietic_stem_cells","common_myeloid_progenitors"),each=5)

metadata_1 <- data.frame(index_1)
rownames(metadata_1) <- ID_1

## CMP->GMP: metadata_2
ID_2 <- c("CMP_ID236", "CMP_ID235", "CMP_ID185", "CMP_ID132", "CMP_ID234", "GMP_ID236",
	 "GMP_ID235", "GMP_ID185", "GMP_ID132", "GMP_ID234")
index_2 <- rep(c("common_myeloid_progenitors","granulocyte-macrophage_progenitors"),each=5)

metadata_2 <- data.frame(index_2)
rownames(metadata_2) <- ID_2

## CMP->MEP: metadata_3
ID_3 <- c("CMP_ID236", "CMP_ID235", "CMP_ID185", "CMP_ID132", "CMP_ID234", "MEP_ID236",
	"MEP_ID235", "MEP_ID185", "MEP_ID132", "MEP_ID234")
index_3 <- rep(c("common_myeloid_progenitors","megakaryocyte-erythrocyte_progenitors"),each=5)

metadata_3 <- data.frame(index_3)
rownames(metadata_3) <- ID_3


## 3. Create DESeqDataSet
library(DESeq2)

## HSC -> CMP
dds1 <- DESeqDataSetFromMatrix(countData = HSC_CMP ,
                               colData = metadata_1,
                               design = ~index_1, tidy = TRUE)
dds_HSC_CMP <- DESeq(dds1)
res_HSC_CMP <- results(dds_HSC_CMP)
head(res_HSC_CMP)

# Check dimension
#res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% dim  #887
res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1.5, padj < 0.05) %>% dim  #457
#res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 2, padj < 0.05) %>% dim  #296
#res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 2.5, padj < 0.05) %>% dim  #204

## CMP -> GMP
dds2 <- DESeqDataSetFromMatrix(countData = CMP_GMP ,
                               colData = metadata_2,
                               design = ~index_2, tidy = TRUE)
dds_CMP_GMP <- DESeq(dds2)
res_CMP_GMP <- results(dds_CMP_GMP)
head(res_CMP_GMP)

# Check dimension
res_CMP_GMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% dim   #334개 
#res_CMP_GMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1.1, padj < 0.05) %>% dim  #297개
#res_CMP_GMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1.5, padj < 0.05) %>% dim  #213개

## CMP -> MEP
dds3 <- DESeqDataSetFromMatrix(countData = CMP_MEP ,
                               colData = metadata_3,
                               design = ~index_3, tidy = TRUE)
dds_CMP_MEP <- DESeq(dds3)
res_CMP_MEP <- results(dds_CMP_MEP)
head(res_CMP_MEP)

# Check dimension
res_CMP_MEP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% dim  #201개
###############


###############
### MA-Plot of deseq data
## HSC -> CMP
plotMA(res_HSC_CMP)
plotMA(res_HSC_CMP, ylim=c(-2,2))

## CMP -> GMP
plotMA(res_CMP_GMP)
plotMA(res_CMP_GMP, ylim=c(-2,2))

## CMP -> MEP 
plotMA(res_CMP_MEP)
plotMA(res_CMP_MEP, ylim=c(-2,2))
###############


###############
### Create TPM dataset
## 1. Prepare data containing gene length
library(tidyverse)

# Upload data
load('/spstorage/USERS/sung/projects/DEGanalysis/Genelength.RData')
head(genelength)

genelength_1 <- rownames_to_column(genelength, var='Feature.ID')
head(genelength_1)

# Combine length columns with existing pcg data
GSE113182_length <- GSE113182_pcg %>% left_join(genelength_1, by = "Feature.ID")

# Make gene data and number data, respectively
GSE113182_genes <- GSE113182_length[,c(1,length(GSE113182_length))]
GSE113182_counts <- GSE113182_length[,-c(1,length(GSE113182_length))]


## 2. Create a function that calculates TPM
tpm <- function(counts, lengths){
	rate <- counts / lengths
	rate / sum(rate) * 1e6
}


## 3. Create TPM dataframe
#Apply to created datas(genes, counts)
GSE113182_tpm_counts <- tpm(GSE113182_counts, GSE113182_genes$length)
head(GSE113182_tpm_counts)

GSE113182_tpm <- data.frame(GSE113182_genes$Feature.ID, GSE113182_tpm_counts)
colnames(GSE113182_tpm)[1] <- 'Feature.ID'


## 4. Filtering by gene-expression-average
GSE113182_tpmcut <- GSE113182_tpm %>% mutate(rowmeans = rowMeans(GSE113182_tpm[,-1])) %>%
	filter(rowmeans > 2) %>% select(-rowmeans)
dim(GSE113182_tpmcut)
###############


###############
### PCA-Plot of TPM data
## 1. PCA(Principal component analysis)
# Data trimming
GSE113182_tpmcut_trim <- GSE113182_tpmcut[,-1] %>% t() 

# PCA
pca_tpmcut <- prcomp(GSE113182_tpmcut_trim, scale. = TRUE)
summary(pca_tpmcut)

pca_x_df1 <- pca_tpmcut$x %>% as.data.frame() #dataframing only x
pca_x_df2 <- rownames_to_column(pca_x_df1, var='index')
rownames(pca_x_df2) <- rownames(pca_x_df1)

# organize PC data (x)
pca_x_df <- pca_x_df2 %>% separate(index, into = c("dex","IDnum"), sep="_") %>% select(-IDnum)


## 2. PCA-Plot
pca_x_df %>% ggplot(aes(x = PC1, y = PC2, color = dex)) +
	geom_point() +
	labs(x = "PC1: % variance", y = "PC2: % variance",
		title = "PCA-plot", subtitle = "Whole stage(HSC -> CMP -> GMP/MEP)")


## * Excluding HSC(hematopoietic stem cells)
## 1. PCA(Principal component analysis)
# Create index
index_excluding_hsc <- c(rep(c("common_myeloid_progenitors","granulocyte-macrophage_progenitors",
	"megakaryocyte-erythrocyte_progenitors"),c(3,3,3)), rep(c(
	"common_myeloid_progenitors","granulocyte-macrophage_progenitors",
	"megakaryocyte-erythrocyte_progenitors"),c(2,2,2)))

# Data trimming
GSE113182_excluding_hsc <- GSE113182_tpmcut[,-c(1,8:10,18:19)] %>% t()

# PCA
pca_excluding_hsc <- prcomp(GSE113182_excluding_hsc, scale. = TRUE)
summary(pca_excluding_hsc)

pca_x_df_excluding_1 <- pca_excluding_hsc$x %>% as.data.frame()
pca_x_df_excluding_2 <- rownames_to_column(pca_x_df_excluding_1, var='index')
rownames(pca_x_df_excluding_2) <- rownames(pca_x_df_excluding_1)

# organize PC data
pca_x_excluding_hsc <- pca_x_df_excluding_2 %>%
	separate(index, into = c("dex","IDnum"), sep="_") %>% select(-IDnum)

## 2. PCA-Plot
pca_x_excluding_hsc %>% ggplot(aes(x = PC1, y = PC2, color = dex)) +
	geom_point() +
	labs(x = "PC1: 45% variance", y = "PC2: 14% variance",
		title = "PCA-plot", subtitle = "except HSC stage")
###############


###############
### Visualizing data using t-SNE
## 1. t-SNE(t-distributed stochastic neighbor embedding)
library(Rtsne)
set.seed(2022)

# Data trimming
GSE113182_whole <- GSE113182_tpmcut %>% select(-Feature.ID) %>% t()

# t-SNE(perplexity: 2~6)
tsne_whole <- Rtsne(GSE113182_whole, perplexity=2, check_duplicates=FALSE)

# organize t-SNE data (Y)
tsne_whole_df <- as.data.frame(tsne_whole$Y)
tsne_whole_df <- data.frame(tsne_whole_df, index) #add index


## 2. Visualize t-SNE Data
tsne_whole_df %>% ggplot(aes(x = V1, y = V2, color=index)) +
	geom_point() +
	ggtitle('t-SNE-Plot', subtitle = 'Whole stage(HSC -> CMP -> GMP/MEP)')


## * Excluding HSC(hematopoietic stem cells)
## 1. t-SNE(*t-*distributed stochastic neighbor embedding)
# Data trimming
GSE113182_excluding_hsc <- GSE113182_tpmcut %>% select(-Feature.ID) %>% select(-contains('HSC')) %>% t()

# t-SNE(perplexity: 2~4)
tsne_excluding_hsc <- Rtsne(GSE113182_excluding_hsc, perplexity=2, check_duplicates=FALSE)

# organize t-SNE data (Y)
tsne_excluding_hsc_df <- as.data.frame(tsne_excluding_hsc$Y)
tsne_excluding_hsc_df <- data.frame(tsne_excluding_hsc_df, index_excluding_hsc) #add index


## 2. Visualize t-SNE Data
tsne_excluding_hsc_df %>% ggplot(aes(x = V1, y = V2, color=index_excluding_hsc)) +
	geom_point() +
	ggtitle('t-SNE-Plot', subtitle = 'except HSC stage')
###############


###############
### Select DEG(differentially expressed gene)
## Extracting data for step
RES_HSC_CMP <- res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1.7, padj < 0.05)
RES_CMP_GMP <- res_CMP_GMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05)
RES_CMP_MEP <- res_CMP_MEP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05)





############
RES_HSC_CMP <- RES_HSC_CMP %>% mutate(Feature.ID = rownames(RES_HSC_CMP)) %>%
	select(Feature.ID, baseMean:padj)
RES_CMP_GMP <- RES_CMP_GMP %>% mutate(Feature.ID = rownames(RES_CMP_GMP)) %>%
	select(Feature.ID, baseMean:padj)

RES_common_1 <- RES_HSC_CMP %>% inner_join(RES_CMP_GMP, by = "Feature.ID") #HSC->CMP->GMP DEG
dim(RES_common_1) #56개

RES_CMP_MEP <- RES_CMP_MEP %>% mutate(Feature.ID = rownames(RES_CMP_MEP)) %>%
	select(Feature.ID, baseMean:padj)

RES_common_2 <- RES_HSC_CMP %>% inner_join(RES_CMP_MEP, by = "Feature.ID") #HSC->CMP->MEP DEG
dim(RES_common_2) #40개

RES_common_3 <- RES_common_1 %>% inner_join(RES_common_2, by = "Feature.ID") #모든 과정 공통 DEG
dim(RES_common_3) #23개
RES_common_3$Feature.ID

RES_HSC_CMP_name <- RES_HSC_CMP$Feature.ID
RES_CMP_GMP_name <- RES_CMP_GMP$Feature.ID

RES_common_1_name <- list(HSC_CMP = RES_HSC_CMP_name, CMP_GMP = RES_CMP_GMP_name)

group.venn(list(HSC_CMP = RES_HSC_CMP_name, CMP_GMP = RES_CMP_GMP_name), label = TRUE,
	lab.cex = 0.5, width = 20, height = 20)


#상위,하위 DEG 200개씩 선별
# HSC -> CMP
updeg_HSC_CMP <- res_HSC_CMP %>% as.data.frame %>% filter(padj < 0.05) %>%
	arrange(desc(log2FoldChange)) %>% head(200)
downdeg_HSC_CMP <- res_HSC_CMP %>% as.data.frame %>% filter(padj < 0.05) %>% 
	arrange(log2FoldChange) %>% head(200)
deg_HSC_CMP <- rbind(updeg_HSC_CMP, downdeg_HSC_CMP)

# CMP -> GMP
updeg_CMP_GMP <- res_CMP_GMP %>% as.data.frame %>% filter(padj < 0.05) %>%
	arrange(desc(log2FoldChange)) %>% head(200)
downdeg_CMP_GMP <- res_CMP_GMP %>% as.data.frame %>% filter(padj < 0.05) %>% 
	arrange(log2FoldChange) %>% head(200)
deg_CMP_GMP <- rbind(updeg_CMP_GMP, downdeg_CMP_GMP)

# CMP -> MEP
updeg_CMP_MEP <- res_CMP_MEP %>% as.data.frame %>% filter(padj < 0.05) %>%
	arrange(desc(log2FoldChange)) %>% head(200)
downdeg_CMP_MEP <- res_CMP_MEP %>% as.data.frame %>% filter(padj < 0.05) %>% 
	arrange(log2FoldChange) %>% head(200)
deg_CMP_MEP <- rbind(updeg_CMP_MEP, downdeg_CMP_MEP)


DEG_HSC_CMP <- deg_HSC_CMP %>% mutate(Feature.ID = rownames(deg_HSC_CMP)) %>%
	select(Feature.ID, baseMean:padj)
DEG_CMP_GMP <- deg_CMP_GMP %>% mutate(Feature.ID = rownames(deg_CMP_GMP)) %>%
	select(Feature.ID, baseMean:padj)
DEG_CMP_MEP <- deg_CMP_MEP %>% mutate(Feature.ID = rownames(deg_CMP_MEP)) %>%
	select(Feature.ID, baseMean:padj)

deg_common_1 <- DEG_HSC_CMP %>% inner_join(DEG_CMP_GMP, by = "Feature.ID") #HSC->CMP->GMP DEG
dim(deg_common_1) #60개
deg_common_2 <- DEG_HSC_CMP %>% inner_join(DEG_CMP_MEP, by = "Feature.ID") #HSC->CMP->MEP DEG
dim(deg_common_2) #44개
deg_common_3 <- deg_common_1 %>% inner_join(DEG_CMP_MEP, by = "Feature.ID") #모든 과정 공통 DEG
dim(deg_common_3) #25개

deg_common_3$Feature.ID
RES_common_3$Feature.ID

#UP, DOWN 구분안한 것과 구분한 것 중복 GENE 찾기
DEG_name <- c(RES_common_3$Feature.ID, deg_common_3$Feature.ID)
DEG_name[which(duplicated(DEG_name) == TRUE)]
length(DEG_name[which(duplicated(DEG_name) == TRUE)])


## log2FoldChange graph
log2FC <- c(deg_common_3$log2FoldChange.x, deg_common_3$log2FoldChange.y, deg_common_3$log2FoldChange)
index_stage <- rep(c("HSC->CMP", "CMP->GMP", "CMP->MEP"), c(25, 25, 25))
index_stage <- factor(index_stage, levels = c("HSC->CMP","CMP->GMP", "CMP->MEP"))
DEG_log2FC <- data.frame(deg_common_3$Feature.ID, log2FC, index_stage)

DEG_log2FC %>% ggplot(aes(x = index_stage, y = log2FC, group = deg_common_3.Feature.ID,
	color = deg_common_3.Feature.ID)) +
	geom_point() +
	geom_line() +
	scale_color_discrete(name = "Feature.ID") + 
	labs(x = 'Differentiation Stage', y = 'log2FoldChange',
		title = 'log2FoldChange of DEG', subtitle = 'by Differentiation stage')