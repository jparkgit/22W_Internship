############
## GSEdata 불러오기
GSE11382 <- read.table("GSE113182_raw_count.txt", header = T, sep='\t')
head(GSE11382)
dim(GSE11382)
############


############
## metadata 만들기
ID <- c("CMP_ID236","CMP_ID235", "CMP_ID185", "GMP_ID236","GMP_ID235","GMP_ID185",
	"HSC_ID236","HSC_ID235","HSC_ID185","MEP_ID236","MEP_ID235","MEP_ID185",
	"CMP_ID132","CMP_ID234","GMP_ID132","GMP_ID234","HSC_ID132", "HSC_ID234","MEP_ID132","MEP_ID234")

index <- c(rep(c("common_myeloid_progenitors","granulocyte-macrophage_progenitors",
	"hematopoietic_stem_cells","megakaryocyte-erythrocyte_progenitors"),c(3,3,3,3)), rep(c(
	"common_myeloid_progenitors","granulocyte-macrophage_progenitors","hematopoietic_stem_cells",
	"megakaryocyte-erythrocyte_progenitors"),c(2,2,2,2)))

metadata <- data.frame(index)
rownames(metadata) <- ID
############


############
## Pre-filtering (Protein coding gene filtering)
# 전체 유전자의 info가 담긴 gmt파일을 불러와 그중 pcg만 필터링
library(dplyr)

ann = rtracklayer::import("/spstorage/DB/NGS_reference/GRCh38_v34/gencode.v34.primary_assembly.annotation.gtf") %>% as.data.frame()
pcg = ann %>% filter(gene_type == "protein_coding") %>% distinct(gene_id)
pcg$gene_id <- gsub("\\.[0-9]+","",pcg$gene_id)

# GSE데이터에서 pcg 필터링
pcglist <- ifelse(GSE11382$Feature %in% pcg$gene_id, "True", "False")
GSE11382_pcg <- GSE11382[which(pcglist == 'True'),]
head(GSE11382_pcg)
dim(GSE11382_pcg)
############


############
## Pre-filtering (유전자 발현 합 평균 5 이상 Filtering)
#HSC: 유전자 발현 합이 평균 5이상 Filter
GSE11382_pcg %>% select(contains("HSC")) %>% head
GSE11382_HSC <- GSE11382_pcg %>% select(Feature.ID, contains("HSC")) %>%
				mutate(rowmean = (HSC_ID236 + HSC_ID235 + HSC_ID185 + HSC_ID132 + HSC_ID234)/5) %>%
				filter(rowmean >= 5) %>% select(Feature.ID:HSC_ID234)
#CMP: 유전자 발현 합이 평균 5이상 Filter
GSE11382_pcg %>% select(contains("CMP")) %>% head
GSE11382_CMP <- GSE11382_pcg %>% select(Feature.ID, contains("CMP")) %>%
				mutate(rowmean = (CMP_ID236 + CMP_ID235 + CMP_ID185 + CMP_ID132 + CMP_ID234)/5) %>%
	 			filter(rowmean >= 5) %>% select(Feature.ID:CMP_ID234)
#GMP: 유전자 발현 합이 평균 5이상 Filter
GSE11382_pcg %>% select(contains("GMP")) %>% head
GSE11382_GMP <- GSE11382_pcg %>% select(Feature.ID, contains("GMP")) %>%
				mutate(rowmean = (GMP_ID236 + GMP_ID235 + GMP_ID185 + GMP_ID132 + GMP_ID234)/5) %>%
				filter(rowmean >= 5) %>% select(Feature.ID:GMP_ID234)
#MEP: 유전자 발현 합이 평균 5이상 Filter
GSE11382_pcg %>% select(contains("MEP")) %>% head
GSE11382_MEP <- GSE11382_pcg %>% select(Feature.ID, contains("MEP")) %>%
				mutate(rowmean = (MEP_ID236 + MEP_ID235 + MEP_ID185 + MEP_ID132 + MEP_ID234)/5) %>%
				filter(rowmean >= 5) %>% select(Feature.ID:MEP_ID234)
############


############
## 각 분화단계 DEG 선별하기
## 두 데이터씩 합치기 
# HSC -> CMP 데이터만들기
HSC_CMP <- GSE11382_HSC %>% inner_join(GSE11382_CMP, by = "Feature.ID")
# CMP -> GMP 데이터만들기
CMP_GMP <- GSE11382_CMP %>% inner_join(GSE11382_GMP, by = "Feature.ID")
# CMP -> MEP 데이터만들기 
CMP_MEP <- GSE11382_CMP %>% inner_join(GSE11382_MEP, by = "Feature.ID")


# metadata_1 만들기 (HSC -> CMP)
ID_1 <- c("HSC_ID236","HSC_ID235","HSC_ID185","HSC_ID132", "HSC_ID234",
	"CMP_ID236","CMP_ID235", "CMP_ID185","CMP_ID132","CMP_ID234")
index_1 <- rep(c("hematopoietic_stem_cells","common_myeloid_progenitors"),each=5)

metadata_1 <- data.frame(index_1)
rownames(metadata_1) <- ID_1 

#DESeqDataSet 만들기
library(DESeq2)

dds1 <- DESeqDataSetFromMatrix(countData = HSC_CMP ,
                               colData = metadata_1,
                               design = ~index_1, tidy = TRUE)

dds_HSC_CMP <- DESeq(dds1)
res_HSC_CMP <- results(dds_HSC_CMP)
head(res_HSC_CMP)

#res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% dim  #887
res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1.5, padj < 0.05) %>% dim  #457
#res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 2, padj < 0.05) %>% dim  #296
#res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 2.5, padj < 0.05) %>% dim  #204

# RES_HSC_CMP <- results(dds_HSC_CMP, contrast=c("index_1", "hematopoietic_stem_cells",
# 	"common_myeloid_progenitors"))
RES_HSC_CMP <- res_HSC_CMP %>% as.data.frame() %>% filter(abs(log2FoldChange) > 1.5, padj < 0.05)

# metadata_2 만들기 (CMP -> GMP)
ID_2 <- c("CMP_ID236", "CMP_ID235", "CMP_ID185", "CMP_ID132", "CMP_ID234", "GMP_ID236",
	 "GMP_ID235", "GMP_ID185", "GMP_ID132", "GMP_ID234")
index_2 <- rep(c("common_myeloid_progenitors","granulocyte-macrophage_progenitors"),each=5)

metadata_2 <- data.frame(index_2)
rownames(metadata_2) <- ID_2 

#DESeqDataSet 만들기
dds2 <- DESeqDataSetFromMatrix(countData = CMP_GMP ,
                               colData = metadata_2,
                               design = ~index_2, tidy = TRUE)

dds_CMP_GMP <- DESeq(dds2)
res_CMP_GMP <- results(dds_CMP_GMP)
head(res_CMP_GMP)

res_CMP_GMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% dim   #334개 
res_CMP_GMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1.1, padj < 0.05) %>% dim  #297개
res_CMP_GMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1.5, padj < 0.05) %>% dim  #213개

#RES_CMP_GMP <- results(dds_CMP_GMP, contrast=c("index_2", "common_myeloid_progenitors", "granulocyte-macrophage_progenitors"))

RES_HSC_CMP <- res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 2.5, padj < 0.05)

# metadata_3 만들기 (CMP -> MEP)
ID_3 <- c("CMP_ID236", "CMP_ID235", "CMP_ID185", "CMP_ID132", "CMP_ID234", "MEP_ID236",
	"MEP_ID235", "MEP_ID185", "MEP_ID132", "MEP_ID234")
index_3 <- rep(c("common_myeloid_progenitors","megakaryocyte-erythrocyte_progenitors"),each=5)

metadata_3 <- data.frame(index_3)
rownames(metadata_3) <- ID_3

#DESeqDataSet 만들기
dds3 <- DESeqDataSetFromMatrix(countData = CMP_MEP ,
                               colData = metadata_3,
                               design = ~index_3, tidy = TRUE)

dds_CMP_MEP <- DESeq(dds3)
res_CMP_MEP <- results(dds_CMP_MEP)
head(res_CMP_MEP)

res_CMP_MEP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% dim  #201개 

RES_CMP_MEP <- res_CMP_MEP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05)
#RES_CMP_MEP <- results(dds_CMP_MEP, contrast=c("index_3", "common_myeloid_progenitors", "megakaryocyte-erythrocyte_progenitors"))
############


############
## MA-Plot 그리기
#HSC -> CMP
plotMA(res_HSC_CMP)
plotMA(res_HSC_CMP, ylim=c(-2,2))

# CMP -> GMP
plotMA(res_CMP_GMP)
plotMA(res_CMP_GMP, ylim=c(-2,2))

# CMP -> MEP 
plotMA(res_CMP_MEP)
plotMA(res_CMP_MEP, ylim=c(-2,2))

# resultsNames(dds_HSC_CMP) #coef 확인
# res_HSC_CMP_LFC <- lfcShrink(dds_HSC_CMP, 
# 	coef="index_1_hematopoietic_stem_cells_vs_common_myeloid_progenitors", type="apeglm")
# res_HSC_CMP_LFC


# ############
# ## PCA plot 그리기 
# #HSC -> CMP
# vsdata_1 <- vst(dds_HSC_CMP, blind = FALSE) 
# plotPCA(vsdata_1, intgroup="index_1")

# #CMP -> GMP
# vsdata_2 <- vst(dds_CMP_GMP, blind = FALSE) 
# plotPCA(vsdata_2, intgroup="index_2")

# #CMP -> MEP
# vsdata_3 <- vst(dds_CMP_MEP, blind = FALSE) 
# plotPCA(vsdata_3, intgroup="index_3")
############


############
## TPM 자료 만들기
library(tibble)
load('/spstorage/USERS/sung/projects/DEGanalysis/Genelength.RData') #유전자 size data upload
head(genelength)

head(GSE11382_pcg)

genelength_1 <- rownames_to_column(genelength, var='Feature.ID')
head(genelength_1)

#기존 pcg 자료에 length 열 합치기
GSE11382_length <- GSE11382_pcg %>% left_join(genelength_1, by = "Feature.ID")

GSE11382_genes <- GSE11382_length[,c(1,length(GSE11382_length))]
head(GSE11382_genes)
GSE11382_counts <- GSE11382_length[,-c(1,length(GSE11382_length))]
head(GSE11382_counts)

#tpm 구하는 함수
tpm <- function(counts, lengths){
	rate <- counts / lengths
	rate / sum(rate) * 1e6
}

GSE11382_tpm_counts <- tpm(GSE11382_counts, GSE11382_genes$length)
head(GSE11382_tpm_counts)

#tpm dataframe 만들기
GSE11382_tpm <- data.frame(GSE11382_genes$Feature.ID, GSE11382_tpm_counts)
colnames(GSE11382_tpm)[1] <- 'Feature.ID'
head(GSE11382_tpm)

#유전자 발현 평균 2 이상 filtering
GSE11382_tpmcut <- GSE11382_tpm %>% mutate(rowmeans = rowMeans(GSE11382_tpm[,-1])) %>%
	filter(rowmeans > 2) %>% select(-rowmeans)
head(GSE11382_tpmcut)
dim(GSE11382_tpmcut)
############


############
# PCA-Plot Plotting
library(tidyverse)
set.seed(2022)

index <- c(rep(c("common_myeloid_progenitors","granulocyte-macrophage_progenitors",
	"hematopoietic_stem_cells","megakaryocyte-erythrocyte_progenitors"),c(3,3,3,3)), rep(c(
	"common_myeloid_progenitors","granulocyte-macrophage_progenitors","hematopoietic_stem_cells",
	"megakaryocyte-erythrocyte_progenitors"),c(2,2,2,2)))
metadata <- data.frame(index)
rownames(metadata) <- colnames(GSE11382_tpm)[-1]

# GSE11382_tpmcut_sample <- GSE11382_tpmcut %>% sample_n(1000)
# tsne_tpm <- Rtsne(GSE11382_tpm_sample, perplexity=10, check_duplicates=FALSE)
# plot(tsne_tpm$Y, col = "black", bg = metadata$index, pch = 21, cex = 1)

# library(tsne)
# GSE11382_tpm_1 <- tsne(as.matrix(GSE11382_tpmcut[,-1]))
# GSE11382_tpm_1_1 <- GSE11382_tpm_1 %>% as.data.frame() %>% tbl_df() %>% mutate()


GSE11382_tpmcut_1 <- GSE11382_tpmcut[,-1] #Feature.ID 제거
GSE11382_tpmcut_2 = GSE11382_tpmcut_1 %>% t() #행, 열 바꾸기

pca_tpmcut <- prcomp(GSE11382_tpmcut_2, scale. = TRUE) #pca 구하기
summary(pca_tpmcut)

pca_x_df1<- pca_tpmcut$x %>% as.data.frame() #pca x 자료만 데이터프레임화
pca_x_df2 <- rownames_to_column(pca_x_df, var='index')
rownames(pca_x_df2) <- rownames(pca_x_df1)

pca_x_df <- pca_x_df2 %>% separate(index, into = c("dex","IDnum"), sep="_") %>% select(-IDnum)

#plot 그리기
pca_x_df %>% ggplot(aes(x = PC1, y = PC2, color = dex)) + geom_point()
############


############
# t-SNE
library(Rtsne)
#head(GSE11382_tpmcut)

#GSE11382_tpmcut_CMP <- GSE11382_tpmcut %>% select(contains("CMP"))

#tpm 자료로 tsne 구하기
set.seed(2022)
GSE11382_tpmcut_t <- GSE11382_tpmcut %>% select(-Feature.ID) %>% t()
tsne_tpm_t <- Rtsne(GSE11382_tpmcut_t, perplexity=2, check_duplicates=FALSE)

tsne_tpm_t_df <- as.data.frame(tsne_tpm_t$Y) #데이터프레임화
tsne_tpm_t_df <- data.frame(tsne_tpm_t_df, index)
head(tsne_tpm_t_df)

tsne_tpm_t_df %>% ggplot(aes(x = V1, y = V2, color=index)) + geom_point()
############


############
## Heatmap 그리기
# 발현량 평균 높은 유전자만 뽑아내기
GSE11382_tpmcut1 <- GSE11382_tpm %>% mutate(rowmeans = rowMeans(GSE11382_tpm[,-1])) %>%
	filter(rowmeans > 10) %>% select(-rowmeans)
set.seed(2022)
GSE11382_random <- GSE11382_tpmcut1 %>% sample_n(300)
GSE11382_random_matrix <- GSE11382_random %>% select(-Feature.ID) %>% as.matrix()
rownames(GSE11382_random_matrix) <- GSE11382_random$Feature.ID

heatmap(GSE11382_random_matrix, scale = "row")
library(pheatmap)
pheatmap(GSE11382_random_matrix, scale="row")


# Grouping by ID
library(pheatmap)
GSE11382_random <- GSE11382_tpmcut1 %>% sample_n(100)

#ID132
GSE11382_random_132 <- GSE11382_random %>% select(contains("ID132"))
rownames(GSE11382_random_132) <- GSE11382_random$Feature.ID
GSE11382_random_matrix_132 <- as.matrix(GSE11382_random_132)

heatmap(GSE11382_random_matrix_132, scale = "row")
pheatmap(GSE11382_random_matrix_132, scale="row")

#ID185
GSE11382_random_185 <- GSE11382_random %>% select(contains("ID185"))
rownames(GSE11382_random_185) <- GSE11382_tpmcut1$Feature.ID
GSE11382_random_matrix_185 <- as.matrix(GSE11382_random_185)

heatmap(GSE11382_random_matrix_185, scale = "row")
pheatmap(GSE11382_random_matrix_185, scale="row")

#ID234
GSE11382_random_234 <- GSE11382_random %>% select(contains("ID234"))
rownames(GSE11382_random_234) <- GSE11382_tpmcut1$Feature.ID
GSE11382_random_matrix_234 <- as.matrix(GSE11382_random_234)

heatmap(GSE11382_random_matrix_234, scale = "row")
pheatmap(GSE11382_random_matrix_234, scale="row")

#ID235
GSE11382_random_235 <- GSE11382_random %>% select(contains("ID235"))
rownames(GSE11382_random_235) <- GSE11382_tpmcut1$Feature.ID
GSE11382_random_matrix_235 <- as.matrix(GSE11382_random_235)

heatmap(GSE11382_random_matrix_235, scale = "row")
pheatmap(GSE11382_random_matrix_235, scale="row")

#ID236
GSE11382_random_236 <- GSE11382_random %>% select(contains("ID236"))
rownames(GSE11382_random_236) <- GSE11382_tpmcut1$Feature.ID
GSE11382_random_matrix_236 <- as.matrix(GSE11382_random_236)

heatmap(GSE11382_random_matrix_236, scale = "row")
pheatmap(GSE11382_random_matrix_236, scale="row")
############


####
#PCA
# tpmcut_round <- round(GSE11382_tpmcut_1)
# GSE11382_tpmround <- data.frame(GSE11382_tpmcut$Feature.ID, tpmcut_round)
# colnames(GSE11382_tpmround)[1] <- 'Feature.ID'

# dds_round <- DESeqDataSetFromMatrix(countData = GSE11382_tpmround,
#                                colData = metadata,
#                                design = ~index, tidy = TRUE)
# dds_tpmround <- DESeq(dds_round)
# res_tpmround <- results(dds_tpmround)
# head(res_tpmround)

# vsdata_round <- vst(dds_round, blind = FALSE) 
# plotPCA(vsdata_round, intgroup="index")
