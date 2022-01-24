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

index <- c("common_myeloid_progenitors","common_myeloid_progenitors","common_myeloid_progenitors",
	"granulocyte-macrophage_progenitors","ggranulocyte-macrophage_progenitors",
	"granulocyte-macrophage_progenitors","hematopoietic_stem_cells","hematopoietic_stem_cells",
	"hematopoietic_stem_cells","megakaryocyte-erythrocyte_progenitors",
	"megakaryocyte-erythrocyte_progenitors","megakaryocyte-erythrocyte_progenitors",
	"common_myeloid_progenitors","common_myeloid_progenitors","granulocyte-macrophage_progenitors",
	"granulocyte-macrophage_progenitors","hematopoietic_stem_cells","hematopoietic_stem_cells",
	"megakaryocyte-erythrocyte_progenitors","megakaryocyte-erythrocyte_progenitors")

metadata <- data.frame(index)
rownames(metadata) <- ID
############


############
## Pre-filtering [PCG filtering]
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

#RES_HSC_CMP <- results(dds_HSC_CMP, contrast=c("index_1", "hematopoietic_stem_cells", "common_myeloid_progenitors"))


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


############
## PCA plot 그리기 
#HSC -> CMP
vsdata_1 <- vst(dds_HSC_CMP, blind = FALSE) 
plotPCA(vsdata_1, intgroup="index_1")

#CMP -> GMP
vsdata_2 <- vst(dds_CMP_GMP, blind = FALSE) 
plotPCA(vsdata_2, intgroup="index_2")

#CMP -> MEP
vsdata_3 <- vst(dds_CMP_MEP, blind = FALSE) 
plotPCA(vsdata_3, intgroup="index_3")
############


############
## t-SNE
library(Rtsne)

#whole data로 DESeq
whole_data <- GSE11382_HSC %>% inner_join(GSE11382_CMP, by = "Feature.ID") %>%
	inner_join(GSE11382_GMP, by = "Feature.ID") %>% inner_join(GSE11382_MEP, by = "Feature.ID")

ID <- colnames(whole_data)[-1]
index <- rep(c("hematopoietic_stem_cells","common_myeloid_progenitors",
	"granulocyte-macrophage_progenitors","megakaryocyte-erythrocyte_progenitors"), each=5)
metadata <- data.frame(index)
rownames(metadata) <- ID

dds <- DESeqDataSetFromMatrix(countData = whole_data ,
                               colData = metadata,
                               design = ~index, tidy = TRUE)
dds_whole <- DESeq(dds)
res_whole <- results(dds_whole)
head(res_whole)

res_whole %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% dim   
RES_WHOLE <- res_whole %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) 

rownames_all <- RES_WHOLE %>% rownames
WHOLE_DEG <- whole_data %>% filter(Feature.ID %in% rownames_all)

tsne_WHOLE <- Rtsne(WHOLE_DEG, perplexity=30, check_duplicates=FALSE)
plot(tsne_WHOLE$Y, col = "black", bg = metadata$index, pch = 21, cex = 1)
tsne_WHOLE <- Rtsne(WHOLE_DEG, perplexity=10, check_duplicates=FALSE)
plot(tsne_WHOLE$Y, col = "black", bg = metadata$index, pch = 21, cex = 1)





### TPM
tpm <- function(counts, lengths){
	rpk <- counts/(lengths/1000)
	coef <- sum(rpk)
	rpk/coef
}
tpms <- apply(counts, 2, function(x) tpm(x, lenghts))



##Counting TPM
summary(whole_data[,2:length(whole_data)])
geneLengths <- as.vector(subset(whole_data, select = c(-Feature.ID)))

rpk <- apply(subset(whole_data, select = c(-Feature.ID)), 2,
	function(x) x/(geneLengths/1000))
tpm <- apply(rpk, 2, function(x) x/sum(as.numeric(x)) * 10^6)

# HSC_CMP_sample <- HSC_CMP %>% sample_n(100)
# tsne_HSC_CMP <- Rtsne(HSC_CMP_sample, perplexity=30, check_duplicates=FALSE)


# RES_HSC_CMP <- res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1.5, padj < 0.05)
# rownames_1 <- RES_HSC_CMP %>% rownames
# HSC_CMP_DEG <- HSC_CMP %>% filter(Feature.ID %in% rownames_1)

# tsne_HSC_CMP <- Rtsne(HSC_CMP_DEG, perplexity=30, check_duplicates=FALSE)

# plot(tsne_HSC_CMP$Y, col = "black", bg = metadata_1$index_1, pch = 21, cex = 1)



# ###
# whole_data <- GSE11382_HSC %>% inner_join(GSE11382_CMP, by = "Feature.ID") %>%
# 	inner_join(GSE11382_GMP, by = "Feature.ID") %>% inner_join(GSE11382_MEP, by = "Feature.ID")
# whole_data_sp <- whole_data %>% sample_n(1000)
# tsne_whole <- Rtsne(whole_data_sp, perplexity=30, check_duplicates=FALSE)
# plot(tsne_whole$Y, col = "black", bg = metadata$index, pch = 21, cex = 1)