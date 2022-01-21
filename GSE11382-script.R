## countdata 불러오기
data <- read.table("GSE113182_raw_count.txt", header = T, sep='\t')
head(data)


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
## PCG filtering
# 전체 유전자의 info가 담긴 gmt파일을 불러와 그중 pcg만 필터링
ann = rtracklayer::import("/spstorage/DB/NGS_reference/GRCh38_v34/gencode.v34.primary_assembly.annotation.gtf") %>% as.data.frame()
pcg = ann %>% filter(gene_type == "protein_coding") %>% distinct(gene_id)
pcg$gene_id <- gsub("\\.[0-9]+","",pcg$gene_id)

# GSE데이터에서 pcg 필터링
pcglist <- ifelse(GSE11382_counts$Feature %in% pcg$gene_id, "True", "False")
GSE11382_pcg <- GSE11382_counts[which(pcglist == 'True'),]
head(GSE11382_pcg)
############


############
## 유전자 발현 합 평균 5 이상 Filtering
#HSC: 유전자 발현 합이 평균 5이상 Filter
GSE11382_pcg %>% select(contains("HSC")) %>% head
GSE11382_HSC <- GSE11382_pcg %>% select(Feature.ID, contains("HSC")) %>%
				mutate(rowSUM = (HSC_ID236 + HSC_ID235 + HSC_ID185 + HSC_ID132 + HSC_ID234)/5) %>%
				filter(rowSUM >= 5) %>% select(Feature.ID:HSC_ID234)
#CMP: 유전자 발현 합이 평균 5이상 Filter
GSE11382_pcg %>% select(contains("CMP")) %>% head
GSE11382_CMP <- GSE11382_pcg %>% select(Feature.ID, contains("CMP")) %>%
				mutate(rowSUM = (CMP_ID236 + CMP_ID235 + CMP_ID185 + CMP_ID132 + CMP_ID234)/5) %>%
	 			filter(rowSUM >= 5) %>% select(Feature.ID:CMP_ID234)
#GMP: 유전자 발현 합이 평균 5이상 Filter
GSE11382_pcg %>% select(contains("GMP")) %>% head
GSE11382_GMP <- GSE11382_pcg %>% select(Feature.ID, contains("GMP")) %>%
				mutate(rowSUM = (GMP_ID236 + GMP_ID235 + GMP_ID185 + GMP_ID132 + GMP_ID234)/5) %>%
				filter(rowSUM >= 5) %>% select(Feature.ID:GMP_ID234)
#MEP: 유전자 발현 합이 평균 5이상 Filter
GSE11382_pcg %>% select(contains("MEP")) %>% head
GSE11382_MEP <- GSE11382_pcg %>% select(Feature.ID, contains("MEP")) %>%
				mutate(rowSUM = (MEP_ID236 + MEP_ID235 + MEP_ID185 + MEP_ID132 + MEP_ID234)/5) %>%
				filter(rowSUM >= 5) %>% select(Feature.ID:MEP_ID234)
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


# meta_data 1 & DEG (HSC -> CMP)
ID_1 <- c("HSC_ID236","HSC_ID235","HSC_ID185","HSC_ID132", "HSC_ID234",
	"CMP_ID236","CMP_ID235", "CMP_ID185","CMP_ID132","CMP_ID234")
index_1 <- rep(c("hematopoietic_stem_cells","common_myeloid_progenitors"),each=5)

metadata_1 <- data.frame(index_1)
rownames(metadata_1) <- ID_1 

dds1 <- DESeqDataSetFromMatrix(countData = HSC_CMP ,
                               colData = metadata_1,
                               design = ~index_1, tidy = TRUE)

dds_HSC_CMP <- DESeq(dds1)
res_HSC_CMP <- results(dds_HSC_CMP)

res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 2, padj < 0.05) %>% dim
res_HSC_CMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 2, padj < 0.05) %>% head


# meta data 2 & DEG (CMP -> GMP)
ID_2 <- c("CMP_ID236", "CMP_ID235", "CMP_ID185", "CMP_ID132", "CMP_ID234", "GMP_ID236",
	 "GMP_ID235", "GMP_ID185", "GMP_ID132", "GMP_ID234")
index_2 <- rep(c("common_myeloid_progenitors","granulocyte-macrophage_progenitors"),each=5)

metadata_2 <- data.frame(index_2)
rownames(metadata_2) <- ID_2 

dds2 <- DESeqDataSetFromMatrix(countData = CMP_GMP ,
                               colData = metadata_2,
                               design = ~index_2, tidy = TRUE)

dds_CMP_GMP <- DESeq(dds2)
res_CMP_GMP <- results(dds_CMP_GMP)

res_CMP_GMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% dim   #334개 
res_CMP_GMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% head
res_CMP_GMP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1.5, padj < 0.05) %>% head #213개 


# meta data 3 & DEG (CMP -> MEP)
ID_3 <- c("CMP_ID236", "CMP_ID235", "CMP_ID185", "CMP_ID132", "CMP_ID234", "MEP_ID236",
	"MEP_ID235", "MEP_ID185", "MEP_ID132", "MEP_ID234")
index_3 <- rep(c("common_myeloid_progenitors","megakaryocyte-erythrocyte_progenitors"),each=5)

metadata_3 <- data.frame(index_3)
rownames(metadata_3) <- ID_3

dds3 <- DESeqDataSetFromMatrix(countData = CMP_MEP ,
                               colData = metadata_3,
                               design = ~index_3, tidy = TRUE)
dds_CMP_MEP <- DESeq(dds3)
res_CMP_MEP <- results(dds_CMP_MEP)

res_CMP_MEP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% dim #201개 
res_CMP_MEP %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% head
res_CMP_MEP %>% as.data.frame %>% filter(abs(log2FoldChange) > 2, padj < 0.05) %>% dim   #64개 
############


############
#PCA plot 
#HSC_CMP
vsdata_1 <- vst(dds_HSC_CMP, blind = FALSE) 
plotPCA(vsdata_1, intgroup="index_1")