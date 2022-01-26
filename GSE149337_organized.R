

gse149377 <- read.delim(file="GSE149377_iPSCderived_myeloid_cells_counts.gct.gz", skip=2)






###covert gene ids with biomart

library(biomaRt)
hu.ensembl = useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl') 

#Reference gene list: containing entrezgene id and ensembl gene id
G_list <- getBM(attributes= c("entrezgene_id","ensembl_gene_id"),values=gse149377$Name, mart = hu.ensembl, useCache = FALSE)

#create new column
gse149377$ensembl_gene_id = ""

#input ensemble gene id
gse149377["ensembl_gene_id"] = lapply("ensembl_gene_id", function(x) G_list[[x]][match(gse149377$Name, G_list$entrezgene_id)])










### Pre-filtering [PCG filtering]
### 전체 유전자의 info가 담긴 gmt파일을 불러와 그중 pcg만 필터링

	library(dplyr)

	ann = rtracklayer::import("/spstorage/DB/NGS_reference/GRCh38_v34/gencode.v34.primary_assembly.annotation.gtf") %>% as.data.frame()
	pcg = ann %>% filter(gene_type == "protein_coding") %>% distinct(gene_id)
	pcg$gene_id <- gsub("\\.[0-9]+","",pcg$gene_id)

	pcglist <- ifelse(gse149377$ensembl_gene_id %in% pcg$gene_id, "True", "False")
	gse149377_pcg <- gse149377[which(pcglist == 'True'),]










### TPM

counts <- select(gse149377_pcg, -ensembl_gene_id)
counts <- gse149377_pcg[,-c(1,2)]

tpm <- function(counts, length){
	rpk <- counts/(length/1000)
	coef <- sum(rpk)/1e6
	rpk/coef
}

tpms <- apply(counts, 2, function(x) tpm(x, gse149377_w_length$length))

head(tpms)
colSums(tpms) #all equal to 1e+06

#tpms does not include ensembl gene id
#need to join


#tpm dataframe 만들기
gse149377_tpm <- data.frame(gse149377_pcg$ensembl_gene_id, tpms)
colnames(gse149377_tpm)[1] <- 'ensembl_gene_id'
head(gse149377_tpm)


#유전자 발현 평균 2 이상 filtering
gse149377_tpmcut <-  gse149377_tpm %>% mutate(rowmeans = rowMeans(gse149377_tpm[,-1])) %>% filter(rowmeans > 2) %>% select(-rowmeans)
head(gse149377_tpmcut)


dim(gse149377_tpm)  # before cut : 17172    51
dim(gse149377_tpmcut) # after cut: 10771    51




# PCA-Plot Plotting 	#############################################################################

### meta data ####

names <- colnames(gse149377_tpmcut_1)
meta_total <- data.frame(names)
meta_total$index1 = ""
meta_total$index1 <- ifelse(grepl("PBMC", meta_total$names), "PBMC", "IPS")
meta_total$index2 = ""
meta_total$index2 <- ifelse(grepl("m0", meta_total$names), "m0","") 
meta_total$index3 <- ifelse(grepl("m1", meta_total$names), "m1","")
meta_total$index4 <- ifelse(grepl("m2", meta_total$names), "m2","")
meta_total$index5 <- ifelse(grepl("mono", meta_total$names), "mono","")
meta_total$index6 <- ifelse(grepl("micro", meta_total$names), "micro","")
my_cols <- c("index2", "index3", "index4", "index5", "index6")
meta_total$index0 <- do.call(paste, c(meta_total[my_cols], sep = ""))
meta_total <- meta_total[ , ! colnames(meta_total) %in% my_cols]
rownames(meta_total) <- meta_total$names
colnames(meta_total)[1] <- 'index'




library(tidyverse)
set.seed(2022)

ips_pcg <- gse149377_tpm %>% select(ensembl_gene_id, contains("IPS"))
pbmc_pcg <- gse149377_tpm %>% select(ensembl_gene_id, contains("PBMC"))


gse149377_tpmcut_1 <- gse149377_tpmcut[,-1] #Feature.ID 제거
dim(gse149377_tpmcut_1) # 10771 50
gse149377_tpmcut_2 = gse149377_tpmcut_1 %>% t() #행, 열 바꾸기
dim(gse149377_tpmcut_2) # 50 17172

pca_tpmcut <- prcomp(gse149377_tpmcut_2) #pca 구하기
summary(pca_tpmcut)


pca_x_df1 <- pca_tpmcut$x %>% as.data.frame() #pca x 자료만 데이터프레임화
pca_x_df2 <- rownames_to_column(pca_x_df1, var='index') #df1의 rownames을 df2에게 index라는 column으로 주기
rownames(pca_x_df2) <- rownames(pca_x_df1)  #df1의 rownames을 df2의 rownames로 주기

pca_x_df <- meta_total %>% inner_join(pca_x_df2, by = "index") #df2와 meta total을 index를 겹쳐서 join시키기

pca_x_df %>% ggplot(aes(x = PC1, y = PC2, color = index0)) + geom_point()
pca_x_df %>% ggplot(aes(x = PC1, y = PC2, color = index1)) + geom_point()



















### IPS vs PMBC #############################################################################
### 여기서부터 PBMC, IPS 나눠서 작업: create separate data for each tye


# m0
temp <- pbmc_pcg %>% select(ensembl_gene_id, contains("m0"))
temp$Mean <- apply(temp[,-1], 1, mean)
pbmc_m0 <- temp %>% filter(Mean >= 5) %>% select(-Mean)
temp <- NULL
dim(pbmc_m0) # 8936 6

temp <- ips_pcg %>% select(ensembl_gene_id, contains("m0"))
temp$Mean <- apply(temp[,-1], 1, mean)
ips_m0 <- temp %>% filter(Mean >= 5) %>% select(-Mean)
temp <- NULL
dim(ips_m0) # 9422 6


# m1
temp <- pbmc_pcg %>% select(ensembl_gene_id, contains("m1"))
temp$Mean <- apply(temp[,-1], 1, mean)
pbmc_m1 <- temp %>% filter(Mean >= 5) %>% select(-Mean)
temp <- NULL
dim(pbmc_m1) #8747 6

temp <- ips_pcg %>% select(ensembl_gene_id, contains("m1"))
temp$Mean <- apply(temp[,-1], 1, mean)
ips_m1 <- temp %>% filter(Mean >= 5) %>% select(-Mean)
temp <- NULL
dim(ips_m1) # 8958 6



# m2
temp <- pbmc_pcg %>% select(ensembl_gene_id, contains("m2"))
temp$Mean <- apply(temp[,-1], 1, mean)
pbmc_m2 <- temp %>% filter(Mean >= 5) %>% select(-Mean)
temp <- NULL
dim(pbmc_m2) #8744     6

temp <- ips_pcg %>% select(ensembl_gene_id, contains("m2"))
temp$Mean <- apply(temp[,-1], 1, mean)
ips_m2 <- temp %>% filter(Mean >= 5) %>% select(-Mean)
temp <- NULL
dim(ips_m2) # 9263 6


# monocyte
temp <- pbmc_pcg %>% select(ensembl_gene_id, contains("mono"))
temp$Mean <- apply(temp[,-1], 1, mean)
pbmc_mono <- temp %>% filter(Mean >= 5) %>% select(-Mean)
temp <- NULL
dim(pbmc_mono) # 8653     9

temp <- ips_pcg %>% select(ensembl_gene_id, contains("mono"))
temp$Mean <- apply(temp[,-1], 1, mean)
ips_mono <- temp %>% filter(Mean >= 5) %>% select(-Mean)
temp <- NULL
dim(ips_mono) #9439 8


# PBMC는 microglia 없음



#	IPS vs PBMC 	#############################################################################
mono <- pbmc_mono %>% inner_join(ips_mono, by = "ensembl_gene_id")
m0 <- pbmc_m0 %>% inner_join(ips_m0, by = "ensembl_gene_id")
m1 <- pbmc_m1 %>% inner_join(ips_m1, by = "ensembl_gene_id")
m2 <- pbmc_m2 %>% inner_join(ips_m2, by = "ensembl_gene_id")


### create metadata
	#1. mono
	names <- NULL
	names <- colnames(mono)
	names <- names[-1]
	meta_mono <- data.frame(names)
	meta_mono$typemono = ""
	meta_mono$typemono <- ifelse(grepl("PBMC", meta_mono$names), "PBMC", "IPS")
	rownames(meta_mono) <- meta_mono$names
	meta_mono <- select(meta_mono, -names)
	meta_mono

	#2. m0
	names <- NULL
	names <- colnames(m0)
	names <- names[-1]
	meta_m0 <- data.frame(names)
	meta_m0$typem0 = ""
	meta_m0$typem0 <- ifelse(grepl("PBMC", meta_m0$names), "PBMC", "IPS")
	rownames(meta_m0) <- meta_m0$names
	meta_m0 <- select(meta_m0, -names)
	meta_m0

	#3. m1
	names <- NULL
	names <- colnames(m1)
	names <- names[-1]
	meta_m1 <- data.frame(names)
	meta_m1$typem1 = ""
	meta_m1$typem1 <- ifelse(grepl("PBMC", meta_m1$names), "PBMC", "IPS")
	rownames(meta_m1) <- meta_m1$names
	meta_m1 <- select(meta_m1, -names)
	meta_m1

	#4. m2
	names <- NULL
	names <- colnames(m2)
	names <- names[-1]
	meta_m2 <- data.frame(names)
	meta_m2$typem2 = ""
	meta_m2$typem2 <- ifelse(grepl("PBMC", meta_m2$names), "PBMC", "IPS")
	rownames(meta_m2) <- meta_m2$names
	meta_m2 <- select(meta_m2, -names)
	meta_m2

### create DESeq Data set 
	library(DESeq2)
	typemono <- meta_mono$typemono
	typem0 <- meta_m0$typem0
	typem1 <- meta_m1$typem1
	typem2 <- meta_m2$typem2









dds1 <- DESeqDataSetFromMatrix(countData = select(mono, -1), 
colData = meta_mono, 
design = ~typemono)
dds_mono <- DESeq(dds1)
res_mono <- results(dds_mono)
head(res_mono)

dds2 <- DESeqDataSetFromMatrix(countData = select(m0, -1), 
colData = meta_m0, 
design = ~typem0)
dds_m0 <- DESeq(dds2)
res_m0 <- results(dds_m0)
head(res_m0)


dds3 <- DESeqDataSetFromMatrix(countData = select(m1, -1), 
colData = meta_m1, 
design = ~typem1)
dds_m1 <- DESeq(dds3)
res_m1 <- results(dds_m1)
head(res_m1)

dds4 <- DESeqDataSetFromMatrix(countData = select(m2, -1), 
colData = meta_m2, 
design = ~typem2)
dds_m2 <- DESeq(dds4)
res_m2 <- results(dds_m2)
head(res_m2)



res1 <- res_mono %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) # 11774 -> 3191
res2 <- res_m0 %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) # 12367 -> 2512
res3 <- res_m1 %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) # 12354 -> 2248
res4 <- res_m2 %>% as.data.frame %>% filter(abs(log2FoldChange) > 1, padj < 0.05) # 12354 -> 2248



##log2foldchange cuf off
res_pbmc_m0_m1 %>% as.data.frame %>% filter(abs(log2FoldChange) > 5, padj < 0.05) %>% dim
#mono_m0 
	# n=3: 638
	# n=4: 322
	# n=5: 179
#m0 m1
	# n=4: 239
	# n=5: 84
#m0 m2
	# n=4: 260
	# n=5: 126




##PBMC 		#############################################################################
### create DESeq Data set 
	library(DESeq2)type1 <- meta_pbmc_mono_m0$type1
	type2 <- meta_pbmc_mono_m0$type2
	type3 <- meta_pbmc_mono_m0$type3

##IPS 		#############################################################################







