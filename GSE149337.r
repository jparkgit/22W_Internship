

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





## Pre-filtering [PCG filtering]
# 전체 유전자의 info가 담긴 gmt파일을 불러와 그중 pcg만 필터링

library(dplyr)

ann = rtracklayer::import("/spstorage/DB/NGS_reference/GRCh38_v34/gencode.v34.primary_assembly.annotation.gtf") %>% as.data.frame()
pcg = ann %>% filter(gene_type == "protein_coding") %>% distinct(gene_id)
pcg$gene_id <- gsub("\\.[0-9]+","",pcg$gene_id)

# GSE데이터에서 pcg 필터링

pcglist <- ifelse(gse149377$ensembl_gene_id %in% pcg$gene_id, "True", "False")
gse149377_pcg <- gse149377[which(pcglist == 'True'),]


## 여기서부터 PBMC, IPS 나눠서 작업
ips_pcg <- gse149377_pcg %>% select(ensembl_gene_id, contains("IPS"))
pbmc_pcg <- gse149377_pcg %>% select(ensembl_gene_id, contains("PBMC"))

# m0
temp <- pbmc_pcg %>% select(ensembl_gene_id, contains("m0"))
temp$Mean <- apply(temp[,-1], 1, mean)
pbmc_m0 <- temp %>% filter(Mean >= 5) %>% select(-Mean)
temp <- NULL
dim(pbmc_m0) # 12657 6

# m1
temp <- pbmc_pcg %>% select(ensembl_gene_id, contains("m1"))
temp$Mean <- apply(temp[,-1], 1, mean)
pbmc_m1 <- temp %>% filter(Mean >= 5) %>% select(-Mean)
temp <- NULL
dim(pbmc_m1) #12906 6

# m2
temp <- pbmc_pcg %>% select(ensembl_gene_id, contains("m2"))
temp$Mean <- apply(temp[,-1], 1, mean)
pbmc_m2 <- temp %>% filter(Mean >= 5) %>% select(-Mean)
temp <- NULL
dim(pbmc_m2) #12958     6

# monocyte
temp <- pbmc_pcg %>% select(ensembl_gene_id, contains("mono"))
temp$Mean <- apply(temp[,-1], 1, mean)
pbmc_mono <- temp %>% filter(Mean >= 5) %>% select(-Mean)
temp <- NULL
dim(pbmc_mono) # 11996     9

# PBMC는 microglia 없음


## 두 데이터씩 합치기
pbmc_mono_m0 <- pbmc_mono %>% inner_join(pbmc_m0, by = "ensembl_gene_id")
pbmc_m0_m1 <- pbmc_m0 %>% inner_join(pbmc_m1, by = "ensembl_gene_id")
pbmc_m0_m2 <- pbmc_m0 %>% inner_join(pbmc_m2, by = "ensembl_gene_id")




# create metadata
#mono m0
names <- colnames(pbmc_mono_m0)
names <- names[-1]
meta_pbmc_mono_m0 <- data.frame(names)
meta_pbmc_mono_m0$type = ""
meta_pbmc_mono_m0$type <- ifelse(grepl("mono", meta_pbmc_mono_m0$names), "mono", "m0")

#m0 m1
names <- colnames(pbmc_m0_m1)
names <- names[-1]
meta_pbmc_m0_m1 <- data.frame(names)
meta_pbmc_m0_m1$type = ""
meta_pbmc_m0_m1$type <- ifelse(grepl("m0", meta_pbmc_m0_m1$names), "m0", "m1")

#m0 m2
names <- colnames(pbmc_m0_m2)
names <- names[-1]
meta_pbmc_m0_m2 <- data.frame(names)
meta_pbmc_m0_m2$type = ""
meta_pbmc_m0_m2$type <- ifelse(grepl("m0", meta_pbmc_m0_m2$names), "m0", "m2")




cp /spstorage/INTERNSHIP/jpark/prefiltering_done.RData /spstorage/INTERNSHIP/hyeseon














### [Unsolved] Codes for autonaming types of cells 
# first row, since 3rd column
gse149377_names <- gse149377[c(1:1), ,c(3:52)]]
#1 row, 50 column

# columns to list
nameslist <- list()

assign("nameslist", NULL, envir = .GlobalEnv)

df_draft <- data.frame(nameslist)

patterns <- c("m0", "m1", "m2", "monocytes", "microglia")

#i for rows (going through names)
#j for patterns (checking through which matching cell type)


for(i in 1:nrow(df_draft)) {
	for(j in 1:length(patterns)) {
		if (grepl(df_draft[i, 1], patterns[j]) == TRUE) {
			df_draft$DiffType <- patterns[j]
		}
		else {NULL}
	}
}

for(i in 1:nrow(df_draft)) {
	if (grepl(df_draft[i, 1], patterns[0]) == TRUE) {
		df_draft$DiffType <- patterns[j]
	}
	else if (grepl(df_draft[i, 1], patterns[1]) == TRUE) {
		df_draft$DiffType <- patterns[j]
	}
	else {NULL}

}


for(i in 1:ncol(gse149377_names)) {
	for(j in 1:length(patterns)) {
		if (str_contains(gse149377_names[, i], patterns[j]) == TRUE) {

		}
	}
	
}

str_contains(, )


nameslist <- colnames(gse149377)


