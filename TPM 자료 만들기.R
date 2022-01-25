############
## TPM 자료 만들기
load('/spstorage/USERS/sung/projects/DEGanalysis/Genelength.RData')
head(genelength)

head(GSE11382_pcg)

genelength_1 <- rownames_to_column(genelength, var='Feature.ID')
head(genelength_1)

GSE11382_length <- GSE11382_pcg %>% left_join(genelength_1, by = "Feature.ID")
head(GSE11382_length)

GSE11382_genes <- GSE11382_length[,c(1,length(GSE11382_length))]
head(GSE11382_genes)
GSE11382_counts <- GSE11382_length[,-c(1,length(GSE11382_length))]
head(GSE11382_counts)

tpm <- function(counts, lengths){
	rate <- counts / lengths
	rate / sum(rate) * 1e6
}

GSE11382_tpm_counts <- tpm(GSE11382_counts, GSE11382_genes$length)
head(GSE11382_tpm_counts)

GSE11382_tpm <- data.frame(GSE11382_genes$Feature.ID, GSE11382_tpm_counts)
colnames(GSE11382_tpm)[1] <- 'Feature.ID'
head(GSE11382_tpm)
############
