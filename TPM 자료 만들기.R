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


############
# t-SNE
index <- rep(c("hematopoietic_stem_cells","common_myeloid_progenitors",
	"granulocyte-macrophage_progenitors","megakaryocyte-erythrocyte_progenitors"), each=5)
metadata <- data.frame(index)
rownames(metadata) <- colnames(GSE11382_tpm)[-1]

# dds <- DESeqDataSetFromMatrix(countData = GSE11382_tpm,
#                                colData = metadata,
#                                design = ~index, tidy = TRUE)
# dds_tpm <- DESeq(dds)
# res_tpm <- results(dds_tpm)
# head(res_tpm)

library(Rtsne)
GSE11382_tpm_sample <- GSE11382_tpm %>% sample_n(1000)
tsne_tpm <- Rtsne(GSE11382_tpm_sample, perplexity=30, check_duplicates=FALSE)
plot(tsne_tpm$Y, col = "black", bg = metadata$index, pch = 21, cex = 1)
