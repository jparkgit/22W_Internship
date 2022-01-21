gse149377 <- read.csv(file="GSE149377_iPSCderived_myeloid_cells_counts.gct.gz")
genes.tbl <- getBM(attributes=c(“ensembl_gene_id_version”, “hgnc_symbol”), filters=“ensembl_gene_id_version”, values=ensembl_ids, mart=ens)

# 1. Install 'biomaRt' package via bioconductor
source(“http://Bioconductor.org/biocLite.R”)
biocLite(c(“biomaRt”, “xml2”))


# 2. Read a count matrix
annot.tbl <- read.csv(“”, header=TRUE)
ensembl_ids <- unlist(lapply(annot.tbl$gene_id, as.character), use.names=FALSE)


# 3. Get a list of 'official gene symbol'
library(biomaRt)
ens <- useEnsembl(biomart=“ensembl”, dataset=“hsapiens_gene_ensembl”)
genes.tbl <- getBM(attributes=c(“ensembl_gene_id_version”, “hgnc_symbol”),
                   filters=“ensembl_gene_id_version”,
                   values=ensembl_ids, mart=ens)


# 4. Convert the ensembl to symbol
merge.tbl <- merge(x=annot.tbl, y=genes.tbl,
                   by.x=“gene_id”, by.y=“ensembl_gene_id_version”,
                   all.x=TRUE, all.y=FALSE)
merge.tbl <- merge.tbl[,c(1,14,seq(2,13))]


# 5. Save the result
write.csv(merge.tbl, file=“output_biomart/merged_gene_count_matrix.csv”, row.names=FALSE)


gse149377 <- read.delim(file="GSE149377_iPSCderived_myeloid_cells_counts.gct.gz", skip=2)

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


#covert gene ids with biomart

# Fit your model on your expressionset and design matrix
fit  <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2)

# get biomart annotation
library(biomaRt)
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
attributes <- c("affy_mogene_1_0_st_v1", "chromosome_name",
                "start_position","end_position", "ensembl_gene_id", 
                "external_gene_id", "description")
attributes <- c("affy_mogene_1_0_st_v1", "chromosome_name",
                "start_position","end_position", "ensembl_gene_id", 
                "external_gene_id", "description")

genes <- getBM(attributes=attributes, filters="affy_mogene_1_0_st_v1", 
               values=topTable$ID, mart=mart, uniqueRows=T)


library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

attributes <- c("ensembl_gene_id", "entrezgene_id", "entrezgene_accession")

genes <- getBM(
	attributes = attributes, 
	filters = "ensembl_gene_id", 
	values= c('ENSG00000228288'), 
	mart = mart, 
	uniqueRows = F)

ensembl_gene_id
entrezgene_id
entrezgene_accession


affyids=c("202763_at","209310_s_at","207500_at")


getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene'), 
      filters = 'affy_hg_u133_plus_2', 
      values = affyids, 
      mart = ensembl)