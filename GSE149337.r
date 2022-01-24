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


