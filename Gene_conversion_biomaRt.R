library(biomaRt)
hu.ensembl = useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl') 
listFilters(hu.ensembl) 

#Reference gene list: containing entrezgene id and ensembl gene id
G_list <- getBM(attributes= c("entrezgene_id","ensembl_gene_id"),values=df$Geneid, mart = hu.ensembl, useCache = FALSE)

#New data frame with a column with our entrezgene id and to store a corresponding ensembl gene id
df <- data.frame(gse149377$Name)
names(df)[1] <- 'Geneid'

#Create an empty column to store ensembl gene id
df$ensembl_gene_id = ""

#Input ensembl gene id to df column if entrez gene id from df matches with entrez gene id from G_list 
df["ensembl_gene_id"] = lapply("ensembl_gene_id", function(x) G_list[[x]][match(df$Geneid, G_list$entrezgene_id)])