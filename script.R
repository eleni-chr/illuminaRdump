##Section 1 - Run this section on its own
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)

rm(list = ls())

mmmart <- useMart(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")

#!!!!!!!!!!!!!!!!CHANGE FILENAME ACCORDINGLY!!!!!!!!!!!!!!!!
data <- read.csv("resMainDys_LFC.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE,check.names=FALSE)
mapping <- getBM(attributes = c('entrezgene_id', 'mirbase_id'), filters = 'mirbase_id', values = data[,1], mart = mmmart)
#mapping <- mapping[-c(50),] 

#!!!!!!!!!!!!!!!!DO MANUAL CHECK HERE!!!!!!!!!!!!!!!!
View(data)
View(mapping)

##Section 2 - Run this section on its own
genes <- as.character(mapping$entrezgene_id)
ego <- enrichGO(gene = genes, OrgDb = org.Mm.eg.db, ont = "ALL", minGSSize = 0, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
egoRes <- as.data.frame(ego)
ordered <-egoRes[order(egoRes$p.adjust),]

#!!!!!!!!!!!!!!!!CHANGE FILENAME ACCORDINGLY!!!!!!!!!!!!!!!!
write.csv(as.data.frame(ordered), file="resMainDys_LFC_ego.csv")