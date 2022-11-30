char_data <- read.csv("counts_miRNA.txt", header=TRUE, sep="\t",row.names="miRNA",stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
num_data <- data.frame(data.matrix(char_data),check.names=FALSE) #finds numerical data.
numeric_columns <- sapply(num_data,function(x){mean(as.numeric(is.na(x)))<0.5}) #converts numeric data to numeric type (assuming columns with less than 50% NAs upon converting to numeric are indeed numeric).
cts <- data.frame(num_data[,numeric_columns], char_data[,!numeric_columns],check.names=FALSE) #merges character data and numeric data (note: it rearranges the column order).
coldata <- read.csv("annotations.csv",row.names=1)
coldata <- coldata[,c("Group","Sex","Genotype","Treatment","Pairs")]
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData=cts, colData=coldata, design =~Group)
dds$Genotype <- factor(dds$Genotype, levels = c("NTg","HetWT","HetMT","HomMT"))
dds$Treatment <- factor(dds$Treatment, levels = c("Vehicle","LPS"))
dds$Group <- factor(paste0(dds$Sex, dds$Genotype, dds$Treatment))
dds$Genotype <- relevel(dds$Genotype, ref = "NTg")
dds$Treatment <- relevel(dds$Treatment, ref = "Vehicle")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
design(dds) <- ~Genotype + Treatment + Genotype:Treatment
dds <- DESeq(dds)

#main effect of treatment
resMain <- results(dds, contrast=c("Treatment","LPS","Vehicle"))

#effect of genotype without treatment
resGenoHomMTvsNTg = results(dds, contrast=c("Genotype","HomMT","NTg"))
resGenoHetMTvsNTg = results(dds, contrast=c("Genotype","HetMT","NTg"))
resGenoHetWTvsNTg = results(dds, contrast=c("Genotype","HetWT","NTg"))

#main effect of treatment plus interaction
resHomMTvsNTg <- results(dds, contrast=list( c("Treatment_LPS_vs_Vehicle","GenotypeHomMT.TreatmentLPS")))
resHetMTvsNTg <- results(dds, contrast=list( c("Treatment_LPS_vs_Vehicle","GenotypeHetMT.TreatmentLPS")))
resHetWTvsNTg <- results(dds, contrast=list( c("Treatment_LPS_vs_Vehicle","GenotypeHetWT.TreatmentLPS")))

#interaction for reference Genotype
resIntHomMTvsNTg <- results(dds, name="GenotypeHomMT.TreatmentLPS")
resIntHetMTvsNTg <- results(dds, name="GenotypeHetMT.TreatmentLPS")
resIntHetWTvsNTg <- results(dds, name="GenotypeHetWT.TreatmentLPS")

#interaction for other Genotypes
resIntHomMTvsHetMT <- results(dds, contrast=list("GenotypeHomMT.TreatmentLPS", "GenotypeHetMT.TreatmentLPS"))
resIntHomMTvsHetWT <- results(dds, contrast=list("GenotypeHomMT.TreatmentLPS", "GenotypeHetWT.TreatmentLPS"))
resIntHetMTvsHetWT <- results(dds, contrast=list("GenotypeHetMT.TreatmentLPS", "GenotypeHetWT.TreatmentLPS"))


#Save all data
write.csv(as.data.frame(resMain), file="resMain.csv")

write.csv(as.data.frame(resGenoHomMTvsNTg), file="resGenoHomMTvsNTg.csv")
write.csv(as.data.frame(resGenoHetMTvsNTg), file="resGenoHetMTvsNTg.csv")
write.csv(as.data.frame(resGenoHetWTvsNTg), file="resGenoHetWTvsNTg.csv")

write.csv(as.data.frame(resHomMTvsNTg), file="resHomMTvsNTg.csv")
write.csv(as.data.frame(resHetMTvsNTg), file="resHetMTvsNTg.csv")
write.csv(as.data.frame(resHetWTvsNTg), file="resHetWTvsNTg.csv")

write.csv(as.data.frame(resIntHomMTvsNTg), file="resIntHomMTvsNTg.csv")
write.csv(as.data.frame(resIntHetMTvsNTg), file="resIntHetMTvsNTg.csv")
write.csv(as.data.frame(resIntHetWTvsNTg), file="resIntHetWTvsNTg.csv")

write.csv(as.data.frame(resIntHomMTvsHetMT), file="resIntHomMTvsHetMT.csv")
write.csv(as.data.frame(resIntHomMTvsHetWT), file="resIntHomMTvsHetWT.csv")
write.csv(as.data.frame(resIntHetMTvsHetWT), file="resIntHetMTvsHetWT.csv")


#Save only the data passing the significance thresholds
resOrdered <- resMain[order(resMain$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resMainDys.csv")

resOrdered <- resGenoHomMTvsNTg[order(resGenoHomMTvsNTg$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resGenoHomMTvsNTgDys.csv")

resOrdered <- resGenoHetMTvsNTg[order(resGenoHetMTvsNTg$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resGenoHetMTvsNTgDys.csv")

resOrdered <- resGenoHetWTvsNTg[order(resGenoHetWTvsNTg$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resGenoHetWTvsNTgDys.csv")

resOrdered <- resHomMTvsNTg[order(resHomMTvsNTg$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resHomMTvsNTgDys.csv")

resOrdered <- resHetMTvsNTg[order(resHetMTvsNTg$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resHetMTvsNTgDys.csv")

resOrdered <- resHetWTvsNTg[order(resHetWTvsNTg$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resHetWTvsNTgDys.csv")

resOrdered <- resIntHomMTvsNTg[order(resIntHomMTvsNTg$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resIntHomMTvsNTgDys.csv")

resOrdered <- resIntHetMTvsNTg[order(resIntHetMTvsNTg$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resIntHetMTvsNTgDys.csv")

resOrdered <- resIntHetWTvsNTg[order(resIntHetWTvsNTg$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resIntHetWTvsNTgDys.csv")

resOrdered <- resIntHomMTvsHetMT[order(resIntHomMTvsHetMT$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resIntHomMTvsHetMTDys.csv")

resOrdered <- resIntHomMTvsHetWT[order(resIntHomMTvsHetWT$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resIntHomMTvsHetWTDys.csv")

resOrdered <- resIntHetMTvsHetWT[order(resIntHetMTvsHetWT$log2FoldChange),]
res05 <- subset(resOrdered, pvalue <0.05)
resDys <- subset(res05, log2FoldChange >0.58 | log2FoldChange <(-0.58))
write.csv(as.data.frame(resDys), file="resIntHetMTvsHetWTDys.csv")


#Save normalised counts
normalised_counts <- counts(dds, normalized=TRUE)
write.csv(as.data.frame(normalised_counts), file="normCounts.csv")