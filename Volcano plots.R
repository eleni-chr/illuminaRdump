#Males
char_data <- read.csv("counts_miRNA_M.txt", header=TRUE, sep="\t",row.names="miRNA",stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
num_data <- data.frame(data.matrix(char_data),check.names=FALSE) #finds numerical data.
numeric_columns <- sapply(num_data,function(x){mean(as.numeric(is.na(x)))<0.5}) #converts numeric data to numeric type (assuming columns with less than 50% NAs upon converting to numeric are indeed numeric).
cts <- data.frame(num_data[,numeric_columns], char_data[,!numeric_columns],check.names=FALSE) #merges character data and numeric data (note: it rearranges the column order).
coldata <- read.csv("annotationsM.csv",row.names=1)
coldata <- coldata[,c("Group","Sex","Genotype","Treatment")]
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


#Volcano plots
library(EnhancedVolcano)
png("resMainMale.png", width = 866, height = 553)
EnhancedVolcano(resMain,
                lab = rownames(resMain),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect of Treatment',
                subtitle = 'LPS vs Vehicle in NTg Males',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resMain[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resGenoHomMTvsNTgMale.png", width = 866, height = 553)
EnhancedVolcano(resGenoHomMTvsNTg,
                lab = rownames(resGenoHomMTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect of Genotype',
                subtitle = 'HomMT vs NTg Males (no LPS)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resGenoHomMTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resGenoHetMTvsNTgMale.png", width = 866, height = 553)
EnhancedVolcano(resGenoHetMTvsNTg,
                lab = rownames(resGenoHetMTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect of Genotype',
                subtitle = 'HetMT vs NTg Males (no LPS)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resGenoHetMTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resGenoHetWTvsNTgMale.png", width = 866, height = 553)
EnhancedVolcano(resGenoHetMTvsNTg,
                lab = rownames(resGenoHetWTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect of Genotype',
                subtitle = 'HetWT vs NTg Males (no LPS)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resGenoHetWTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resHomMTvsNTgMale.png", width = 866, height = 553)
EnhancedVolcano(resHomMTvsNTg,
                lab = rownames(resHomMTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect plus interaction',
                subtitle = 'HomMT vs NTg Males (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resHomMTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resHetMTvsNTgMale.png", width = 866, height = 553)
EnhancedVolcano(resHetMTvsNTg,
                lab = rownames(resHetMTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect plus interaction',
                subtitle = 'HetMT vs NTg Males (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resHetMTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resHetWTvsNTgMale.png", width = 866, height = 553)
EnhancedVolcano(resHetWTvsNTg,
                lab = rownames(resHetWTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect plus interaction',
                subtitle = 'HetWT vs NTg Males (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resHetWTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resIntHomMTvsNTgMale.png", width = 866, height = 553)
EnhancedVolcano(resIntHomMTvsNTg,
                lab = rownames(resIntHomMTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Interaction',
                subtitle = 'HomMT vs NTg Males (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resIntHomMTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resIntHetMTvsNTgMale.png", width = 866, height = 553)
EnhancedVolcano(resIntHetMTvsNTg,
                lab = rownames(resIntHetMTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Interaction',
                subtitle = 'HetMT vs NTg Males (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resIntHetMTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resIntHetWTvsNTgMale.png", width = 866, height = 553)
EnhancedVolcano(resIntHetWTvsNTg,
                lab = rownames(resIntHetWTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Interaction',
                subtitle = 'HetWT vs NTg Males (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resIntHetWTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resIntHomMTvsHetMTMale.png", width = 866, height = 553)
EnhancedVolcano(resIntHomMTvsHetMT,
                lab = rownames(resIntHomMTvsHetMT),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Interaction',
                subtitle = 'HomMT vs HetMT Males (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resIntHomMTvsHetMT[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resIntHomMTvsHetWTMale.png", width = 866, height = 553)
EnhancedVolcano(resIntHomMTvsHetWT,
                lab = rownames(resIntHomMTvsHetWT),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Interaction',
                subtitle = 'HomMT vs HetWT Males (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resIntHomMTvsHetWT[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resIntHetMTvsHetWTMale.png", width = 866, height = 553)
EnhancedVolcano(resIntHetMTvsHetWT,
                lab = rownames(resIntHetMTvsHetWT),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Interaction',
                subtitle = 'HetMT vs HetWT Males (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resIntHetMTvsHetWT[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()



#Females
char_data <- read.csv("counts_miRNA_F.txt", header=TRUE, sep="\t",row.names="miRNA",stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
num_data <- data.frame(data.matrix(char_data),check.names=FALSE) #finds numerical data.
numeric_columns <- sapply(num_data,function(x){mean(as.numeric(is.na(x)))<0.5}) #converts numeric data to numeric type (assuming columns with less than 50% NAs upon converting to numeric are indeed numeric).
cts <- data.frame(num_data[,numeric_columns], char_data[,!numeric_columns],check.names=FALSE) #merges character data and numeric data (note: it rearranges the column order).
coldata <- read.csv("annotationsF.csv",row.names=1)
coldata <- coldata[,c("Group","Sex","Genotype","Treatment")]
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


#Volcano plots
png("resMainFemale.png", width = 866, height = 553)
EnhancedVolcano(resMain,
                lab = rownames(resMain),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect of Treatment',
                subtitle = 'LPS vs Vehicle in NTg Females',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resMain[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resGenoHomMTvsNTgFemale.png", width = 866, height = 553)
EnhancedVolcano(resGenoHomMTvsNTg,
                lab = rownames(resGenoHomMTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect of Genotype',
                subtitle = 'HomMT vs NTg Females (no LPS)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resGenoHomMTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resGenoHetMTvsNTgFemale.png", width = 866, height = 553)
EnhancedVolcano(resGenoHetMTvsNTg,
                lab = rownames(resGenoHetMTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect of Genotype',
                subtitle = 'HetMT vs NTg Females (no LPS)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resGenoHetMTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resGenoHetWTvsNTgFemale.png", width = 866, height = 553)
EnhancedVolcano(resGenoHetMTvsNTg,
                lab = rownames(resGenoHetWTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect of Genotype',
                subtitle = 'HetWT vs NTg Females (no LPS)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resGenoHetWTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resHomMTvsNTgFemale.png", width = 866, height = 553)
EnhancedVolcano(resHomMTvsNTg,
                lab = rownames(resHomMTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect plus interaction',
                subtitle = 'HomMT vs NTg Females (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resHomMTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resHetMTvsNTgFemale.png", width = 866, height = 553)
EnhancedVolcano(resHetMTvsNTg,
                lab = rownames(resHetMTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect plus interaction',
                subtitle = 'HetMT vs NTg Females (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resHetMTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resHetWTvsNTgFemale.png", width = 866, height = 553)
EnhancedVolcano(resHetWTvsNTg,
                lab = rownames(resHetWTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Main effect plus interaction',
                subtitle = 'HetWT vs NTg Females (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resHetWTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resIntHomMTvsNTgFemale.png", width = 866, height = 553)
EnhancedVolcano(resIntHomMTvsNTg,
                lab = rownames(resIntHomMTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Interaction',
                subtitle = 'HomMT vs NTg Females (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resIntHomMTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resIntHetMTvsNTgFemale.png", width = 866, height = 553)
EnhancedVolcano(resIntHetMTvsNTg,
                lab = rownames(resIntHetMTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Interaction',
                subtitle = 'HetMT vs NTg Females (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resIntHetMTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resIntHetWTvsNTgFemale.png", width = 866, height = 553)
EnhancedVolcano(resIntHetWTvsNTg,
                lab = rownames(resIntHetWTvsNTg),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Interaction',
                subtitle = 'HetWT vs NTg Females (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resIntHetWTvsNTg[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resIntHomMTvsHetMTFemale.png", width = 866, height = 553)
EnhancedVolcano(resIntHomMTvsHetMT,
                lab = rownames(resIntHomMTvsHetMT),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Interaction',
                subtitle = 'HomMT vs HetMT Females (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resIntHomMTvsHetMT[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resIntHomMTvsHetWTFemale.png", width = 866, height = 553)
EnhancedVolcano(resIntHomMTvsHetWT,
                lab = rownames(resIntHomMTvsHetWT),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Interaction',
                subtitle = 'HomMT vs HetWT Females (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resIntHomMTvsHetWT[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()


png("resIntHetMTvsHetWTFemale.png", width = 866, height = 553)
EnhancedVolcano(resIntHetMTvsHetWT,
                lab = rownames(resIntHetMTvsHetWT),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Interaction',
                subtitle = 'HetMT vs HetWT Females (LPS vs Vehicle)',
                caption = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 3.0,
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                col = c("grey30", "forestgreen", "royalblue", "purple"),
                ylim = c(0, max(-log10(resIntHetMTvsHetWT[['pvalue']]), na.rm = TRUE)),
                legendPosition = 'right',
                drawConnectors = TRUE,
                legendLabels=c('Not significant','Log2FC>|0.58|','p<0.05',
                               'p<0.05 & Log2FC>|0.58|'))
dev.off()