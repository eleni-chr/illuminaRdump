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
image=EnhancedVolcano(resMain,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'LPS vs Vehicle in NTg Males',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resMain[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resMainMale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resGenoHomMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HomMT vs NTg Males (no LPS)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resGenoHomMTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resGenoHomMTvsNTgMale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resGenoHetMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetMT vs NTg Males (no LPS)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resGenoHetMTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resGenoHetMTvsNTgMale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resGenoHetMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetWT vs NTg Males (no LPS)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resGenoHetWTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resGenoHetWTvsNTgMale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resHomMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HomMT vs NTg Males (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resHomMTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resHomMTvsNTgMale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resHetMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetMT vs NTg Males (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resHetMTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resHetMTvsNTgMale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resHetWTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetWT vs NTg Males (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resHetWTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resHetWTvsNTgMale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resIntHomMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HomMT vs NTg Males (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resIntHomMTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resIntHomMTvsNTgMale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resIntHetMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetMT vs NTg Males (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resIntHetMTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resIntHetMTvsNTgMale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resIntHetWTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetWT vs NTg Males (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resIntHetWTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resIntHetWTvsNTgMale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resIntHomMTvsHetMT,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HomMT vs HetMT Males (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resIntHomMTvsHetMT[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resIntHomMTvsHetMTMale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resIntHomMTvsHetWT,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HomMT vs HetWT Males (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resIntHomMTvsHetWT[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resIntHomMTvsHetWTMale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resIntHetMTvsHetWT,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetMT vs HetWT Males (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resIntHetMTvsHetWT[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resIntHetMTvsHetWTMale.svg", plot=image, width=5, height=5)


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
image=EnhancedVolcano(resMain,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'LPS vs Vehicle in NTg Females',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resMain[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resMainFemale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resGenoHomMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HomMT vs NTg Females (no LPS)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resGenoHomMTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resGenoHomMTvsNTgFemale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resGenoHetMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetMT vs NTg Females (no LPS)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resGenoHetMTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resGenoHetMTvsNTgFemale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resGenoHetMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetWT vs NTg Females (no LPS)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resGenoHetWTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resGenoHetWTvsNTgFemale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resHomMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HomMT vs NTg Females (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resHomMTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resHomMTvsNTgFemale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resHetMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetMT vs NTg Females (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resHetMTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resHetMTvsNTgFemale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resHetWTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetWT vs NTg Females (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resHetWTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resHetWTvsNTgFemale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resIntHomMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HomMT vs NTg Females (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resIntHomMTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resIntHomMTvsNTgFemale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resIntHetMTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetMT vs NTg Females (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resIntHetMTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resIntHetMTvsNTgFemale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resIntHetWTvsNTg,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetWT vs NTg Females (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resIntHetWTvsNTg[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resIntHetWTvsNTgFemale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resIntHomMTvsHetMT,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HomMT vs HetMT Females (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resIntHomMTvsHetMT[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resIntHomMTvsHetMTFemale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resIntHomMTvsHetWT,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HomMT vs HetWT Females (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resIntHomMTvsHetWT[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resIntHomMTvsHetWTFemale.svg", plot=image, width=5, height=5)


image=EnhancedVolcano(resIntHetMTvsHetWT,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = NULL,
                      subtitle = 'HetMT vs HetWT Females (LPS vs Vehicle)',
                      subtitleLabSize = 18,
                      caption = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = 3.0,
                      gridlines.major=FALSE,
                      gridlines.minor=FALSE,
                      col = c("grey30", "forestgreen", "royalblue", "purple"),
                      ylim = c(0, max(-log10(resIntHetMTvsHetWT[['pvalue']]), na.rm = TRUE)),
                      legendPosition = "none")
ggsave(file="resIntHetMTvsHetWTFemale.svg", plot=image, width=5, height=5)
