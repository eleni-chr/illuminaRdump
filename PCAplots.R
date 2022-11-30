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

#Transform the data
ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

#PCA plots
#ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Genotype, label=Sex)) +
#  geom_point(size=3) +
#  geom_label_repel(size=3) +
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#  coord_fixed()

library("ggplot2")
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
image=ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Genotype)) +
  geom_point(size=3) +
  theme(axis.title = element_text(size = 16)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
ggsave(file="PCA_vsd_Male.svg", plot=image, width=5, height=5)


library("ggplot2")
pcaData <- plotPCA(ntd, intgroup=c("Treatment", "Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
image=ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Genotype)) +
  geom_point(size=3) +
  theme(axis.title = element_text(size = 16)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
ggsave(file="PCA_ntd_Male.svg", plot=image, width=5, height=5)


library("ggplot2")
pcaData <- plotPCA(rld, intgroup=c("Treatment", "Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
image=ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Genotype)) +
  geom_point(size=3) +
  theme(axis.title = element_text(size = 16)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
ggsave(file="PCA_rld_Male.svg", plot=image, width=5, height=5)


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

#Transform the data
ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

#PCA plots
library("ggplot2")
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
image=ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Genotype)) +
  geom_point(size=3) +
  theme(axis.title = element_text(size = 16)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
ggsave(file="PCA_vsd_Female.svg", plot=image, width=5, height=5)


library("ggplot2")
pcaData <- plotPCA(ntd, intgroup=c("Treatment", "Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
image=ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Genotype)) +
  geom_point(size=3) +
  theme(axis.title = element_text(size = 16)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
ggsave(file="PCA_ntd_Female.svg", plot=image, width=5, height=5)


library("ggplot2")
pcaData <- plotPCA(rld, intgroup=c("Treatment", "Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
image=ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Genotype)) +
  geom_point(size=3) +
  theme(axis.title = element_text(size = 16)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
ggsave(file="PCA_rld_Female.svg", plot=image, width=5, height=5)
