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

#Transformations
ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

#Heatmap of the count matrix
library("pheatmap")
select <- order(rowVars(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20] #select top 20 miRNAs with highest variance across samples
df <- as.data.frame(colData(dds)[,c("Genotype","Treatment")])

image=pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_colnames=FALSE,
        cluster_cols=TRUE, annotation_col=df)
ggsave(file="Heatmap_rld_Male.svg", plot=image, width=7, height=5)

image=pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_colnames=FALSE,
               cluster_cols=TRUE, annotation_col=df)
ggsave(file="Heatmap_vsd_Male.svg", plot=image, width=7, height=5)

image=pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_colnames=FALSE,
               cluster_cols=TRUE, annotation_col=df)
ggsave(file="Heatmap_ntd_Male.svg", plot=image, width=7, height=5)

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

#Transformations
ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

#Heatmap of the count matrix
library("pheatmap")
select <- order(rowVars(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20] #select top 20 miRNAs with highest variance across samples
df <- as.data.frame(colData(dds)[,c("Genotype","Treatment")])

image=pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_colnames=FALSE,
               cluster_cols=TRUE, annotation_col=df)
ggsave(file="Heatmap_rld_Female.svg", plot=image, width=7, height=5)

image=pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_colnames=FALSE,
               cluster_cols=TRUE, annotation_col=df)
ggsave(file="Heatmap_vsd_Female.svg", plot=image, width=7, height=5)

image=pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_colnames=FALSE,
               cluster_cols=TRUE, annotation_col=df)
ggsave(file="Heatmap_ntd_Female.svg", plot=image, width=7, height=5)


#Both sexes
char_data <- read.csv("counts_miRNA.txt", header=TRUE, sep="\t",row.names="miRNA",stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
num_data <- data.frame(data.matrix(char_data),check.names=FALSE) #finds numerical data.
numeric_columns <- sapply(num_data,function(x){mean(as.numeric(is.na(x)))<0.5}) #converts numeric data to numeric type (assuming columns with less than 50% NAs upon converting to numeric are indeed numeric).
cts <- data.frame(num_data[,numeric_columns], char_data[,!numeric_columns],check.names=FALSE) #merges character data and numeric data (note: it rearranges the column order).
coldata <- read.csv("annotations.csv",row.names=1)
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

#Transformations
ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

#Heatmap of the count matrix
library("pheatmap")
select <- order(rowVars(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20] #select top 20 miRNAs with highest variance across samples
df <- as.data.frame(colData(dds)[,c("Genotype","Treatment","Sex")])

image=pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_colnames=FALSE,
               cluster_cols=TRUE, annotation_col=df)
ggsave(file="Heatmap_rld_Both.svg", plot=image, width=7, height=5)

image=pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_colnames=FALSE,
               cluster_cols=TRUE, annotation_col=df)
ggsave(file="Heatmap_vsd_Both.svg", plot=image, width=7, height=5)

image=pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_colnames=FALSE,
               cluster_cols=TRUE, annotation_col=df)
ggsave(file="Heatmap_ntd_Both.svg", plot=image, width=7, height=5)