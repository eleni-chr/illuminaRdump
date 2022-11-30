#Both
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

#Count plots
d <- plotCounts(dds, "mmu-miR-128-3p", intgroup=c("Treatment","Genotype","Sex","Pairs"), returnData=TRUE)
library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
ggplot(d, aes(x=Treatment, y=count)) +
  geom_point() +
  geom_line(aes(group = Pairs, color="red"), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), strip.text = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16), plot.title = element_text(size = 16)) +
  ggtitle("mmu-miR-128-3p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-128-3p.svg", width=7, height=5)

d <- plotCounts(dds, "mmu-miR-32-5p", intgroup=c("Treatment","Genotype","Sex","Pairs"), returnData=TRUE)
library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
ggplot(d, aes(x=Treatment, y=count)) +
  geom_point() +
  geom_line(aes(group = Pairs, color="red"), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), strip.text = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16), plot.title = element_text(size = 16)) +
  ggtitle("mmu-miR-32-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-32-5p.svg", width=7, height=5)

d <- plotCounts(dds, "mmu-miR-326-3p", intgroup=c("Treatment","Genotype","Sex","Pairs"), returnData=TRUE)
library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
ggplot(d, aes(x=Treatment, y=count)) +
  geom_point() +
  geom_line(aes(group = Pairs, color="red"), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), strip.text = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16), plot.title = element_text(size = 16)) +
  ggtitle("mmu-miR-326-3p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-326-3p.svg", width=7, height=5)

d <- plotCounts(dds, "mmu-let-7b-5p", intgroup=c("Treatment","Genotype","Sex","Pairs"), returnData=TRUE)
library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
ggplot(d, aes(x=Treatment, y=count)) +
  geom_point() +
  geom_line(aes(group = Pairs, color="red"), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), strip.text = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16), plot.title = element_text(size = 16)) +
  ggtitle("mmu-let-7b-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-let-7b-5p.svg", width=7, height=5)

d <- plotCounts(dds, "mmu-miR-28a-3p", intgroup=c("Treatment","Genotype","Sex","Pairs"), returnData=TRUE)
library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
ggplot(d, aes(x=Treatment, y=count)) +
  geom_point() +
  geom_line(aes(group = Pairs, color="red"), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), strip.text = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16), plot.title = element_text(size = 16)) +
  ggtitle("mmu-miR-28a-3p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-28a-3p.svg", width=7, height=5)

d <- plotCounts(dds, "mmu-miR-187-5p", intgroup=c("Treatment","Genotype","Sex","Pairs"), returnData=TRUE)
library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
ggplot(d, aes(x=Treatment, y=count)) +
  geom_point() +
  geom_line(aes(group = Pairs, color="red"), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), strip.text = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16), plot.title = element_text(size = 16)) +
  ggtitle("mmu-miR-187-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-187-5p.svg", width=7, height=5)

d <- plotCounts(dds, "mmu-miR-190a-5p", intgroup=c("Treatment","Genotype","Sex","Pairs"), returnData=TRUE)
library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
ggplot(d, aes(x=Treatment, y=count)) +
  geom_point() +
  geom_line(aes(group = Pairs, color="red"), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), strip.text = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16), plot.title = element_text(size = 16)) +
  ggtitle("mmu-miR-190a-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-190a-5p.svg", width=7, height=5)

d <- plotCounts(dds, "mmu-miR-221-5p", intgroup=c("Treatment","Genotype","Sex","Pairs"), returnData=TRUE)
library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
ggplot(d, aes(x=Treatment, y=count)) +
  geom_point() +
  geom_line(aes(group = Pairs, color="red"), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), strip.text = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16), plot.title = element_text(size = 16)) +
  ggtitle("mmu-miR-221-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-221-5p.svg", width=7, height=5)