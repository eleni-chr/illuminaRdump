qPCR_data <- read.csv("mmu-miR-128-3p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.

library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
qPCR_data$Genotype <- factor(qPCR_data$Genotype, levels=c("NTg", "HetWT", "HomWT", "HetMT", "HomMT"))
ggplot(qPCR_data, aes(x=factor(Treatment, level=treatment_order), y=DCq, color=Type)) +
  scale_color_manual(values = c("S"="magenta3", "V"="orange3", "Mean"="black"), labels=c("Seq", "Val", "Mean")) +
  geom_point() + ylim(-6,3) +
  geom_line(aes(group = Pairs, color=Type), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-128-3p") + xlab("") + ylab("ΔCq") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="DCq_mmu-miR-128-3p.svg", width=7, height=5)


qPCR_data <- read.csv("mmu-miR-221-5p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.

library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
qPCR_data$Genotype <- factor(qPCR_data$Genotype, levels=c("NTg", "HetWT", "HomWT", "HetMT", "HomMT"))
ggplot(qPCR_data, aes(x=factor(Treatment, level=treatment_order), y=DCq, color=Type)) +
  scale_color_manual(values = c("S"="magenta3", "V"="orange3", "Mean"="black"), labels=c("Seq", "Val", "Mean")) +
  geom_point() + ylim(-6,3) +
  geom_line(aes(group = Pairs, color=Type), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-221-5p") + xlab("") + ylab("ΔCq") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="DCq_mmu-miR-221-5p.svg", width=7, height=5)


qPCR_data <- read.csv("mmu-miR-326-3p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.

library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
qPCR_data$Genotype <- factor(qPCR_data$Genotype, levels=c("NTg", "HetWT", "HomWT", "HetMT", "HomMT"))
ggplot(qPCR_data, aes(x=factor(Treatment, level=treatment_order), y=DCq, color=Type)) +
  scale_color_manual(values = c("S"="magenta3", "V"="orange3", "Mean"="black"), labels=c("Seq", "Val", "Mean")) +
  geom_point() + ylim(-9,1) +
  geom_line(aes(group = Pairs, color=Type), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-326-3p") + xlab("") + ylab("ΔCq") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="DCq_mmu-miR-326-3p.svg", width=7, height=5)


qPCR_data <- read.csv("hsa-let-7b-5p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.

library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
qPCR_data$Genotype <- factor(qPCR_data$Genotype, levels=c("NTg", "HetWT", "HomWT", "HetMT", "HomMT"))
ggplot(qPCR_data, aes(x=factor(Treatment, level=treatment_order), y=DCq, color=Type)) +
  scale_color_manual(values = c("S"="magenta3", "V"="orange3", "Mean"="black"), labels=c("Seq", "Val", "Mean")) +
  geom_point() + ylim(-7,5) +
  geom_line(aes(group = Pairs, color=Type), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("hsa-let-7b-5p") + xlab("") + ylab("ΔCq") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="DCq_hsa-let7b-5p.svg", width=7, height=5)
