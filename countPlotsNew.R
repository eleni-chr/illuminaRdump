d <- read.csv("mmu-miR-128-3p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.

library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
d$Genotype <- factor(d$Genotype, levels=c("NTg", "HetWT", "HetMT", "HomMT"))
ggplot(d, aes(x=factor(Treatment, level=treatment_order), y=count, color=Type)) +
  scale_color_manual(values = c("Sample"="magenta3", "Mean"="black"), labels=c("Replicates", "Mean")) +
  geom_point() + ylim(0,200) +
  geom_line(aes(group = Pairs), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-128-3p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-128-3p.svg", width=7, height=5)


d <- read.csv("mmu-miR-28a-3p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.

library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
d$Genotype <- factor(d$Genotype, levels=c("NTg", "HetWT", "HetMT", "HomMT"))
ggplot(d, aes(x=factor(Treatment, level=treatment_order), y=count, color=Type)) +
  scale_color_manual(values = c("Sample"="magenta3", "Mean"="black"), labels=c("Replicates", "Mean")) +
  geom_point() +
  geom_line(aes(group = Pairs), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-28a-3p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-28a-3p.svg", width=7, height=5)


d <- read.csv("mmu-miR-32-5p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.

library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
d$Genotype <- factor(d$Genotype, levels=c("NTg", "HetWT", "HetMT", "HomMT"))
ggplot(d, aes(x=factor(Treatment, level=treatment_order), y=count, color=Type)) +
  scale_color_manual(values = c("Sample"="magenta3", "Mean"="black"), labels=c("Replicates", "Mean")) +
  geom_point() + ylim(0,200) +
  geom_line(aes(group = Pairs), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-32-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-32-5p.svg", width=7, height=5)


d <- read.csv("mmu-miR-187-5p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.

library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
d$Genotype <- factor(d$Genotype, levels=c("NTg", "HetWT", "HetMT", "HomMT"))
ggplot(d, aes(x=factor(Treatment, level=treatment_order), y=count, color=Type)) +
  scale_color_manual(values = c("Sample"="magenta3", "Mean"="black"), labels=c("Replicates", "Mean")) +
  geom_point() + ylim(0,200) +
  geom_line(aes(group = Pairs), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-187-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-187-5p.svg", width=7, height=5)


d <- read.csv("mmu-miR-190a-5p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.

library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
d$Genotype <- factor(d$Genotype, levels=c("NTg", "HetWT", "HetMT", "HomMT"))
ggplot(d, aes(x=factor(Treatment, level=treatment_order), y=count, color=Type)) +
  scale_color_manual(values = c("Sample"="magenta3", "Mean"="black"), labels=c("Replicates", "Mean")) +
  geom_point() + ylim(0,200) +
  geom_line(aes(group = Pairs), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-190a-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-190a-5p.svg", width=7, height=5)


d <- read.csv("mmu-miR-221-5p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.

library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
d$Genotype <- factor(d$Genotype, levels=c("NTg", "HetWT", "HetMT", "HomMT"))
ggplot(d, aes(x=factor(Treatment, level=treatment_order), y=count, color=Type)) +
  scale_color_manual(values = c("Sample"="magenta3", "Mean"="black"), labels=c("Replicates", "Mean")) +
  geom_point() + ylim(0,200) +
  geom_line(aes(group = Pairs), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-221-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-221-5p.svg", width=7, height=5)


d <- read.csv("mmu-miR-326-3p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.

library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
d$Genotype <- factor(d$Genotype, levels=c("NTg", "HetWT", "HetMT", "HomMT"))
ggplot(d, aes(x=factor(Treatment, level=treatment_order), y=count, color=Type)) +
  scale_color_manual(values = c("Sample"="magenta3", "Mean"="black"), labels=c("Replicates", "Mean")) +
  geom_point() + ylim(0,200) +
  geom_line(aes(group = Pairs), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-326-3p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-326-3p.svg", width=7, height=5)


d <- read.csv("mmu-let-7b-5p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.

library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
d$Genotype <- factor(d$Genotype, levels=c("NTg", "HetWT", "HetMT", "HomMT"))
ggplot(d, aes(x=factor(Treatment, level=treatment_order), y=count, color=Type)) +
  scale_color_manual(values = c("Sample"="magenta3", "Mean"="black"), labels=c("Replicates", "Mean")) +
  geom_point() + ylim(0,200) +
  geom_line(aes(group = Pairs), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-let-7b-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-let-7b-5p.svg", width=7, height=5)




#NEW

d <- read.csv("mmu-miR-425-5p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
d$Genotype <- factor(d$Genotype, levels=c("NTg", "HetWT", "HetMT", "HomMT"))
ggplot(d, aes(x=factor(Treatment, level=treatment_order), y=count, color=Type)) +
  scale_color_manual(values = c("Sample"="magenta3", "Mean"="black"), labels=c("Replicates", "Mean")) +
  geom_point() + ylim(0,200) +
  geom_line(aes(group = Pairs), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-425-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-425-5p.svg", width=7, height=5)


d <- read.csv("mmu-miR-339-5p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
d$Genotype <- factor(d$Genotype, levels=c("NTg", "HetWT", "HetMT", "HomMT"))
ggplot(d, aes(x=factor(Treatment, level=treatment_order), y=count, color=Type)) +
  scale_color_manual(values = c("Sample"="magenta3", "Mean"="black"), labels=c("Replicates", "Mean")) +
  geom_point() + ylim(0,200) +
  geom_line(aes(group = Pairs), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-339-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-339-5p.svg", width=7, height=5)


d <- read.csv("mmu-miR-92a-3p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
d$Genotype <- factor(d$Genotype, levels=c("NTg", "HetWT", "HetMT", "HomMT"))
ggplot(d, aes(x=factor(Treatment, level=treatment_order), y=count, color=Type)) +
  scale_color_manual(values = c("Sample"="magenta3", "Mean"="black"), labels=c("Replicates", "Mean")) +
  geom_point() + ylim(0,200) +
  geom_line(aes(group = Pairs), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-92a-3p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-92a-3p.svg", width=7, height=5)


d <- read.csv("mmu-miR-142a-5p.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
library("ggplot2")
Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")
treatment_order <- c('Vehicle', 'LPS')
d$Genotype <- factor(d$Genotype, levels=c("NTg", "HetWT", "HetMT", "HomMT"))
ggplot(d, aes(x=factor(Treatment, level=treatment_order), y=count, color=Type)) +
  scale_color_manual(values = c("Sample"="magenta3", "Mean"="black"), labels=c("Replicates", "Mean")) +
  geom_point() + ylim(0,200) +
  geom_line(aes(group = Pairs), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype, labeller=labeller(Sex=Sex.labs)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5), 
        strip.text = element_text(size = 16), axis.text = element_text(size = 16),
        axis.title = element_text(size = 16), plot.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title=element_blank()) +
  ggtitle("mmu-miR-142a-5p") + xlab("") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(file="counts_mmu-miR-142a-5p.svg", width=7, height=5)