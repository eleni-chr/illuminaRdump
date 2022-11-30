d <- plotCounts(dds, gene=which.max(resHomMTvsNTg$log2FoldChange), intgroup="Genotype", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=Genotype, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))



d <- plotCounts(dds, gene=which.max(resHomMTvsNTg$log2FoldChange), intgroup=c("Treatment","Genotype"), 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=interaction(Genotype,Treatment), y=count, color=Genotype, shape=Treatment)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))




d <- plotCounts(dds, "mmu-miR-344f-5p", intgroup=c("Treatment","Genotype"), returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=interaction(Genotype,Treatment), y=count, color=Genotype, shape=Treatment)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


d <- plotCounts(dds, "mmu-miR-182-5p", intgroup=c("Treatment","Genotype"), returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=interaction(Genotype,Treatment), y=count, color=Genotype, shape=Treatment)) + 
  geom_point(position=position_jitter(w=0.1,h=0))


d <- plotCounts(dds, "mmu-miR-182-5p", intgroup=c("Treatment","Genotype","Sex"), returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=interaction(Genotype,Treatment), y=count, color=Genotype, shape=Treatment)) + 
  geom_point(position=position_jitter(w=0.1,h=0))


d <- plotCounts(dds, "mmu-miR-182-5p", intgroup=c("Treatment","Genotype","Sex"), returnData=TRUE)
library("ggplot2")
p <- ggplot(d, aes(x=interaction(Genotype,Treatment), y=count, color=Genotype, shape=Treatment)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
               facet_wrap(~ Sex)
p + ggtitle("mmu-miR-182-5p") + xlab("Interaction") + ylab("Normalised counts")


d <- plotCounts(dds, "mmu-miR-182-5p", intgroup=c("Treatment","Genotype","Sex"), returnData=TRUE)
library("ggplot2")
p <- ggplot(d, aes(x=interaction(Genotype,Treatment), y=count, color=Genotype, shape=Treatment)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  facet_wrap(~ Sex) +
  theme(axis.text.x = element_text(angle = 90, hjust=1))
p + ggtitle("mmu-miR-182-5p") + xlab("Interaction(Genotype,Treatment)") + ylab("Normalised counts")


d <- plotCounts(dds, "mmu-miR-182-5p", intgroup=c("Treatment","Genotype","Sex"), returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=Genotype, y=count, color=Genotype, shape=Treatment)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  facet_grid(Sex ~ Treatment) +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ggtitle("mmu-miR-182-5p") + xlab("Genotype") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))


d <- plotCounts(dds, "mmu-miR-128-3p", intgroup=c("Treatment","Genotype","Sex","Pairs"), returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=Treatment, y=count)) +
  geom_point() +
  geom_line(aes(group = Pairs, color="red"), show.legend=FALSE) +
  facet_grid(Sex ~ Genotype) +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ggtitle("mmu-miR-128-3p") + xlab("Treatment") + ylab("Normalised counts") +
  theme(plot.title = element_text(hjust=0.5))