install.packages("fitdistrplus")
library(fitdistrplus)

install.packages("logspline")
library(logspline)

install.packages("tidyr")
library(tidyr)

qPCR_data <- read.csv("mmu-miR-128-3p_M.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
qPCR_data <- qPCR_data %>% drop_na()
x <- qPCR_data$DCq
descdist(x, discrete = FALSE)
fit.norm <- fitdist(x, "norm")
plot(fit.norm)
fit.norm$aic

qPCR_data <- read.csv("mmu-miR-128-3p_F.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
qPCR_data <- qPCR_data %>% drop_na()
x <- qPCR_data$DCq
fit.norm <- fitdist(x, "norm")
plot(fit.norm)

qPCR_data <- read.csv("mmu-miR-326-3p_M.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
qPCR_data <- qPCR_data %>% drop_na()
x <- qPCR_data$DCq
fit.norm <- fitdist(x, "norm")
plot(fit.norm)

qPCR_data <- read.csv("mmu-miR-326-3p_F.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
qPCR_data <- qPCR_data %>% drop_na()
x <- qPCR_data$DCq
fit.norm <- fitdist(x, "norm")
plot(fit.norm)

qPCR_data <- read.csv("mmu-miR-221-5p_M.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
qPCR_data <- qPCR_data %>% drop_na()
x <- qPCR_data$DCq
fit.norm <- fitdist(x, "norm")
plot(fit.norm)

qPCR_data <- read.csv("mmu-miR-221-5p_F.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
qPCR_data <- qPCR_data %>% drop_na()
x <- qPCR_data$DCq
fit.norm <- fitdist(x, "norm")
plot(fit.norm)

qPCR_data <- read.csv("hsa-let-7b-5p_M.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
qPCR_data <- qPCR_data %>% drop_na()
x <- qPCR_data$DCq
fit.norm <- fitdist(x, "norm")
plot(fit.norm)

qPCR_data <- read.csv("hsa-let-7b-5p_F.csv", header=TRUE, sep=",",row.names=NULL,stringsAsFactors=FALSE,check.names=FALSE) #imports all data as character type.
qPCR_data <- qPCR_data %>% drop_na()
x <- qPCR_data$DCq
fit.norm <- fitdist(x, "norm")
plot(fit.norm)