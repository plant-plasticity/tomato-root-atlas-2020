#Updated on 11/10/18
#Updated on 12/14/18
#Updated on 04/24/19
#install.packages("matrixStats", lib="/share/brady/people/AlexMason/R_packages")
#require(matrixStats, lib.loc="/share/brady/people/AlexMason/R_packages")
#E:/TEMP/DiffPeak_betaTesting
require(matrixStats)
#working_directories = c("COR","EN","EP","EXO","MCO","MZ","PH","V","WOX","X35S","XY")
#working_directories = c("EXO")

#COR_reps=c("K06","O07","R07","S06")
#EN_reps=c("L05","N06","R06","S09")
#EP_reps=c("J05","N04","Q07","S07")
#EXO_reps=c("J06","M05","P04","S05")
#MCO_reps=c("I08","K05","M06","P05")
#MZ_reps=c("L07","O09","R09","S03")
#PH_reps=c("I06","L06","M07","Q05")
#V_reps=c("I05","J07","M08","R08")
#WOX_reps=c("K07","O08","P06","S08")
#X35S_reps=c("J08","O06","Q08","S02")
#XY_reps=c("I07","K08","N05","Q06")

#####
#COR#
#####
setwd("E:/TEMP/DiffPeak_betaTesting/COR")
#Make counts matrix of all the samples you're going to look across
Rep1=read.delim("COR_K06.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep2=read.delim("COR_O07.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep3=read.delim("COR_R07.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep4=read.delim("COR_S06.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)

#for later
Rep1[,1] = as.character(Rep1[,1])

#cbind(..., deparse.level = 1)
counts.data= cbind(Rep1$V4, Rep2$V4, Rep3$V4, Rep4$V4)
counts.data[is.na(counts.data)] <- 0
head(counts.data)

counts.data[,1] = as.numeric(as.character(counts.data[,1]))
counts.data[,2] = as.numeric(as.character(counts.data[,2]))
counts.data[,3] = as.numeric(as.character(counts.data[,3]))
counts.data[,4] = as.numeric(as.character(counts.data[,4]))

#A large number to normalize cut counts with
total = 10000000

#If data not there, mark as 0
counts.data[which(is.nan(counts.data[,1])),1] = 0
counts.data[which(is.nan(counts.data[,2])),2] = 0
counts.data[which(is.nan(counts.data[,3])),3] = 0
counts.data[which(is.nan(counts.data[,4])),4] = 0


#Normalize within sample cut counts based on the sum of cuts across all the urepTHSs
counts.data[,1] = total * counts.data[,1] / sum(counts.data[,1])
counts.data[,2] = total * counts.data[,2] / sum(counts.data[,2])
counts.data[,3] = total * counts.data[,3] / sum(counts.data[,3])
counts.data[,4] = total * counts.data[,4] / sum(counts.data[,4])

##How I determined what reps were outliers
setwd("E:/BRADY_LAB/ATLAS/ATLAS_ATAC_figures/For paper")

pdf(file='scatter_plot_replicates_repUnion_THSs_COR_042419.pdf', height=8, width=22)
par(mfrow=c(1,6))
plot(log(counts.data[,1]), log(counts.data[,2]), pch=20, cex=3, col=adjustcolor("#2E5EAC", alpha.f=0.25))
plot(log(counts.data[,1]), log(counts.data[,3]), pch=20, cex=3, col=adjustcolor("#2E5EAC", alpha.f=0.25))
plot(log(counts.data[,1]), log(counts.data[,4]), pch=20, cex=3, col=adjustcolor("#2E5EAC", alpha.f=0.25))
plot(log(counts.data[,2]), log(counts.data[,3]), pch=20, cex=3, col=adjustcolor("#2E5EAC", alpha.f=0.25))
plot(log(counts.data[,2]), log(counts.data[,4]), pch=20, cex=3, col=adjustcolor("#2E5EAC", alpha.f=0.25))
plot(log(counts.data[,3]), log(counts.data[,4]), pch=20, cex=3, col=adjustcolor("#2E5EAC", alpha.f=0.25))
dev.off()

#needs to go across rows
means = rowMeans(counts.data)
SD = rowSds(counts.data)
COR.CV = (SD/means)*100

pdf(file='scatter_plot_SD_mean_CV_repUnion_THSs_COR_042419.pdf', height=11, width=8)
#par(mfrow=c(1,2))
plot(means, SD, pch=20, cex=3, col=adjustcolor("#2E5EAC", alpha.f=0.25))
abline(fit <- lm(SD ~ means),col="#225ea8", lwd=3)
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
#plot(COR.CV, means, pch=20, cex=3, col=adjustcolor("#2E5EAC", alpha.f=0.25))
dev.off()

#Looks like a CV of 50 maybe a good cutoff for COR
par(mfrow=c(1,2))
plot(log(means), log(SD), pch=20, cex=3, col=adjustcolor("#fb6a4a", alpha.f=0.25))
plot(log(means), log(COR.CV), pch=20, cex=3, col=adjustcolor("#fb6a4a", alpha.f=0.25))


#paste columns in data frame to do stuff with
COR.data= cbind(Rep1$V1, Rep1$V2, Rep1$V3, COR.CV)
COR.data <- as.data.frame(COR.data)
names(COR.data) <- c("chr","start","end","cv")
COR.data$cv = as.numeric(as.character(COR.data$cv))
#find top X% and subset against
top = 15
COR.top=COR.data[COR.data$cv > quantile(COR.data$cv,prob=1-top/100),]
COR.bottom=COR.data[COR.data$cv < quantile(COR.data$cv,prob=1-top/100),]
max(COR.bottom$cv)
[1] 73.32582
#write subsetted data to a bed file
write.table(COR.bottom, file = "COR.masked.25mil.localSize10kb.repTHSs.bed", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#plots


COR.h <- hist(COR.CV, breaks = 15, density = 10, col = "lightgray", xlab = "CV", main = "COR") 
COR.xfit <- seq(min(COR.CV), max(COR.CV), length = 40) 
COR.yfit <- dnorm(COR.xfit, mean = mean(COR.CV), sd = sd(COR.CV)) 
COR.yfit <- COR.yfit * diff(COR.h$mids[1:2]) * length(COR.CV) 
lines(COR.xfit, COR.yfit, col = "black", lwd = 2)
abline(v=max(COR.bottom$cv))
legend(0,4000,legend=paste("CV: ",round(max(COR.bottom$cv), digits=3),sep=''))

####
#EN#
####
setwd("E:/TEMP/DiffPeak_betaTesting/EN")
#Make counts matrix of all the samples you're going to look across
Rep1=read.delim("EN_L05.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep2=read.delim("EN_N06.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep3=read.delim("EN_R06.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep4=read.delim("EN_S09.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)

#for later
Rep1[,1] = as.character(Rep1[,1])

#cbind(..., deparse.level = 1)
counts.data= cbind(Rep1$V4, Rep2$V4, Rep3$V4, Rep4$V4)
counts.data[is.na(counts.data)] <- 0
head(counts.data)

counts.data[,1] = as.numeric(as.character(counts.data[,1]))
counts.data[,2] = as.numeric(as.character(counts.data[,2]))
counts.data[,3] = as.numeric(as.character(counts.data[,3]))
counts.data[,4] = as.numeric(as.character(counts.data[,4]))

#A large number to normalize cut counts with
total = 10000000

#If data not there, mark as 0
counts.data[which(is.nan(counts.data[,1])),1] = 0
counts.data[which(is.nan(counts.data[,2])),2] = 0
counts.data[which(is.nan(counts.data[,3])),3] = 0
counts.data[which(is.nan(counts.data[,4])),4] = 0

#Normalize within sample cut counts based on the sum of cuts across all the urepTHSs
counts.data[,1] = total * counts.data[,1] / sum(counts.data[,1])
counts.data[,2] = total * counts.data[,2] / sum(counts.data[,2])
counts.data[,3] = total * counts.data[,3] / sum(counts.data[,3])
counts.data[,4] = total * counts.data[,4] / sum(counts.data[,4])

##How I determined what reps were outliers
setwd("E:/BRADY_LAB/ATLAS/ATLAS_ATAC_figures/For paper")

pdf(file='scatter_plot_replicates_repUnion_THSs_EN_042419.pdf', height=11, width=8)
par(mfrow=c(1,6))
plot(counts.data[,1], counts.data[,2], pch=20, cex=3, col=adjustcolor("#016E82", alpha.f=0.25))
plot(counts.data[,1], counts.data[,3], pch=20, cex=3, col=adjustcolor("#016E82", alpha.f=0.25))
plot(counts.data[,1], counts.data[,4], pch=20, cex=3, col=adjustcolor("#016E82", alpha.f=0.25))
plot(counts.data[,2], counts.data[,3], pch=20, cex=3, col=adjustcolor("#016E82", alpha.f=0.25))
plot(counts.data[,2], counts.data[,4], pch=20, cex=3, col=adjustcolor("#016E82", alpha.f=0.25))
plot(counts.data[,3], counts.data[,4], pch=20, cex=3, col=adjustcolor("#016E82", alpha.f=0.25))
dev.off()

#needs to go across rows
means = rowMeans(counts.data)
SD = rowSds(counts.data)
EN.CV = (SD/means)*100

pdf(file='scatter_plot_SD_mean_CV_repUnion_THSs_EN_042419.pdf', height=11, width=8)
#par(mfrow=c(1,2))
plot(means, SD, pch=20, cex=3, col=adjustcolor("#016E82", alpha.f=0.25))
abline(fit <- lm(SD ~ means),col="#225ea8", lwd=3)
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
#plot(EN.CV, means, pch=20, cex=3, col=adjustcolor("#016E82", alpha.f=0.25))
dev.off()


#paste columns in data frame to do stuff with
EN.data= cbind(Rep1$V1, Rep1$V2, Rep1$V3, EN.CV)
EN.data <- as.data.frame(EN.data)
names(EN.data) <- c("chr","start","end","cv")
EN.data$cv = as.numeric(as.character(EN.data$cv))
#find top X% and subset against
top = 15
EN.top=EN.data[EN.data$cv > quantile(EN.data$cv,prob=1-top/100),]
EN.bottom=EN.data[EN.data$cv < quantile(EN.data$cv,prob=1-top/100),]
#write subsetted data to a bed file
write.table(EN.bottom, file = "EN.masked.25mil.localSize10kb.repTHSs.bed", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#plots
EN.h <- hist(EN.CV, breaks = 10, density = 10, col = "lightgray", xlab = "CV", main = "EN") 
EN.xfit <- seq(min(EN.CV), max(EN.CV), length = 40) 
EN.yfit <- dnorm(EN.xfit, mean = mean(EN.CV), sd = sd(EN.CV)) 
EN.yfit <- EN.yfit * diff(EN.h$mids[1:2]) * length(EN.CV) 
lines(EN.xfit, EN.yfit, col = "black", lwd = 2)
abline(v=max(EN.bottom$cv))
legend(0,4000,legend=paste("CV: ",round(max(EN.bottom$cv), digits=3),sep=''))

####
#EP#
####
setwd("E:/TEMP/DiffPeak_betaTesting/EP")
#Make counts matrix of all the samples you're going to look across
Rep1=read.delim("EP_J05.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep2=read.delim("EP_N04.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep3=read.delim("EP_Q07.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep4=read.delim("EP_S07.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)

#for later
Rep1[,1] = as.character(Rep1[,1])

#cbind(..., deparse.level = 1)
counts.data= cbind(Rep1$V4, Rep2$V4, Rep3$V4, Rep4$V4)
counts.data[is.na(counts.data)] <- 0
head(counts.data)

counts.data[,1] = as.numeric(as.character(counts.data[,1]))
counts.data[,2] = as.numeric(as.character(counts.data[,2]))
counts.data[,3] = as.numeric(as.character(counts.data[,3]))
counts.data[,4] = as.numeric(as.character(counts.data[,4]))

#A large number to normalize cut counts with
total = 10000000

#If data not there, mark as 0
counts.data[which(is.nan(counts.data[,1])),1] = 0
counts.data[which(is.nan(counts.data[,2])),2] = 0
counts.data[which(is.nan(counts.data[,3])),3] = 0
counts.data[which(is.nan(counts.data[,4])),4] = 0

#Normalize within sample cut counts based on the sum of cuts across all the urepTHSs
counts.data[,1] = total * counts.data[,1] / sum(counts.data[,1])
counts.data[,2] = total * counts.data[,2] / sum(counts.data[,2])
counts.data[,3] = total * counts.data[,3] / sum(counts.data[,3])
counts.data[,4] = total * counts.data[,4] / sum(counts.data[,4])

##How I determined what reps were outliers
setwd("E:/BRADY_LAB/ATLAS/ATLAS_ATAC_figures/For paper")

pdf(file='scatter_plot_replicates_repUnion_THSs_EP_042419.pdf', height=22, width=8)
par(mfrow=c(1,6))
plot(counts.data[,1], counts.data[,2], pch=20, cex=3, col=adjustcolor("#44C3D0", alpha.f=0.25))
plot(counts.data[,1], counts.data[,3], pch=20, cex=3, col=adjustcolor("#44C3D0", alpha.f=0.25))
plot(counts.data[,1], counts.data[,4], pch=20, cex=3, col=adjustcolor("#44C3D0", alpha.f=0.25))
plot(counts.data[,2], counts.data[,3], pch=20, cex=3, col=adjustcolor("#44C3D0", alpha.f=0.25))
plot(counts.data[,2], counts.data[,4], pch=20, cex=3, col=adjustcolor("#44C3D0", alpha.f=0.25))
plot(counts.data[,3], counts.data[,4], pch=20, cex=3, col=adjustcolor("#44C3D0", alpha.f=0.25))
dev.off()

#needs to go across rows
means = rowMeans(counts.data)
SD = rowSds(counts.data)
EP.CV = (SD/means)*100

pdf(file='scatter_plot_SD_mean_CV_repUnion_THSs_EP_042419.pdf', height=11, width=8)
#par(mfrow=c(1,2))
plot(means, SD, pch=20, cex=3, col=adjustcolor("#44C3D0", alpha.f=0.25))
abline(fit <- lm(SD ~ means),col="#44C3D0", lwd=3)
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
#plot(EP.CV, means, pch=20, cex=3, col=adjustcolor("#44C3D0", alpha.f=0.25))
dev.off()

#paste columns in data frame to do stuff with
EP.data= cbind(Rep1$V1, Rep1$V2, Rep1$V3, EP.CV)
EP.data <- as.data.frame(EP.data)
names(EP.data) <- c("chr","start","end","cv")
EP.data$cv = as.numeric(as.character(EP.data$cv))
#find top X% and subset against
top = 15
EP.top=EP.data[EP.data$cv > quantile(EP.data$cv,prob=1-top/100),]
EP.bottom=EP.data[EP.data$cv < quantile(EP.data$cv,prob=1-top/100),]
#write subsetted data to a bed file
write.table(EP.bottom, file = "EP.masked.25mil.localSize10kb.repTHSs.bed", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#PLOTS
EP.h <- hist(EP.CV, breaks = 10, density = 10, col = "lightgray", xlab = "CV", main = "EP") 
EP.xfit <- seq(min(EP.CV), max(EP.CV), length = 40) 
EP.yfit <- dnorm(EP.xfit, mean = mean(EP.CV), sd = sd(EP.CV)) 
EP.yfit <- EP.yfit * diff(EP.h$mids[1:2]) * length(EP.CV) 
lines(EP.xfit, EP.yfit, col = "black", lwd = 2)
abline(v=max(EP.bottom$cv))
legend(0,3000,legend=paste("CV: ",round(max(EP.bottom$cv), digits=3),sep=''))

#####
#EXO#
#####
setwd("E:/TEMP/DiffPeak_betaTesting/EXO")
#Make counts matrix of all the samples you're going to look across
Rep1=read.delim("EXO_J06.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep2=read.delim("EXO_M05.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep3=read.delim("EXO_P04.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep4=read.delim("EXO_S05.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)

#for later
Rep1[,1] = as.character(Rep1[,1])

#cbind(..., deparse.level = 1)
counts.data= cbind(Rep1$V4, Rep2$V4, Rep3$V4, Rep4$V4)
counts.data[is.na(counts.data)] <- 0
head(counts.data)

counts.data[,1] = as.numeric(as.character(counts.data[,1]))
counts.data[,2] = as.numeric(as.character(counts.data[,2]))
counts.data[,3] = as.numeric(as.character(counts.data[,3]))
counts.data[,4] = as.numeric(as.character(counts.data[,4]))

#A large number to normalize cut counts with
total = 10000000

#If data not there, mark as 0
counts.data[which(is.nan(counts.data[,1])),1] = 0
counts.data[which(is.nan(counts.data[,2])),2] = 0
counts.data[which(is.nan(counts.data[,3])),3] = 0
counts.data[which(is.nan(counts.data[,4])),4] = 0

#Normalize within sample cut counts based on the sum of cuts across all the urepTHSs
counts.data[,1] = total * counts.data[,1] / sum(counts.data[,1])
counts.data[,2] = total * counts.data[,2] / sum(counts.data[,2])
counts.data[,3] = total * counts.data[,3] / sum(counts.data[,3])
counts.data[,4] = total * counts.data[,4] / sum(counts.data[,4])

##How I determined what reps were outliers
setwd("E:/BRADY_LAB/ATLAS/ATLAS_ATAC_figures/For paper")

pdf(file='scatter_plot_replicates_repUnion_THSs_EXO_042419.pdf', height=11, width=8)
par(mfrow=c(1,6))
plot(counts.data[,1], counts.data[,2], pch=20, cex=3, col=adjustcolor("#7D2C91", alpha.f=0.25))
plot(counts.data[,1], counts.data[,3], pch=20, cex=3, col=adjustcolor("#7D2C91", alpha.f=0.25))
plot(counts.data[,1], counts.data[,4], pch=20, cex=3, col=adjustcolor("#7D2C91", alpha.f=0.25))
plot(counts.data[,2], counts.data[,3], pch=20, cex=3, col=adjustcolor("#7D2C91", alpha.f=0.25))
plot(counts.data[,2], counts.data[,4], pch=20, cex=3, col=adjustcolor("#7D2C91", alpha.f=0.25))
plot(counts.data[,3], counts.data[,4], pch=20, cex=3, col=adjustcolor("#7D2C91", alpha.f=0.25))
dev.off()


#needs to go across rows
means = rowMeans(counts.data)
SD = rowSds(counts.data)
EXO.CV = (SD/means)*100

pdf(file='scatter_plot_SD_mean_CV_repUnion_THSs_EXO_042419.pdf', height=11, width=8)
#par(mfrow=c(1,2))
plot(means, SD, pch=20, cex=3, col=adjustcolor("#7D2C91", alpha.f=0.25))
abline(fit <- lm(SD ~ means),col="#225ea8", lwd=3)
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
#plot(EXO.CV, means, pch=20, cex=3, col=adjustcolor("#fb6a4a", alpha.f=0.25))
dev.off()


#paste columns in data frame to do stuff with
EXO.data= cbind(Rep1$V1, Rep1$V2, Rep1$V3, EXO.CV)
EXO.data <- as.data.frame(EXO.data)
names(EXO.data) <- c("chr","start","end","cv")
EXO.data$cv = as.numeric(as.character(EXO.data$cv))
#find top X% and subset against
top = 15
EXO.top=EXO.data[EXO.data$cv > quantile(EXO.data$cv,prob=1-top/100),]
EXO.bottom=EXO.data[EXO.data$cv < quantile(EXO.data$cv,prob=1-top/100),]
#write subsetted data to a bed file
write.table(EXO.bottom, file = "EXO.masked.25mil.localSize10kb.repTHSs.bed", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#PLOTS
EXO.h <- hist(EXO.CV, breaks = 10, density = 10, col = "lightgray", xlab = "CV", main = "EXO") 
EXO.xfit <- seq(min(EXO.CV), max(EXO.CV), length = 40) 
EXO.yfit <- dnorm(EXO.xfit, mean = mean(EXO.CV), sd = sd(EXO.CV)) 
EXO.yfit <- EXO.yfit * diff(EXO.h$mids[1:2]) * length(EXO.CV) 
lines(EXO.xfit, EXO.yfit, col = "black", lwd = 2)
abline(v=max(EXO.bottom$cv))
legend(0,5000,legend=paste("CV: ",round(max(EXO.bottom$cv), digits=3),sep=''))

#####
#MCO#
#####
setwd("E:/TEMP/DiffPeak_betaTesting/MCO")
#Make counts matrix of all the samples you're going to look across
Rep1=read.delim("MCO_I08.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep2=read.delim("MCO_K05.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep3=read.delim("MCO_M06.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep4=read.delim("MCO_P05.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)

#for later
Rep1[,1] = as.character(Rep1[,1])

#cbind(..., deparse.level = 1)
counts.data= cbind(Rep1$V4, Rep2$V4, Rep3$V4, Rep4$V4)
counts.data[is.na(counts.data)] <- 0
head(counts.data)

counts.data[,1] = as.numeric(as.character(counts.data[,1]))
counts.data[,2] = as.numeric(as.character(counts.data[,2]))
counts.data[,3] = as.numeric(as.character(counts.data[,3]))
counts.data[,4] = as.numeric(as.character(counts.data[,4]))

#A large number to normalize cut counts with
total = 10000000

#If data not there, mark as 0
counts.data[which(is.nan(counts.data[,1])),1] = 0
counts.data[which(is.nan(counts.data[,2])),2] = 0
counts.data[which(is.nan(counts.data[,3])),3] = 0
counts.data[which(is.nan(counts.data[,4])),4] = 0

#Normalize within sample cut counts based on the sum of cuts across all the urepTHSs
counts.data[,1] = total * counts.data[,1] / sum(counts.data[,1])
counts.data[,2] = total * counts.data[,2] / sum(counts.data[,2])
counts.data[,3] = total * counts.data[,3] / sum(counts.data[,3])
counts.data[,4] = total * counts.data[,4] / sum(counts.data[,4])


##How I determined what reps were outliers
setwd("E:/BRADY_LAB/ATLAS/ATLAS_ATAC_figures/For paper")

pdf(file='scatter_plot_replicates_repUnion_THSs_MCO_042419.pdf', height=11, width=8)
par(mfrow=c(1,6))
plot(counts.data[,1], counts.data[,2], SD, pch=20, cex=3, col=adjustcolor("#CD85B9", alpha.f=0.25))
plot(counts.data[,1], counts.data[,3], SD, pch=20, cex=3, col=adjustcolor("#CD85B9", alpha.f=0.25))
plot(counts.data[,1], counts.data[,4], SD, pch=20, cex=3, col=adjustcolor("#CD85B9", alpha.f=0.25))
plot(counts.data[,2], counts.data[,3], SD, pch=20, cex=3, col=adjustcolor("#CD85B9", alpha.f=0.25))
plot(counts.data[,2], counts.data[,4], SD, pch=20, cex=3, col=adjustcolor("#CD85B9", alpha.f=0.25))
plot(counts.data[,3], counts.data[,4], SD, pch=20, cex=3, col=adjustcolor("#CD85B9", alpha.f=0.25))
dev.off()

#needs to go across rows
means = rowMeans(counts.data)
SD = rowSds(counts.data)
MCO.CV = (SD/means)*100

pdf(file='scatter_plot_SD_mean_CV_repUnion_THSs_MCO_042419.pdf', height=11, width=8)
#par(mfrow=c(1,2))
plot(means, SD, pch=20, cex=3, col=adjustcolor("#CD85B9", alpha.f=0.25))
abline(fit <- lm(SD ~ means),col="#225ea8", lwd=3)
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
#plot(MCO.CV, means, pch=20, cex=3, col=adjustcolor("#CD85B9", alpha.f=0.25))
dev.off()


#paste columns in data frame to do stuff with
MCO.data= cbind(Rep1$V1, Rep1$V2, Rep1$V3, MCO.CV)
MCO.data <- as.data.frame(MCO.data)
names(MCO.data) <- c("chr","start","end","cv")
MCO.data$cv = as.numeric(as.character(MCO.data$cv))
#find top X% and subset against
top = 15
MCO.top=MCO.data[MCO.data$cv > quantile(MCO.data$cv,prob=1-top/100),]
MCO.bottom=MCO.data[MCO.data$cv < quantile(MCO.data$cv,prob=1-top/100),]
#write subsetted data to a bed file
write.table(MCO.bottom, file = "MCO.masked.25mil.localSize10kb.repTHSs.bed", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#PLOTS
MCO.h <- hist(MCO.CV, breaks = 10, density = 10, col = "lightgray", xlab = "CV", main = "MCO") 
MCO.xfit <- seq(min(MCO.CV), max(MCO.CV), length = 40) 
MCO.yfit <- dnorm(MCO.xfit, mean = mean(MCO.CV), sd = sd(MCO.CV)) 
MCO.yfit <- MCO.yfit * diff(MCO.h$mids[1:2]) * length(MCO.CV) 
lines(MCO.xfit, MCO.yfit, col = "black", lwd = 2)
abline(v=max(MCO.bottom$cv))
legend(0,3000,legend=paste("CV: ",round(max(MCO.bottom$cv), digits=3),sep=''))

####
#MZ#
####
setwd("E:/TEMP/DiffPeak_betaTesting/MZ")
#Make counts matrix of all the samples you're going to look across
Rep1=read.delim("MZ_L07.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep2=read.delim("MZ_O09.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep3=read.delim("MZ_R09.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep4=read.delim("MZ_S03.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)

#for later
Rep1[,1] = as.character(Rep1[,1])

#cbind(..., deparse.level = 1)
counts.data= cbind(Rep1$V4, Rep2$V4, Rep3$V4, Rep4$V4)
counts.data[is.na(counts.data)] <- 0
head(counts.data)

counts.data[,1] = as.numeric(as.character(counts.data[,1]))
counts.data[,2] = as.numeric(as.character(counts.data[,2]))
counts.data[,3] = as.numeric(as.character(counts.data[,3]))
counts.data[,4] = as.numeric(as.character(counts.data[,4]))

#A large number to normalize cut counts with
total = 10000000

#If data not there, mark as 0
counts.data[which(is.nan(counts.data[,1])),1] = 0
counts.data[which(is.nan(counts.data[,2])),2] = 0
counts.data[which(is.nan(counts.data[,3])),3] = 0
counts.data[which(is.nan(counts.data[,4])),4] = 0

#Normalize within sample cut counts based on the sum of cuts across all the urepTHSs
counts.data[,1] = total * counts.data[,1] / sum(counts.data[,1])
counts.data[,2] = total * counts.data[,2] / sum(counts.data[,2])
counts.data[,3] = total * counts.data[,3] / sum(counts.data[,3])
counts.data[,4] = total * counts.data[,4] / sum(counts.data[,4])

##How I determined what reps were outliers
setwd("E:/BRADY_LAB/ATLAS/ATLAS_ATAC_figures/For paper")

pdf(file='scatter_plot_replicates_repUnion_THSs_MZ_042419.pdf', height=11, width=8)
par(mfrow=c(1,6))
plot(counts.data[,1], counts.data[,2], pch=20, cex=3, col=adjustcolor("#AA1D3F", alpha.f=0.25))
plot(counts.data[,1], counts.data[,3], pch=20, cex=3, col=adjustcolor("#AA1D3F", alpha.f=0.25))
plot(counts.data[,1], counts.data[,4], pch=20, cex=3, col=adjustcolor("#AA1D3F", alpha.f=0.25))
plot(counts.data[,2], counts.data[,3], pch=20, cex=3, col=adjustcolor("#AA1D3F", alpha.f=0.25))
plot(counts.data[,2], counts.data[,4], pch=20, cex=3, col=adjustcolor("#AA1D3F", alpha.f=0.25))
plot(counts.data[,3], counts.data[,4], pch=20, cex=3, col=adjustcolor("#AA1D3F", alpha.f=0.25))
dev.off()

#needs to go across rows
means = rowMeans(counts.data)
SD = rowSds(counts.data)
MZ.CV = (SD/means)*100

pdf(file='scatter_plot_SD_mean_CV_repUnion_THSs_MZ_042419.pdf', height=11, width=8)
#par(mfrow=c(1,2))
plot(means, SD, pch=20, cex=3, col=adjustcolor("#AA1D3F", alpha.f=0.25))
abline(fit <- lm(SD ~ means),col="#225ea8", lwd=3)
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
#plot(MZ.CV, means, pch=20, cex=3, col=adjustcolor("#AA1D3F", alpha.f=0.25))
dev.off()


#paste columns in data frame to do stuff with
MZ.data= cbind(Rep1$V1, Rep1$V2, Rep1$V3, MZ.CV)
MZ.data <- as.data.frame(MZ.data)
names(MZ.data) <- c("chr","start","end","cv")
MZ.data$cv = as.numeric(as.character(MZ.data$cv))
#find top X% and subset against
top = 15
MZ.top=MZ.data[MZ.data$cv > quantile(MZ.data$cv,prob=1-top/100),]
MZ.bottom=MZ.data[MZ.data$cv < quantile(MZ.data$cv,prob=1-top/100),]
#write subsetted data to a bed file
write.table(MZ.bottom, file = "MZ.masked.25mil.localSize10kb.repTHSs.bed", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#PLOTS
MZ.h <- hist(MZ.CV, breaks = 10, density = 10, col = "lightgray", xlab = "CV", main = "MZ") 
MZ.xfit <- seq(min(MZ.CV), max(MZ.CV), length = 40) 
MZ.yfit <- dnorm(MZ.xfit, mean = mean(MZ.CV), sd = sd(MZ.CV)) 
MZ.yfit <- MZ.yfit * diff(MZ.h$mids[1:2]) * length(MZ.CV) 
lines(MZ.xfit, MZ.yfit, col = "black", lwd = 2)
abline(v=max(MZ.bottom$cv))
legend(0,5000,legend=paste("CV: ",round(max(MZ.bottom$cv), digits=3),sep=''))

####
#PH#
####
setwd("E:/TEMP/DiffPeak_betaTesting/PH")
#Make counts matrix of all the samples you're going to look across
Rep1=read.delim("PH_I06.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep2=read.delim("PH_L06.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep3=read.delim("PH_M07.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep4=read.delim("PH_Q05.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)

#for later
Rep1[,1] = as.character(Rep1[,1])

#cbind(..., deparse.level = 1)
counts.data= cbind(Rep1$V4, Rep2$V4, Rep3$V4, Rep4$V4)
counts.data[is.na(counts.data)] <- 0
head(counts.data)

counts.data[,1] = as.numeric(as.character(counts.data[,1]))
counts.data[,2] = as.numeric(as.character(counts.data[,2]))
counts.data[,3] = as.numeric(as.character(counts.data[,3]))
counts.data[,4] = as.numeric(as.character(counts.data[,4]))

#A large number to normalize cut counts with
total = 10000000

#If data not there, mark as 0
counts.data[which(is.nan(counts.data[,1])),1] = 0
counts.data[which(is.nan(counts.data[,2])),2] = 0
counts.data[which(is.nan(counts.data[,3])),3] = 0
counts.data[which(is.nan(counts.data[,4])),4] = 0

#Normalize within sample cut counts based on the sum of cuts across all the urepTHSs
counts.data[,1] = total * counts.data[,1] / sum(counts.data[,1])
counts.data[,2] = total * counts.data[,2] / sum(counts.data[,2])
counts.data[,3] = total * counts.data[,3] / sum(counts.data[,3])
counts.data[,4] = total * counts.data[,4] / sum(counts.data[,4])

##How I determined what reps were outliers
setwd("E:/BRADY_LAB/ATLAS/ATLAS_ATAC_figures/For paper")


pdf(file='scatter_plot_replicates_repUnion_THSs_PH_042419.pdf', height=22, width=8)
par(mfrow=c(1,6))
plot(counts.data[,1], counts.data[,2], pch=20, cex=3, col=adjustcolor("#ABD379", alpha.f=0.25))
plot(counts.data[,1], counts.data[,3], pch=20, cex=3, col=adjustcolor("#ABD379", alpha.f=0.25))
plot(counts.data[,1], counts.data[,4], pch=20, cex=3, col=adjustcolor("#ABD379", alpha.f=0.25))
plot(counts.data[,2], counts.data[,3], pch=20, cex=3, col=adjustcolor("#ABD379", alpha.f=0.25))
plot(counts.data[,2], counts.data[,4], pch=20, cex=3, col=adjustcolor("#ABD379", alpha.f=0.25))
plot(counts.data[,3], counts.data[,4], pch=20, cex=3, col=adjustcolor("#ABD379", alpha.f=0.25))
dev.off()

#needs to go across rows
means = rowMeans(counts.data)
SD = rowSds(counts.data)
PH.CV = (SD/means)*100

pdf(file='scatter_plot_SD_mean_CV_repUnion_THSs_PH_042419.pdf', height=11, width=8)
#par(mfrow=c(1,2))
plot(means, SD, pch=20, cex=3, col=adjustcolor("#ABD379", alpha.f=0.25))
abline(fit <- lm(SD ~ means),col="#ABD379", lwd=3)
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
#plot(PH.CV, means, pch=20, cex=3, col=adjustcolor("#ABD379", alpha.f=0.25))
dev.off()


#paste columns in data frame to do stuff with
PH.data= cbind(Rep1$V1, Rep1$V2, Rep1$V3, PH.CV)
PH.data <- as.data.frame(PH.data)
names(PH.data) <- c("chr","start","end","cv")
PH.data$cv = as.numeric(as.character(PH.data$cv))
#find top X% and subset against
top = 15
PH.top=PH.data[PH.data$cv > quantile(PH.data$cv,prob=1-top/100),]
PH.bottom=PH.data[PH.data$cv < quantile(PH.data$cv,prob=1-top/100),]
#write subsetted data to a bed file
write.table(PH.bottom, file = "PH.masked.25mil.localSize10kb.repTHSs.bed", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#PLOTS
PH.h <- hist(PH.CV, breaks = 10, density = 10, col = "lightgray", xlab = "CV", main = "PH") 
PH.xfit <- seq(min(PH.CV), max(PH.CV), length = 40) 
PH.yfit <- dnorm(PH.xfit, mean = mean(PH.CV), sd = sd(PH.CV)) 
PH.yfit <- PH.yfit * diff(PH.h$mids[1:2]) * length(PH.CV) 
lines(PH.xfit, PH.yfit, col = "black", lwd = 2)
abline(v=max(PH.bottom$cv))
legend(0,4000,legend=paste("CV: ",round(max(PH.bottom$cv), digits=3),sep=''))

###
#V#
###
#REP1/I05 is super weird was thus excluded
setwd("E:/TEMP/DiffPeak_betaTesting/V")
#Make counts matrix of all the samples you're going to look across
Rep1=read.delim("V_I05.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep2=read.delim("V_J07.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep3=read.delim("V_M08.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep4=read.delim("V_R08.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)

#for later
Rep1[,1] = as.character(Rep1[,1])

#cbind(..., deparse.level = 1)
#counts.data= cbind(Rep1$V4, Rep2$V4, Rep3$V4, Rep4$V4)
counts.data= cbind(Rep2$V4, Rep3$V4, Rep4$V4)
counts.data[is.na(counts.data)] <- 0
head(counts.data)

counts.data[,1] = as.numeric(as.character(counts.data[,1]))
counts.data[,2] = as.numeric(as.character(counts.data[,2]))
counts.data[,3] = as.numeric(as.character(counts.data[,3]))
#counts.data[,4] = as.numeric(as.character(counts.data[,4]))

#A large number to normalize cut counts with
total = 10000000

#If data not there, mark as 0
counts.data[which(is.nan(counts.data[,1])),1] = 0
counts.data[which(is.nan(counts.data[,2])),2] = 0
counts.data[which(is.nan(counts.data[,3])),3] = 0
#counts.data[which(is.nan(counts.data[,4])),4] = 0


#Normalize within sample cut counts based on the sum of cuts across all the urepTHSs
counts.data[,1] = total * counts.data[,1] / sum(counts.data[,1])
counts.data[,2] = total * counts.data[,2] / sum(counts.data[,2])
counts.data[,3] = total * counts.data[,3] / sum(counts.data[,3])
#counts.data[,4] = total * counts.data[,4] / sum(counts.data[,4])

##How I determined what reps were outliers
setwd("E:/BRADY_LAB/ATLAS/ATLAS_ATAC_figures/For paper")

pdf(file='scatter_plot_replicates_repUnion_THSs_V_042419.pdf', height=11, width=8)
par(mfrow=c(1,6))
plot(counts.data[,1], counts.data[,2], pch=20, cex=3, col=adjustcolor("#19B35A", alpha.f=0.25))
plot(counts.data[,1], counts.data[,3], pch=20, cex=3, col=adjustcolor("#19B35A", alpha.f=0.25))
plot(counts.data[,1], counts.data[,4], pch=20, cex=3, col=adjustcolor("#19B35A", alpha.f=0.25))
plot(counts.data[,2], counts.data[,3], pch=20, cex=3, col=adjustcolor("#19B35A", alpha.f=0.25))
plot(counts.data[,2], counts.data[,4], pch=20, cex=3, col=adjustcolor("#19B35A", alpha.f=0.25))
plot(counts.data[,3], counts.data[,4], pch=20, cex=3, col=adjustcolor("#19B35A", alpha.f=0.25))
dev.off()

#####################################################################################################################
#rep1/I05 appears to be an outlier, so I'm removing it and thus altered the previous lines of code for the V samples#
#####################################################################################################################

#needs to go across rows
means = rowMeans(counts.data)
SD = rowSds(counts.data)
V.CV = (SD/means)*100

pdf(file='scatter_plot_SD_mean_CV_repUnion_THSs_V_042419.pdf', height=11, width=8)
#par(mfrow=c(1,2))
plot(means, SD, pch=20, cex=3, col=adjustcolor("#19B35A", alpha.f=0.25))
abline(fit <- lm(SD ~ means),col="#225ea8", lwd=3)
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
#plot(V.CV, means, pch=20, cex=3, col=adjustcolor("#19B35A", alpha.f=0.25))
dev.off()

V.data= cbind(Rep1$V1, Rep1$V2, Rep1$V3, V.CV)
V.data <- as.data.frame(V.data)
names(V.data) <- c("chr","start","end","cv")
V.data$cv = as.numeric(as.character(V.data$cv))
#find top X% and subset against
top = 15
V.top=V.data[V.data$cv > quantile(V.data$cv,prob=1-top/100),]
V.bottom=V.data[V.data$cv < quantile(V.data$cv,prob=1-top/100),]
#write subsetted data to a bed file
write.table(V.bottom, file = "V.masked.25mil.localSize10kb.repTHSs.bed", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#PLOTS

V.h <- hist(V.CV, breaks = 10, density = 10, col = "lightgray", xlab = "CV", main = "V") 
V.xfit <- seq(min(V.CV), max(V.CV), length = 40) 
V.yfit <- dnorm(V.xfit, mean = mean(V.CV), sd = sd(V.CV)) 
V.yfit <- V.yfit * diff(V.h$mids[1:2]) * length(V.CV) 
lines(V.xfit, V.yfit, col = "black", lwd = 2)
abline(v=max(V.bottom$cv))
legend(0,2000,legend=paste("CV: ",round(max(V.bottom$cv), digits=3),sep=''))


#####
#WOX#
#####
setwd("E:/TEMP/DiffPeak_betaTesting/WOX")
#Make counts matrix of all the samples you're going to look across
Rep1=read.delim("WOX_K07.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep2=read.delim("WOX_O08.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep3=read.delim("WOX_P06.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep4=read.delim("WOX_S08.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)

#for later
Rep1[,1] = as.character(Rep1[,1])

#cbind(..., deparse.level = 1)
counts.data= cbind(Rep1$V4, Rep2$V4, Rep3$V4, Rep4$V4)
counts.data[is.na(counts.data)] <- 0
head(counts.data)

counts.data[,1] = as.numeric(as.character(counts.data[,1]))
counts.data[,2] = as.numeric(as.character(counts.data[,2]))
counts.data[,3] = as.numeric(as.character(counts.data[,3]))
counts.data[,4] = as.numeric(as.character(counts.data[,4]))

#A large number to normalize cut counts with
total = 10000000

#If data not there, mark as 0
counts.data[which(is.nan(counts.data[,1])),1] = 0
counts.data[which(is.nan(counts.data[,2])),2] = 0
counts.data[which(is.nan(counts.data[,3])),3] = 0
counts.data[which(is.nan(counts.data[,4])),4] = 0

#Normalize within sample cut counts based on the sum of cuts across all the urepTHSs
counts.data[,1] = total * counts.data[,1] / sum(counts.data[,1])
counts.data[,2] = total * counts.data[,2] / sum(counts.data[,2])
counts.data[,3] = total * counts.data[,3] / sum(counts.data[,3])
counts.data[,4] = total * counts.data[,4] / sum(counts.data[,4])

##How I determined what reps were outliers
setwd("E:/BRADY_LAB/ATLAS/ATLAS_ATAC_figures/For paper")

pdf(file='scatter_plot_replicates_repUnion_THSs_WOX_042419.pdf', height=11, width=8)
par(mfrow=c(1,6))
plot(counts.data[,1], counts.data[,2], pch=20, cex=3, col=adjustcolor("#F47752", alpha.f=0.25))
plot(counts.data[,1], counts.data[,3], pch=20, cex=3, col=adjustcolor("#F47752", alpha.f=0.25))
plot(counts.data[,1], counts.data[,4], pch=20, cex=3, col=adjustcolor("#F47752", alpha.f=0.25))
plot(counts.data[,2], counts.data[,3], pch=20, cex=3, col=adjustcolor("#F47752", alpha.f=0.25))
plot(counts.data[,2], counts.data[,4], pch=20, cex=3, col=adjustcolor("#F47752", alpha.f=0.25))
plot(counts.data[,3], counts.data[,4], pch=20, cex=3, col=adjustcolor("#F47752", alpha.f=0.25))
dev.off()


#needs to go across rows
means = rowMeans(counts.data)
SD = rowSds(counts.data)
WOX.CV = (SD/means)*100

pdf(file='scatter_plot_SD_mean_CV_repUnion_THSs_WOX_042419.pdf', height=11, width=8)
#par(mfrow=c(1,2))
plot(means, SD, pch=20, cex=3, col=adjustcolor("#F47752", alpha.f=0.25))
abline(fit <- lm(SD ~ means),col="#225ea8", lwd=3)
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
#plot(WOX.CV, means, pch=20, cex=3, col=adjustcolor("#F47752", alpha.f=0.25))
dev.off()

#paste columns in data frame to do stuff with
WOX.data= cbind(Rep1$V1, Rep1$V2, Rep1$V3, WOX.CV)
WOX.data <- as.data.frame(WOX.data)
names(WOX.data) <- c("chr","start","end","cv")
WOX.data$cv = as.numeric(as.character(WOX.data$cv))
#find top X% and subset against
top = 15
WOX.top=WOX.data[WOX.data$cv > quantile(WOX.data$cv,prob=1-top/100),]
WOX.bottom=WOX.data[WOX.data$cv < quantile(WOX.data$cv,prob=1-top/100),]
#write subsetted data to a bed file
write.table(WOX.bottom, file = "WOX.masked.25mil.localSize10kb.repTHSs.bed", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#PLOTS
WOX.h <- hist(WOX.CV, breaks = 10, density = 10, col = "lightgray", xlab = "CV", main = "WOX") 
WOX.xfit <- seq(min(WOX.CV), max(WOX.CV), length = 40) 
WOX.yfit <- dnorm(WOX.xfit, mean = mean(WOX.CV), sd = sd(WOX.CV)) 
WOX.yfit <- WOX.yfit * diff(WOX.h$mids[1:2]) * length(WOX.CV) 
lines(WOX.xfit, WOX.yfit, col = "black", lwd = 2)
abline(v=max(WOX.bottom$cv))
legend(0,4000,legend=paste("CV: ",round(max(WOX.bottom$cv), digits=3),sep=''))


####
#XY#
####
setwd("E:/TEMP/DiffPeak_betaTesting/XY")
#Make counts matrix of all the samples you're going to look across
Rep1=read.delim("XY_I07.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep2=read.delim("XY_K08.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep3=read.delim("XY_N05.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
#Rep4=read.delim("XY_Q06.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)

#for later
Rep1[,1] = as.character(Rep1[,1])

#cbind(..., deparse.level = 1)
counts.data= cbind(Rep1$V4, Rep2$V4, Rep3$V4)
counts.data[is.na(counts.data)] <- 0
head(counts.data)

counts.data[,1] = as.numeric(as.character(counts.data[,1]))
counts.data[,2] = as.numeric(as.character(counts.data[,2]))
counts.data[,3] = as.numeric(as.character(counts.data[,3]))
#counts.data[,4] = as.numeric(as.character(counts.data[,4]))

#A large number to normalize cut counts with
total = 10000000

#If data not there, mark as 0
counts.data[which(is.nan(counts.data[,1])),1] = 0
counts.data[which(is.nan(counts.data[,2])),2] = 0
counts.data[which(is.nan(counts.data[,3])),3] = 0
#counts.data[which(is.nan(counts.data[,4])),4] = 0

#Normalize within sample cut counts based on the sum of cuts across all the urepTHSs
counts.data[,1] = total * counts.data[,1] / sum(counts.data[,1])
counts.data[,2] = total * counts.data[,2] / sum(counts.data[,2])
counts.data[,3] = total * counts.data[,3] / sum(counts.data[,3])
#counts.data[,4] = total * counts.data[,4] / sum(counts.data[,4])

##How I determined what reps were outliers
setwd("E:/BRADY_LAB/ATLAS/ATLAS_ATAC_figures/For paper")

pdf(file='scatter_plot_replicates_repUnion_THSs_XY_042419.pdf', height=11, width=8)
par(mfrow=c(1,3))
plot(counts.data[,1], counts.data[,2], pch=20, cex=3, col=adjustcolor("#EDE83B", alpha.f=0.25))
plot(counts.data[,1], counts.data[,3], pch=20, cex=3, col=adjustcolor("#EDE83B", alpha.f=0.25))
plot(counts.data[,2], counts.data[,3], pch=20, cex=3, col=adjustcolor("#EDE83B", alpha.f=0.25))
dev.off()
##Rep3 is a little weird and I thus taken it out

#needs to go across rows
means = rowMeans(counts.data)
SD = rowSds(counts.data)
XY.CV = (SD/means)*100

pdf(file='scatter_plot_SD_mean_CV_repUnion_THSs_XY_042419.pdf', height=11, width=8)
#par(mfrow=c(1,2))
plot(means, SD, pch=20, cex=3, col=adjustcolor("#EDE83B", alpha.f=0.25))
abline(fit <- lm(SD ~ means),col="#225ea8", lwd=3)
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
#plot(XY.CV, means, pch=20, cex=3, col=adjustcolor("#EDE83B", alpha.f=0.25))
dev.off()


#paste columns in data frame to do stuff with
XY.data= cbind(Rep1$V1, Rep1$V2, Rep1$V3, XY.CV)
XY.data <- as.data.frame(XY.data)
names(XY.data) <- c("chr","start","end","cv")
XY.data$cv = as.numeric(as.character(XY.data$cv))
#find top X% and subset against
top = 15
XY.top=XY.data[XY.data$cv > quantile(XY.data$cv,prob=1-top/100),]
XY.bottom=XY.data[XY.data$cv < quantile(XY.data$cv,prob=1-top/100),]
#write subsetted data to a bed file
write.table(XY.bottom, file = "XY.masked.25mil.localSize10kb.repTHSs.bed", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#PLOTS
XY.h <- hist(XY.CV, breaks = 10, density = 10, col = "lightgray", xlab = "CV", main = "XY") 
XY.xfit <- seq(min(XY.CV), max(XY.CV), length = 40) 
XY.yfit <- dnorm(XY.xfit, mean = mean(XY.CV), sd = sd(XY.CV)) 
XY.yfit <- XY.yfit * diff(XY.h$mids[1:2]) * length(XY.CV) 
lines(XY.xfit, XY.yfit, col = "black", lwd = 2)
abline(v=max(XY.bottom$cv))
legend(0,2000,legend=paste("CV: ",round(max(XY.bottom$cv), digits=3),sep=''))

#back to the cluster after transferring files to their proper locations on the brady share space

##FANCY PLOT TO LOOK AT OVERLAP
### see CV_diff_peak_plots_without_inputs_files.R for actual summary plots
cd /share/brady/people/AlexMason/ATLAS_ATAC_data_analysis/BWA_mem_based_bams_and_analyses

#Use cat, bedtools sort, and  bed tools merge make union set of THSs ##X35S not included from here on out
cat COR/COR.masked.25mil.localSize10kb.repTHSs.bed EN/EN.masked.25mil.localSize10kb.repTHSs.bed EP/EP.masked.25mil.localSize10kb.repTHSs.bed EXO/EXO.masked.25mil.localSize10kb.repTHSs.bed MCO/MCO.masked.25mil.localSize10kb.repTHSs.bed MZ/MZ.masked.25mil.localSize10kb.repTHSs.bed PH/PH.masked.25mil.localSize10kb.repTHSs.bed V/V.masked.25mil.localSize10kb.repTHSs.bed WOX/WOX.masked.25mil.localSize10kb.repTHSs.bed XY/XY.masked.25mil.localSize10kb.repTHSs.bed > dummy.bed

bedtools sort -i dummy.bed > dummy.sort.bed

bedtools merge -d 10 -i dummy.sort.bed > uTHS.masked.25mil.localSize10kb.bed

rm dummy.bed
rm dummy.sort.bed


#ADD ID'S
R
my.data=read.delim("uTHS.masked.25mil.localSize10kb.bed", sep="\t", header=FALSE)
my.data$V4 <- paste("ID-",seq.int(nrow(my.data)),sep='')
write.table(my.data, file = "IDs.uTHS.masked.25mil.localSize10kb.bed", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


my.data=read.delim("uTHS.masked.25mil.localSize10kb.bed", header=FALSE, sep='\t')
names(my.data)<-c("chr", "start", "end")
my.data$length = my.data$end - my.data$start
#fancy histogram
pdf(file="histogram_unionTHS_lengths_042419.pdf", height=11, width=8)
h <- hist(my.data$length, breaks = 10, density = 10, col = "lightgray", xlab = "Length (bp)", main = "uTHSs across cell types") 
xfit <- seq(min(my.data$length), max(my.data$length), length = 40) 
yfit <- dnorm(xfit, mean = mean(my.data$length), sd = sd(my.data$length)) 
yfit <- yfit * diff(h$mids[1:2]) * length(my.data$length) 
lines(xfit, yfit, col = "black", lwd = 2)
#abline(v=max(XY.bottom$cv))
#legend(0,2000,legend=paste("CV: ",round(max(XY.bottom$cv), digits=3),sep=''))
dev.off()
q()

#now turn into GFF file
python /share/brady/people/AlexMason/Alex_scripts/BedtoGff.py IDs.uTHS.masked.25mil.localSize10kb.bed IDs.uTHS.masked.25mil.localSize10kb.cv30.gff

python /share/brady/people/AlexMason/Alex_scripts/BedtoGff.py all_unionpeaks.bed all_unionpeaks.gff
#python /net/gs/vol2/home/gmason88/Alex_scripts/Python_scripts/BedtoGff.py IDs.uTHS.masked.25mil.localSize10kb.cv15.bed IDs.uTHS.masked.25mil.localSize10kb.cv15.gff
EVIDENCE OF ISNURANCE/ BINDER
#now use excel to make sure in this format

##gff-version 3
##sequence-region   SL3.0ch00 1 20852292
##sequence-region   SL3.0ch01 1 98455869
##sequence-region   SL3.0ch02 1 55977580
##sequence-region   SL3.0ch03 1 72290146
##sequence-region   SL3.0ch04 1 66557038
##sequence-region   SL3.0ch05 1 66723567
##sequence-region   SL3.0ch06 1 49794276
##sequence-region   SL3.0ch07 1 68175699
##sequence-region   SL3.0ch08 1 65987440
##sequence-region   SL3.0ch09 1 72906345
##sequence-region   SL3.0ch10 1 65633393
##sequence-region   SL3.0ch11 1 56597135
##sequence-region   SL3.0ch12 1 68126176
#SL3.0ch01	bed2gff	uTHS	77	227	0	+	.	ID00001
#SL3.0ch01	bed2gff	uTHS	21771	21921	0	+	.	ID00002
#SL3.0ch01	bed2gff	uTHS	83219	83369	0	+	.	ID00003

#now use bedmap to counts tags in uTHSs
#COR DONE
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed COR/COR_K06.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > COR/COR_K06.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed COR/COR_O07.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > COR/COR_O07.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed COR/COR_R07.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > COR/COR_R07.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed COR/COR_S06.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > COR/COR_S06.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &

#EN DONE
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed EN/EN_L05.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > EN/EN_L05.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed EN/EN_N06.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > EN/EN_N06.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed EN/EN_R06.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > EN/EN_R06.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed EN/EN_S09.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > EN/EN_S09.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &

#EP DONE
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed EP/EP_J05.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > EP/EP_J05.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed EP/EP_N04.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > EP/EP_N04.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed EP/EP_Q07.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > EP/EP_Q07.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed EP/EP_S07.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > EP/EP_S07.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &

#EXO DONE
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed EXO/EXO_J06.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > EXO/EXO_J06.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed EXO/EXO_M05.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > EXO/EXO_M05.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed EXO/EXO_P04.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > EXO/EXO_P04.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed EXO/EXO_S05.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > EXO/EXO_S05.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &

#MCO 
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed MCO/MCO_I08.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > MCO/MCO_I08.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed MCO/MCO_K05.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > MCO/MCO_K05.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed MCO/MCO_M06.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > MCO/MCO_M06.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed MCO/MCO_P05.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > MCO/MCO_P05.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &

#MZ 
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed MZ/MZ_L07.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > MZ/MZ_L07.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed MZ/MZ_O09.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > MZ/MZ_O09.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed MZ/MZ_R09.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > MZ/MZ_R09.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed MZ/MZ_S03.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > MZ/MZ_S03.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &

#PH
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed PH/PH_I06.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > PH/PH_I06.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed PH/PH_L06.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > PH/PH_L06.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed PH/PH_M07.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > PH/PH_M07.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed PH/PH_Q05.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > PH/PH_Q05.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &

#V replicate I05 removed because it's weird
#bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed V/V_I05.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > V/V_I05.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed V/V_J07.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > V/V_J07.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed V/V_M08.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > V/V_M08.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed V/V_R08.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > V/V_R08.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &

#WOX
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed WOX/WOX_K07.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > WOX/WOX_K07.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed WOX/WOX_O08.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > WOX/WOX_O08.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed WOX/WOX_P06.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > WOX/WOX_P06.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed WOX/WOX_S08.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > WOX/WOX_S08.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &

#X35S
#bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.bed X35S/X35S_J08.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > X35S/X35S_J08.bwa.sorted.q2.noCh0.masked.25mil.counts.per.uTHS.bed &
#bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.bed X35S/X35S_O06.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > X35S/X35S_O06.bwa.sorted.q2.noCh0.masked.25mil.counts.per.uTHS.bed &
#bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.bed X35S/X35S_Q08.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > X35S/X35S_Q08.bwa.sorted.q2.noCh0.masked.25mil.counts.per.uTHS.bed &
#bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.bed X35S/X35S_S02.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > X35S/X35S_S02.bwa.sorted.q2.noCh0.masked.25mil.counts.per.uTHS.bed &

#XY
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed XY/XY_I07.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > XY/XY_I07.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed XY/XY_K08.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > XY/XY_K08.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed XY/XY_N05.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > XY/XY_N05.bwa.sorted.q2.noCh0.masked.25mil.counts.per.cv15.uTHS.bed &
#bedmap --echo --sum --delim "\t" IDs.uTHS.masked.25mil.localSize10kb.cv15.bed XY/XY_Q06.bwa.sorted.q2.noCh0.masked.25mil.perbase.bed > XY/XY_Q06.bwa.sorted.q2.noCh0.masked.25mil.counts.per.uTHS.bed &



#################################################
#X35S, not repeated with new cutoffs on 04/24/19#
#################################################
setwd("E:/TEMP/DiffPeak_betaTesting/X35S")
#Make counts matrix of all the samples you're going to look across
Rep1=read.delim("X35S_J08.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep2=read.delim("X35S_O06.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep3=read.delim("X35S_Q08.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)
Rep4=read.delim("X35S_S02.bwa.sorted.q2.noCh0.masked.25mil.repUnion.countsuTHS.bed", sep="\t", header=FALSE)

#for later
Rep1[,1] = as.character(Rep1[,1])

#cbind(..., deparse.level = 1)
counts.data= cbind(Rep1$V4, Rep2$V4, Rep3$V4, Rep4$V4)
counts.data[is.na(counts.data)] <- 0
head(counts.data)

counts.data[,1] = as.numeric(as.character(counts.data[,1]))
counts.data[,2] = as.numeric(as.character(counts.data[,2]))
counts.data[,3] = as.numeric(as.character(counts.data[,3]))
counts.data[,4] = as.numeric(as.character(counts.data[,4]))

#A large number to normalize cut counts with
total = 10000000

#If data not there, mark as 0
counts.data[which(is.nan(counts.data[,1])),1] = 0
counts.data[which(is.nan(counts.data[,2])),2] = 0
counts.data[which(is.nan(counts.data[,3])),3] = 0
counts.data[which(is.nan(counts.data[,4])),4] = 0

#Normalize within sample cut counts based on the sum of cuts across all the urepTHSs
counts.data[,1] = total * counts.data[,1] / sum(counts.data[,1])
counts.data[,2] = total * counts.data[,2] / sum(counts.data[,2])
counts.data[,3] = total * counts.data[,3] / sum(counts.data[,3])
counts.data[,4] = total * counts.data[,4] / sum(counts.data[,4])

##How I determined what reps were outliers
pdf(file='scatter_plot_replicates_repUnion_THSs.pdf', height=11, width=8)
par(mfrow=c(2,3))
plot(counts.data[,1], counts.data[,2])
plot(counts.data[,1], counts.data[,3])
plot(counts.data[,1], counts.data[,4])
plot(counts.data[,2], counts.data[,3])
plot(counts.data[,2], counts.data[,4])
plot(counts.data[,3], counts.data[,4])
dev.off()


#needs to go across rows
means = rowMeans(counts.data)
SD = rowSds(counts.data)
X35S.CV = (SD/means)*100

#paste columns in data frame to do stuff with
X35S.data= cbind(Rep1$V1, Rep1$V2, Rep1$V3, X35S.CV)
X35S.data <- as.data.frame(X35S.data)
names(X35S.data) <- c("chr","start","end","cv")
X35S.data$cv = as.numeric(as.character(X35S.data$cv))
#find top 15% and subset against
top = 15
X35S.top=X35S.data[X35S.data$cv > quantile(X35S.data$cv,prob=1-top/100),]
X35S.bottom=X35S.data[X35S.data$cv < quantile(X35S.data$cv,prob=1-top/100),]
#write subsetted data to a bed file
write.table(X35S.bottom, file = "X35S.masked.25mil.localSize10kb.repTHSs.bed", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#PLOTS
X35S.h <- hist(X35S.CV, breaks = 10, density = 10, col = "lightgray", xlab = "CV", main = "X35S") 
X35S.xfit <- seq(min(X35S.CV), max(X35S.CV), length = 40) 
X35S.yfit <- dnorm(X35S.xfit, mean = mean(X35S.CV), sd = sd(X35S.CV)) 
X35S.yfit <- X35S.yfit * diff(X35S.h$mids[1:2]) * length(X35S.CV) 
lines(X35S.xfit, X35S.yfit, col = "black", lwd = 2)
abline(v=max(X35S.bottom$cv))
legend(0,2000,legend=paste("CV: ",round(max(X35S.bottom$cv), digits=3),sep=''))