##### question 4 ######
#remove NAs
#mt.maxT command in mlttest package for multiple testing
#plug p-val or abs(test-statistic) in because its a single cut-off and not a double cutoff
#ROC curves are predicated on single cutoff discrimination
#pROC package for the ROC curves
#"expresso", bgcorrect, normalization, PM correction, summarization
rm(list=ls())
library(affy)
library(affyPLM)
# this provides an AffyBatch object
load("PSpikeAffyBatch.RData")
load("question4.RData")

spikeDF <- read.table(file="AffyProbeSpikeValues.csv",sep="\t")
levels(spikeDF[,2])
summary(spikeDF[,2])

#grab Spike fold changes for all entries as numeric ... so anything that was not a number is now an NA
SpikeFC <- as.numeric(levels(spikeDF[,2])[spikeDF[,2]])

#grab IDs
names(SpikeFC) <- spikeDF$V1

#remove NAs
nonZeroDX <- which(!is.na(SpikeFC) & (SpikeFC != 0)) 
spikeFC.clean <- SpikeFC[nonZeroDX]

#check how many genes are left in dataset
length(spikeFC.clean)

#the eight routes
bgcorrect.mtd <- c("rma", "mas", "rma", "mas", "rma", "none", "mas", "mas")
normalized.mtd <- c("constant", "quantiles", "quantiles", "loess", "loess", "constant", "qspline", "qspline")
pmcorrect.mtd <- c("pmonly", "pmonly", "subtractmm", "mas", "mas", "pmonly", "subtractmm", "mas")
summary.mtd <- c("mas", "mas", "avgdiff", "medianpolish", "medianpolish", "avgdiff", "mas", "mas")

route.expsets <- list()
for (i in 1:length(bgcorrect.mtd)){
        routes <- expresso(affydata, 
                           bgcorrect.method = bgcorrect.mtd[i], 
                           normalize.method = normalized.mtd[i],
                           pmcorrect.method = pmcorrect.mtd[i], 
                           summary.method = summary.mtd[i])
        route.expsets[[i]] <- exprs(routes)[nonZeroDX,]
}
length(route.expsets)

#create labels for multitesting class labels -- 18 samples, 9 control, 9 experimental
labels <- factor(c(rep(0, 9), rep(1, 9))) #0 is control, 1 is experimental

#source("http://bioconductor.org/biocLite.R")
#biocLite("multtest")				
library(multtest)
library(pROC)
stats <- list()
for (i in 1:length(route.expsets)){
        `testing.routes <- mt.maxT(route.expsets[[i]], classlabel = labels, B = 10000) #conduct multiple t test for 10000 permutations  
        stats[[i]] <- testing.routes$teststat[order(testing.routes$index)]
}

nrow <- nrow(route.expsets[[1]])
myresponse <- rep(NA, nrow)
myresponse[which(spikeFC.clean == 1)] = 1 #if the spike value is 1 ... assign as 0s in matrix ... 1s are the control	
myresponse[which(spikeFC.clean != 1)] = 0 #if spike value is not 1 ... assign as 1s in the matrix ... 0s are exp.
roc.fnct <- function(x){roc(response = myresponse, predictor = abs(x), plot=TRUE, print.auc=TRUE)} #apply roc function on input 
roc <- lapply(stats, roc.fnct) #apply roc function on all the ordered test statistics in the list stats

#time to make roc curves
pdf('homework1-question4-ROCcurves.pdf')
rainbow <-  palette(rainbow(length(route.expsets)))
plot(roc[[1]], main = "ROC curves for different normalization routes using expresso()", col=rainbow[1])
legendText <- c()
for(i in 1:length(route.expsets)){
        plot(roc[[i]], add=TRUE, col=rainbow[i])
        legendText[i] <- paste(bgcorrect.mtd[i],"/",normalized.mtd[i],"/", 
                               pmcorrect.mtd[i],"/",summary.mtd[i], "   AUC: ", round(as.numeric(roc[[i]]$auc),3), sep="")
}
legend("bottomright",
       legendText,
       lty=c(rep(1,length(route.expsets))), 				
       lwd=c(rep(3,length(route.expsets))),
       col=rainbow)
dev.off()


##### CHECKING TO SEE IF SAME AS GROUP ONE -- LOOKS LIKE THEIRS COULD BE WRONG -- THEY KEPT THE NON-NUMERICAL from spikeDF in their counts
#this is group 1s code
route1 <- expresso(affydata, bgcorrect.method = "rma", normalize.method = "constant", pmcorrect.method = "pmonly",
                   summary.method = "mas")
data2 = which(( spikeDF[, 2]	 != 0) & (!is.na( spikeDF[, 2]	)))
v =  spikeDF[, 2][data2]				
data3 = which(v == 1)				
data4 = which(v != 1)	
grp1.rt1 = exprs(route1)[data2, ]				
ln = dim(grp1.rt1)[1]
myresponse = rep(NA, ln)
myresponse[data3] = 1	
myresponse[data4] = 0
grp1.rt1.testing <- mt.maxT(grp1.rt1, classlabel = labels, B = 10000) #conduct multiple t test for 10000 permutations
grp1.stats.rt1 <- grp1.rt1.testing$teststat[order(grp1.rt1.testing$index)]
grp1.roc1 = roc(response = myresponse, predictor = abs(grp1.stats.rt1), plot=TRUE, print.auc=TRUE)
#Data: abs(grp1.stats.rt1) in 3426 controls (myresponse 0) < 2189 cases (myresponse 1).
#Area under the curve: 0.8996

#plot both now
#pdf("group1vgroup2_sameroute.pdf")
plot(roc[[1]], main = "ROC for route1", xlab="1 - Specificity", col="BLUE")
plot(grp1.roc1, add=TRUE, col='GREEN' )
legend("bottomright", c(paste("Group 2 - rma/constant/pmonly/mas, AUC: ", round(as.numeric(roc[[1]]$auc),4)), 
                        paste("Group 1 - rma/constant/pmonly/mas, AUC:", sep=" ", round(as.numeric(grp1.roc1$auc),4))), 
       lty=c(1,1), lwd=c(1,1),col=c("BLUE","GREEN")) 
#dev.off()
#

save.image("question4.RData")
