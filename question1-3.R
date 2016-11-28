rm(list=ls())
#install.packages("knitr")
library(knitr)
library(GEOquery)

######## homework 1 #########
######## question 1 #########
# orignally I ran this in the HW1 directory
# gse19439 <- getGEO("GSE19439",GSEMatrix=T)
# gse19439 <- gse19439[[1]]
# save(file="gse19439.RData",gse19439)

#so I am just going to load it now
load(file="gse19439.RData")
gse19439
experimentData(gse19439)
feat <- featureData(gse19439)
varLabels(feat)

# make a factor object for Control, latent TB and Positive TB
pData(phenoData(gse19439))[,1]

head(pData(phenoData(gse19439)))

tmp <- as.character(pData(phenoData(gse19439))[,1])
J <- length(tmp) # J=number of samples
TBgroup <- rep("",J)
for(j in 1:J) TBgroup[j]=substring(tmp[j],1,3)
# make a factor for TBgroup
FTB <- factor(TBgroup,levels=c("CON","LTB","PTB"))

# get our expression set
X <- exprs(gse19439)

# do a quick kruskal-wallis scan
myKrusk <- function(i){
        cat(i,"...",fill=F)
        kruskal.test(x=X[i,], g=FTB)$p.value
}

# originally, I ran this in the HW1 directory:
#myPvals <- mapply(myKrusk,1:(dim(X)[1]))

#save(file="myPvals.RData",myPvals)
# so that now I can save time and just load:
load("myPvals.RData")

# populate vector with last names of the groups in the class.
# note: the code is written this way, because it used to be student names and not Group names
GroupLabels <- c("Group I","Group II","Group III","Group IV")

# pick the best 4 p-values and assign them to the students.
best4 <- order(myPvals)[1:4]
# print out list of best 4
print(best4)

for(i in 1:length(GroupLabels)){
        cat("Group Label:", GroupLabels[i], "\t\t row assignment:", best4[i], fill=T)
}
myrow <- gse19439[10685,]
dim(myrow)
featureNames(myrow)
gene.info <- pData(featureData(myrow))

#10685 GENE IS LACTB
#LACTB (Lactamase, Beta) is a Protein Coding gene. 
#Diseases associated with LACTB include lung abscess and bacterial conjunctivitis. 
#GO annotations related to this gene include hydrolase activity.

######### question 2 ########
sampleNames(myrow)
boxplot(log2(X[10685,sampleNames(myrow)[FTB == 'CON']]),
        log2(X[10685,sampleNames(myrow)[FTB == 'LTB']]),
        log2(X[10685,sampleNames(myrow)[FTB == 'PTB']]),
        main = "Boxplot of Row Data as Function of TB Phenotype",
        names = c("CON", "LTB", "PTB"))
#boxplots disguise multimodal data 

#install.packages('vioplot')
library(vioplot)
vioplot(log2(X[10685,sampleNames(myrow)[FTB == 'CON']]),
        log2(X[10685,sampleNames(myrow)[FTB == 'LTB']]),
        log2(X[10685,sampleNames(myrow)[FTB == 'PTB']]),
        names = c("CON", "LTB", "PTB"))
title("Violin Plot of Row Data as Function of TB Phenotype")

#install.packages('beanplot')
library(beanplot)
beanplot(log2(X[10685,sampleNames(myrow)[FTB == 'CON']]),
        log2(X[10685,sampleNames(myrow)[FTB == 'LTB']]),
        log2(X[10685,sampleNames(myrow)[FTB == 'PTB']]),
        names = c("CON", "LTB", "PTB"))
title("Bean Plot of Row Data as Function of TB Phenotype")

##### question 3 #####
#install.packages('gplots')
library(gplots)
best20 <- order(myPvals)[1:20]
x.best <- X[best20,]
heatmap.2(X[best20,], main="20 student features", scale = "row", trace="none")
heatmap.2(X[best20,], main="20 student features", scale = "row", trace="none", labCol = FTB) #change sample names into group name
#Does a decent job discriminating the groups
#it would be cool if someone got the other heat map to work that labels the columns by color with the group names

