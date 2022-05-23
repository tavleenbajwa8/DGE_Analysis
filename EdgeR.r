#Differential Gene Analysis (EdgeR)

# Install the following libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library(edgeR)
library(gplots) 
library(dplyr)
library(tidyverse)

# set the current working directory
#setwd("")  #Where the raw_counts.csv file is present

#input file
countdata<- read.csv("raw_counts_HCT116.csv", sep = ",", header = T, row.names = 1)
head(countdata)
nrow((countdata))

#levels
group <- c(rep("Untreated",2),rep("Treated",2))
head(group)

#DGElist creation
y<- DGEList(counts = countdata, group = group)

#Filtering
head(y)
myCPM<- cpm(y) #CPM calculation
head(myCPM)
head(group)
design<- model.matrix(~0+group)
#head(design)
keep<- filterByExpr(y, design = design)
#head(keep)
y<- y[keep, ,keep.lib.sizes= FALSE]
head(y)
nrow(y)
#annotation using biomaRt
library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
list<-c(row.names(y), quote(""))
head(list)
length(list)
write.table(list, file="ENSIDS.csv",quote = FALSE, sep = ",")
extraInfo <- getBM(attributes = c("ensembl_gene_id","gene_biotype","external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = list,
                   mart = ensembl, useCache = FALSE)
head(extraInfo)
nrow(extraInfo)
annotatedcounts <- merge(y$counts, extraInfo, by.x = 0, by.y = "ensembl_gene_id", sort=FALSE)
head(annotatedcounts)
nrow(annotatedcounts)
# only choosing protein coding genes
finalanno<-annotatedcounts[grep("protein_coding", annotatedcounts$gene_biotype),] #only choosing protein coding genes
nrow(finalanno)
head(finalanno)
# Remove duplicates
finalanno2<-finalanno[!duplicated(finalanno$external_gene_name),]
nrow(finalanno2)
head(finalanno2)
row.names(finalanno2)<- finalanno2$external_gene_name
head(finalanno2)
finalanno2<-finalanno2[,2:5]
head(finalanno2)
nrow(finalanno2)

#normalization

dgeObj <- DGEList(counts = finalanno2, group = group)
head(dgeObj)
dgeObj$samples$group
dgeObj <- calcNormFactors(dgeObj, method = "upperquartile") #normalizing the RNA composition using TMM
head(dgeObj)

summary(dgeObj)
head(dgeObj$samples)
dgeObj <- estimateCommonDisp(dgeObj, verbose=TRUE)
dgeObj <- estimateGLMTrendedDisp(dgeObj, verbose=TRUE)
dgeObj <- estimateTagwiseDisp(dgeObj,verbose=TRUE)

head(dgeObj)
normlogcounts <- cpm(dgeObj,log=TRUE)
head(normlogcounts)
write.table(normlogcounts, file = "normcount.csv",
            quote = FALSE, sep = ",", col.names = NA)


#differential analysis
head(dgeObj)
dgeanalysistest <- exactTest(dgeObj,pair = c("Untreated","Treated"))#taking WT as control. now +ve logFC means that it is upregulated in null
head(dgeanalysistest)
nrow(dgeanalysistest)
dgeanalysistest$table$PValue_fdr <- p.adjust(method="fdr",p=dgeanalysistest$table$PValue)
write.table(dgeanalysistest, file = "Differentiallyexpressed.csv",
            quote = FALSE, sep = ",")

#summary(dgeanalysistest)
resul<-decideTestsDGE(dgeanalysistest, adjust.method= "fdr", p.value = 0.01,lfc = 2)
nrow(resul)
summary(resul)
