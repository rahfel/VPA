
library(UpSetR)
library(ggplot2)
library(reshape2)
library(dplyr)
library(WebGestaltR)


rm(list=ls())
tmp = list.files(pattern="*.txt")
Allfiles = lapply(tmp,read.delim,sep="\t")
filenames=list.files(pattern="*.txt", full.names=TRUE)
filenames
names(Allfiles) = gsub(".*/(.*)\\..*", "\\1", filenames)
names(Allfiles)
Allfiles2 <- Allfiles[ grepl("*Regulated$", names(Allfiles)) ]
names(Allfiles2)


Inpath="/Users/rafel137/Documents/"

for (i in 1:length(Allfiles2)){
WebGestaltR(enrichMethod="ORA", organism="rnorvegicus",
enrichDatabase="geneontology_Biological_Process_noRedundant", interestGeneFile=paste0(Inpath,names(Allfiles2)[i],".txt"),
interestGeneType="ensembl_gene_id", referenceGeneFile=paste0(Inpath,names(Allfiles2)[i],".AllGenes.txt"),
referenceGeneType="ensembl_gene_id", isOutput=TRUE,
outputDirectory=Inpath, projectName=paste0("BP2.",names(Allfiles2)[i]))}
