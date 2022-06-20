################### Dampened weighted least squares (DWLS) method for gene expression deconvolution ############

data <- read.csv("/Users/rafel137/Documents/Aus_Samples/DECON/GSE124952_expression_matrix.csv",header=T,row.names=1)
dim(data) # observe number of rows(genes) and columns(cells)
meta <- read.csv("/Users/rafel137/Documents/Aus_Samples/DECON/GSE124952_meta_data.csv",header=T,row.names=1)
dim(meta)
colnames(meta)

scdata=data

pval.cutoff=0.01
diff.cutoff=0.5
path="/Users/rafel137/Documents/Aus_Samples/DECON2/DWLS.results"
id=c("Astro","Endo","Excitatory","Microglia","Oligo","OPC","Inhibitory")
#for each cell type, choose genes in which FDR adjusted p-value is less than 0.01 and the estimated fold-change
#is greater than 0.5
numberofGenes<-c()
for (i in unique(id)){
    if(file.exists(paste(path,"/",i,"_MIST.RData", sep=""))){
        load(file=paste(path,"/",i,"_MIST.RData", sep=""))
        pvalue_adjusted<-p.adjust(cluster_lrTest.table[,3], method = "fdr", n = length(cluster_lrTest.table[,3]))
        cluster_lrTest.table<-cbind(cluster_lrTest.table,pvalue_adjusted)
        DEGenes<-cluster_lrTest.table$Gene[intersect(which(pvalue_adjusted<pval.cutoff),which(cluster_lrTest.table$log2fold_change>diff.cutoff))]
        nonMir = grep("MIR|Mir", DEGenes, invert = T)  # because Mir gene is usually not accurate
        assign(paste("cluster_lrTest.table.",i,sep=""),cluster_lrTest.table[which(cluster_lrTest.table$Gene%in%DEGenes[nonMir]),])
        numberofGenes<-c(numberofGenes,length(DEGenes[nonMir]))
    }
}
conditionNumbers<-c()
for(G in 50:200){
    Genes<-c()
    j=1
    for (i in unique(id)){
        if(numberofGenes[j]>0){
            temp<-paste("cluster_lrTest.table.",i,sep="")
            temp<-as.name(temp)
            temp<-eval(parse(text = temp))
            temp<-temp[order(temp$log2fold_change,decreasing=TRUE),]
            Genes<-c(Genes,varhandle::unfactor(temp$Gene[1:min(G,numberofGenes[j])]))
        }
        j=j+1
    }
    Genes<-unique(Genes)
    #make signature matrix
    ExprSubset<-scdata[Genes,]
    Sig<-NULL
    for (i in unique(id)){
        Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
    }
    colnames(Sig)<-unique(id)
    conditionNumbers<-c(conditionNumbers,kappa(Sig))
}
G<-which.min(conditionNumbers)+min(49,numberofGenes-1)
Genes<-c()
j=1
for (i in unique(id)){
    if(numberofGenes[j]>0){
        temp<-paste("cluster_lrTest.table.",i,sep="")
        temp<-as.name(temp)
        temp<-eval(parse(text = temp))
        temp<-temp[order(temp$log2fold_change,decreasing=TRUE),]
        Genes<-c(Genes,varhandle::unfactor(temp$Gene[1:min(G,numberofGenes[j])]))
    }
    j=j+1
}
Genes<-unique(Genes)
ExprSubset<-scdata[Genes,]
Sig<-NULL
for (i in unique(id)){
    Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
}
colnames(Sig)<-unique(id)
save(Sig,file=paste(path,"/Sig.RData",sep=""))
return(Sig)
}

Signature=Sig



setwd("/Users/rafel137/Documents/AusSamples")
library(e1071)
library(quadprog)



#solve using OLS, constrained such that cell type numbers>0
solveOLS<-function(S,B){
    D<-t(S)%*%S
    d<-t(S)%*%B
    A<-cbind(diag(dim(S)[2]))
    bzero<-c(rep(0,dim(S)[2]))
    solution<-solve.QP(D,d,A,bzero)$solution
    names(solution)<-colnames(S)
    print(round(solution/sum(solution),5))
    return(solution/sum(solution))
}

#return cell number, not proportion
#do not print output
solveOLSInternal<-function(S,B){
    D<-t(S)%*%S
    d<-t(S)%*%B
    A<-cbind(diag(dim(S)[2]))
    bzero<-c(rep(0,dim(S)[2]))
    solution<-solve.QP(D,d,A,bzero)$solution
    names(solution)<-colnames(S)
    return(solution)
}

#solve using WLS with weights dampened by a certain dampening constant
solveDampenedWLS<-function(S,B){
    #first solve OLS, use this solution to find a starting point for the weights
    solution<-solveOLSInternal(S,B)
    #now use dampened WLS, iterate weights until convergence
    iterations<-0
    changes<-c()
    #find dampening constant for weights using cross-validation
    j<-findDampeningConstant(S,B,solution)
    change<-1
    while(change>.01 & iterations<1000){
        newsolution<-solveDampenedWLSj(S,B,solution,j)
        #decrease step size for convergence
        solutionAverage<-rowMeans(cbind(newsolution,matrix(solution,nrow = length(solution),ncol = 4)))
        change<-norm(as.matrix(solutionAverage-solution))
        solution<-solutionAverage
        iterations<-iterations+1
        changes<-c(changes,change)
    }
    print(round(solution/sum(solution),5))
    return(solution/sum(solution))
}

#solve WLS given a dampening constant
solveDampenedWLSj<-function(S,B,goldStandard,j){
    multiplier<-1*2^(j-1)
    sol<-goldStandard
    ws<-as.vector((1/(S%*%sol))^2)
    wsScaled<-ws/min(ws)
    wsDampened<-wsScaled
    wsDampened[which(wsScaled>multiplier)]<-multiplier
    W<-diag(wsDampened)
    D<-t(S)%*%W%*%S
    d<- t(S)%*%W%*%B
    A<-cbind(diag(dim(S)[2]))
    bzero<-c(rep(0,dim(S)[2]))
    sc <- norm(D,"2")
    solution<-solve.QP(D/sc,d/sc,A,bzero)$solution
    names(solution)<-colnames(S)
    return(solution)
}

#find a dampening constant for the weights using cross-validation
findDampeningConstant<-function(S,B,goldStandard){
    solutionsSd<-NULL
    #goldStandard is used to define the weights
    sol<-goldStandard
    ws<-as.vector((1/(S%*%sol))^2)
    wsScaled<-ws/min(ws)
    wsScaledMinusInf<-wsScaled
    #ignore infinite weights
    if(max(wsScaled)=="Inf"){
        wsScaledMinusInf<-wsScaled[-which(wsScaled=="Inf")]
    }
    #try multiple values of the dampening constant (multiplier)
    #for each, calculate the variance of the dampened weighted solution for a subset of genes
    for (j in 1:ceiling(log2(max(wsScaledMinusInf)))){
        multiplier<-1*2^(j-1)
        wsDampened<-wsScaled
        wsDampened[which(wsScaled>multiplier)]<-multiplier
        solutions<-NULL
        seeds<-c(1:100)
        for (i in 1:100){
            set.seed(seeds[i]) #make nondeterministic
            subset<-sample(length(ws),size=length(ws)*0.5) #randomly select half of gene set
            #solve dampened weighted least squares for subset
            fit = lm (B[subset] ~ -1+S[subset,],weights=wsDampened[subset])
            sol<-fit$coef*sum(goldStandard)/sum(fit$coef)
            solutions<-cbind(solutions,sol)
        }
        solutionsSd<-cbind(solutionsSd,apply(solutions,1,sd))
    }
    #choose dampening constant that results in least cross-validation variance
    j<-which.min(colMeans(solutionsSd^2))
    return(j)
}

solveSVR<-function(S,B){
    #scaling
    ub=max(c(as.vector(S),B)) #upper bound
    lb=min(c(as.vector(S),B)) #lower bound
    Bs=((B-lb)/ub)*2-1
    Ss=((S-lb)/ub)*2-1

    #perform SVR
    model<-svm(Ss,Bs, nu=0.5,scale = TRUE, type = "nu-regression",kernel ="linear",cost = 1)
    coef <- t(model$coefs) %*% model$SV
    coef[which(coef<0)]<-0
    coef<-as.vector(coef)
    names(coef)<-colnames(S)
    print(round(coef/sum(coef),5))
    return(coef/sum(coef))
}

#perform DE analysis using Seurat
DEAnalysis<-function(scdata,id,path){
    exprObj<-CreateSeuratObject(raw.data=as.data.frame(scdata), project = "DE")
    exprObj2<-SetIdent(exprObj,ident.use=as.vector(id))
    print("Calculating differentially expressed genes:")
    for (i in unique(id)){
        de_group <- FindMarkers(object=exprObj2, ident.1 = i, ident.2 = NULL,
                                only.pos = TRUE, test.use = "bimod")
        save(de_group,file=paste(path,"/de_",i,".RData",sep=""))
    }
}

load("/Users/rafel137/Documents/Aus_Samples/DECON2/DWLS.results/Sig.RData")

Signature=Sig

COUNTS <- read.csv("/Users/rafel137/Documents/AusSamples/mat_merged_matrix.E.GAERS.csv",header=T,row.names=1)
COUNTS[1:2,1:2]
colnames(COUNTS)
COUNTS=COUNTS[,-4]
ensembl = useMart("ensembl", dataset="rnorvegicus_gene_ensembl")
values <- unique(rownames(COUNTS))

data <- getBM(

    attributes=c(
        "ensembl_gene_id",
        "mmusculus_homolog_ensembl_gene",
        "mmusculus_homolog_associated_gene_name"),

    filters = "ensembl_gene_id",
    values = values,
    mart= ensembl
)

data <- data[duplicated(data$ensembl_gene_id)==F &
                 duplicated(data$mmusculus_homolog_associated_gene_name)==F ,]

rownames(data) <- data$ensembl_gene_id

data <- data[!data$mmusculus_homolog_associated_gene_name=="",]

mouse.COUNTS <- COUNTS[rownames(data),]

CPM=mouse.COUNTS

# 2 step process of making your mouse gene symbols (to avoding a formatting error)
Name <- data$mmusculus_homolog_associated_gene_name
rownames(CPM) <- Name
dim(CPM)

CPM.trim <- CPM[rownames(Signature),] # filters your counts for genes maintained in the Signature Matrix

samples<- colnames(CPM.trim) # grabbing a vector of your subject names

#trim bulk and single-cell data to contain the same genes
trimData<-function(Signature,bulkData){
    Genes<-intersect(rownames(Signature),names(bulkData))
    B<-bulkData[Genes]
    S<-Signature[Genes,]
    return(list("sig"=S,"bulk"=B))
}

# Initializing Results Data Frames for Each Algorithm that DWLS employs (DWLS being the superstar)
DWLS <- data.frame()
SVR <- data.frame()
OLS <- data.frame()

# Computing Cell Fractions for Each Sample and Adding to Results Data Frames
for(sample in samples){
    bulk <- CPM.trim[,sample]
    names(bulk) <- rownames(CPM.trim)
    tr<-trimData(Signature,bulk)
    solDWLS<-solveDampenedWLS(tr$sig,tr$bulk)
    solSVR <- solveSVR(tr$sig,tr$bulk)
    solOLS <- solveOLS(tr$sig,tr$bulk)
    DWLS <-rbind(DWLS,solDWLS)
    SVR <- rbind(SVR,solSVR)
    OLS <- rbind(OLS,solOLS)
}

model<- list(DWLS,SVR,OLS) # List of dataframes

# Adding the rownames and columns back to each dataframe in list
L <- lapply(model, function(df)
{
    colnames(df) <- colnames(Signature)
    rownames(df) <- samples
    df
}
)



write.csv(L[1],file="CellFractions_PFCref_DWLS.E.GAERS.csv")
write.csv(L[2],file="CellFractions_PFCref_SVR.E.GAERS.csv")
write.csv(L[3],file="CellFractions_PFCref_OLS.E.GAERS.csv")


L[1][[1]]

library(ggplot2)
library(reshape2)
N.NEC.DWLS=melt(L[1][[1]])
N.NEC.DWLS

colnames(N.NEC.DWLS)=c("Celltypes","Estimated.Celltype.Fraction")

DWLS Cell Fraction Estimates

Estimated Celltype Fraction

Celltypes



setwd("/Users/rafel137/Documents/Aus_Samples/DECON_MOUSE/DWLS.results")




rm(list=ls())
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)


tmp = list.files(pattern="CellFractions_PFCref_DWLS.*MOUSEEMB18.csv")
names(tmp)=tmp
Allfiles = lapply(tmp,read.csv,sep=",",row.names=1)
names(Allfiles)
names(Allfiles)=gsub("CellFractions_PFCref_DWLS.","",names(Allfiles))
names(Allfiles)=gsub(".csv","",names(Allfiles))

names(Allfiles)
Allfiles2 = lapply(Allfiles,melt)



A=Reduce(rbind,Allfiles2)


Fills <- c("#293241","#ee6c4d","#e0fbfc","#98c1d9")
DECON.PLOT=ggplot(A, aes(x = Celltypes, y = Estimated.Celltype.Fraction)) + geom_boxplot(aes(fill = Group)) + scale_fill_manual(values = Fills) + theme_bw()



DeconResults=A %>%
  group_by(Celltypes) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Celltypes, values_from = Estimated.Celltype.Fraction) %>%
  select(-row)
dim(DeconResults)


Wilcoxon.Results <- lapply(2:9, function(x) pairwise.wilcox.test(DeconResults[[x]], DeconResults$Group))
names(Wilcoxon.Results) <- names(DeconResults)[2:9]
Wilcoxon.Results
Wilcoxon.Results.PVALUE=sapply(Wilcoxon.Results, function(x) {
  pValue <- x$p.value
  n <- outer(rownames(pValue), colnames(pValue), paste, sep='v')
  pValue <- as.vector(pValue)
  names(pValue) <- n
  pValue
})


Wilcoxon.Results.names=c("N.NECvN.GAERS","N.NECvE.NEC","N.NECvE.GAERS")


Wilcoxon.Results.B=melt(Wilcoxon.Results.PVALUE)
Wilcoxon.Results.B$FDR=p.adjust(Wilcoxon.Results.B$value, method = "BH", n = length(Wilcoxon.Results.B$value))

colnames(Wilcoxon.Results.B)=c("Comparison","Celltypes","pvalue","FDR")
Wilcoxon.Results.B$log_fdr=-log10(Wilcoxon.Results.B$FDR)

Rownames=c(rownames(Wilcoxon.Results.B)[which(Wilcoxon.Results.B$Comparison == "N.NECvN.GAERS")],rownames(Wilcoxon.Results.B)[which(Wilcoxon.Results.B$Comparison == "N.NECvE.NEC")],
rownames(Wilcoxon.Results.B)[which(Wilcoxon.Results.B$Comparison == "N.NECvE.GAERS")])

Wilcoxon.Results.B2=Wilcoxon.Results.B[intersect(rownames(Wilcoxon.Results.B),Rownames),]
