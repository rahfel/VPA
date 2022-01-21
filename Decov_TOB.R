########Final from loom to creating single cell e18 dataset to generating the deconvolution signature matrix #########

etwd("/rds/general/user/rf1116/ephemeral/DECON/LOOMR/")

library(loomR)
MOUSEBRAIN <- connect(filename = "dev_all.loom", mode = "r+",skip.validate = TRUE)

MOUSEBRAIN

MOUSEBRAIN[["matrix"]]
MOUSEBRAIN[["col_attrs"]]

A=MOUSEBRAIN[["matrix"]][,]
dim(A)

# Pull three bits of metadata from the column attributes
attrs <- c("nUMI", "nGene", "orig.ident")
attr.df <- MOUSEBRAIN$get.attribute.df(MARGIN = 2, attribute.names = attrs)
head(x = attr.df)


dim(A)

Age=MOUSEBRAIN[["col_attrs/CellID"]][]
CellID=MOUSEBRAIN[["col_attrs/CellID"]][]
Class=MOUSEBRAIN[["col_attrs/Class"]][]


Age=MOUSEBRAIN[["col_attrs/CellID"]][]
CellID=MOUSEBRAIN[["col_attrs/CellID"]][]
Class=MOUSEBRAIN[["col_attrs/Class"]][]
Clusters=MOUSEBRAIN[["col_attrs/Clusters"]][]
Tissue=MOUSEBRAIN[["col_attrs/Tissue"]]
table(Tissue)
MOUSEBRAIN.seurat <- as.Seurat(MOUSEBRAIN)


#####################################

setwd("/rds/general/user/rf1116/ephemeral/DECON/LOOMR/")

library(Seurat)
SO=readRDS("MOUSEBRAIN.seurat.rds")

SO.Age.Gene.Region = subset(SO, subset = Region == 'Forebrain' & Age == 'e18.0')
SO.Age.Gene.Region[["percent.mt"]] = PercentageFeatureSet(SO.Age.Gene.Region, pattern = "^mt-")
SO.Age.Gene.Region
SO.Age.Gene.Region2 = subset(SO.Age.Gene.Region, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
SO.Age.Gene.Region2

SO.Matrix=GetAssayData(object = SO.Age.Gene.Region2, assay = "RNA", slot = "data")

GenesDetected = apply(SO.Matrix, 1, function(x) sum(x>0))
table(GenesDetected>=5) # a lot of genes are not detected in 5 or more cells
 # all cells have at least 200 detected genes
keep = GenesDetected>=5
SO.Matrix.B= SO.Matrix[keep,]
dim(SO.Matrix.B)
> dim(SO.Matrix.B)
15237  5530

SO.Final.B= CreateSeuratObject(counts =as.matrix(SO.Matrix.B))




SO.Final.B= NormalizeData(SO.Final.B, normalization.method = "LogNormalize", scale.factor = 10000)
## Find variable features,runPCA etc
SO.Final.B = FindVariableFeatures(SO.Final.B, selection.method = "mean.var.plot", nfeatures = 2000)
All.genes = rownames(SO.Final.B)
SO.Final.B = ScaleData(SO.Final.B, features = All.genes)
SO.Final.B = RunPCA(SO.Final.B,npcs = 30, verbose = FALSE)
SO.Final.B = FindNeighbors(SO.Final.B, dims = 1:10)
SO.Final.B= FindClusters(SO.Final.B, resolution = 2)
SO.Final.B = RunUMAP(SO.Final.B, dims = 1:10)
pdf("SO.Final.UMAP.pdf")
DimPlot(SO.Final.B, reduction = "umap", label = TRUE, pt.size = 1,label.size=5) + NoLegend()
dev.off()

library(dplyr)
SO.Final.B.markers= FindAllMarkers(SO.Final.B, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SO.Final.B.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
SO.Final.B.markers.LIST=unstack(SO.Final.B.markers, SO.Final.B.markers$gene ~ SO.Final.B.markers$cluster)
saveRDS(SO.Final.B.markers.LIST,"SO.Final.B.markers.rds")

SO.Final.B.RawCountsMatrix=GetAssayData(object = SO.Final.B, assay = "RNA", slot = "data")



FisherTest.zeisel=function(x)
{
  TMP=matrix(ncol=9,nrow=1)
  Overlap.Genes=intersect(Filtered.Genes,x)
  for(j in 1:9)
  {
    TMP.Genes=intersect(Filtered.Genes,Zeisel_Mouse.NAMES.List2[[j]])
    TMP.MAT=matrix(ncol=2,nrow=2)
    TMP.MAT[1,1]=length(intersect(TMP.Genes,Overlap.Genes))
    TMP.MAT[1,2]=length(setdiff(TMP.Genes,Overlap.Genes))
    TMP.MAT[2,1]=length(setdiff(Overlap.Genes,TMP.Genes))
    TMP.MAT[2,2]=length(Filtered.Genes)-TMP.MAT[1,1]-TMP.MAT[1,2]-TMP.MAT[2,1]
    TMP[1,j]=fisher.test(TMP.MAT,alternative="greater")$p.value
  }
  TMP
}

CellTypeTest.zeisel=function(x)
{
  CellType=apply(x,1,function(y){colnames(x)[which(y == min(y) & min(y) < 0.05)]})
  CellType.Names=c(rep("Unclassified",length(SO.Final.B.FisherTest.res.zeisel)))
  CellType.Names[which(CellType=="S1.Pyramidal")]="S1.Pyramidal"
  CellType.Names[which(CellType=="CA1.Pyramidal")]="CA1.Pyramidal"
  CellType.Names[which(CellType=="Interneuron")]="Interneuron"
  CellType.Names[which(CellType=="Astrocyte")]="Astrocyte"
  CellType.Names[which(CellType=="Oligodendrocyte")]="Oligodendrocyte"
  CellType.Names[which(CellType=="Microglia")]="Microglia"
  CellType.Names[which(CellType=="Endothelial")]="Endothelial"
  CellType.Names[which(CellType=="Ependymal")]="Ependymal"
  CellType.Names[which(CellType=="Mural")]="Mural"
  CellType.Names
}


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
Zeisel_Mouse.NAMES.List=as.list(read.delim("aaa1934_TableS1.csv",sep=",",header=T))
str(Zeisel_Mouse.NAMES.List)

Zeisel_Mouse.NAMES.List2 <- lapply(Zeisel_Mouse.NAMES.List, function(x) as.factor(firstup(tolower(as.character(x)))))
str(Zeisel_Mouse.NAMES.List2)




# # # #
# Create an object for wang annotation
SO.Final.B.zeisel=SO.Final.B
Filtered.Genes=rownames(SO.Final.B)

SO.Final.B.FisherTest.res.zeisel=lapply(SO.Final.B.markers.LIST,FisherTest.zeisel)

TMP.zeisel = matrix(unlist(SO.Final.B.FisherTest.res.zeisel), ncol = 9, byrow = TRUE)
colnames(TMP.zeisel)=names(Zeisel_Mouse.NAMES.List2)
write.csv(TMP.zeisel,"SO.Final.B.TMP.zeisel.csv")
TMP.LABELS.zeisel=CellTypeTest.zeisel(TMP.zeisel)
names(TMP.LABELS.zeisel)=names(SO.Final.B.FisherTest.res.zeisel)
wang.cluster.combined.ids = TMP.LABELS.zeisel
names(wang.cluster.combined.ids) = levels(SO.Final.B.zeisel)
SO.Final.B.zeisel = RenameIdents(SO.Final.B.zeisel,wang.cluster.combined.ids)
levels(SO.Final.B.zeisel)
pdf("SO.Final.B.LABELLED.pdf")
DimPlot(SO.Final.B.zeisel, reduction = "umap", label = TRUE, pt.size = 1,label.size=5) + NoLegend()
dev.off()
table(Idents(SO.Final.B.zeisel))
levels(SO.Final.B.zeisel)

SO.Final.B.zeisel[["ClusterIdent"]] = Idents(object = SO.Final.B.zeisel)

Meta=SO.Final.B.zeisel@meta.data

Meta2=Meta[!grepl("Unclassified", Meta$ClusterIdent),]

ClusterIdent=as.character(SO.Final.B.zeisel@meta.data$ClusterIdent)
names(ClusterIdent)=rownames(SO.Final.B.zeisel@meta.data)

saveRDS(as.factor(ClusterIdent),"SO.Final.B.CellID_Type.rds")

#Extract counts only for filtered matrix reads
SO.Final.B.FinalMatrix=GetAssayData(object = SO.Final.B.zeisel, assay = "RNA", slot = "data")

SO.Final.B.FinalMatrix1=SO.Final.B.FinalMatrix[,intersect(colnames(SO.Final.B.FinalMatrix),rownames(Meta2))]
SO.Final.B.FinalMatrix2=SO.Final.B.RawCountsMatrix[,colnames(SO.Final.B.FinalMatrix1)]
dim(SO.Final.B.FinalMatrix2)
dim(SO.Final.B.FinalMatrix)

saveRDS(SO.Final.B.FinalMatrix2,"SO.Final.B.FinalMatrix2.nonNormalised.dgCMatrix.rds")

saveRDS(Meta2,"SO.Final.B.FinalMatrix2.MetaData.rds")


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################


#Deconvolution Functions
#install these packages if necessary...

#load packages
library(quadprog)
library(reshape)
library(e1071)
library(Seurat)
library(ROCR)
library(varhandle)
library(MAST)
library(stringr)


#trim bulk and single-cell data to contain the same genes
trimData<-function(Signature,bulkData){
	Genes<-intersect(rownames(Signature),names(bulkData))
  	B<-bulkData[Genes]
  	S<-Signature[Genes,]
  	return(list("sig"=S,"bulk"=B))
}


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

#build signature matrix using genes identified by DEAnalysis()
buildSignatureMatrixUsingSeurat<-function(scdata,id,path,diff.cutoff=0.5,pval.cutoff=0.01){

  #perform differential expression analysis
  DEAnalysis(scdata,id,path)

  numberofGenes<-c()
  for (i in unique(id)){
    load(file=paste(path,"/de_",i,".RData",sep=""))
    DEGenes<-rownames(de_group)[intersect(which(de_group$p_val_adj<pval.cutoff),which(de_group$avg_logFC>diff.cutoff))]
    nonMir = grep("MIR|Mir", DEGenes, invert = T)
    assign(paste("cluster_lrTest.table.",i,sep=""),de_group[which(rownames(de_group)%in%DEGenes[nonMir]),])
    numberofGenes<-c(numberofGenes,length(DEGenes[nonMir]))
  }

  #need to reduce number of genes
  #for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
  #choose matrix with lowest condition number
  conditionNumbers<-c()
  for(G in 50:200){
    Genes<-c()
    j=1
    for (i in unique(id)){
      if(numberofGenes[j]>0){
        temp<-paste("cluster_lrTest.table.",i,sep="")
        temp<-as.name(temp)
        temp<-eval(parse(text = temp))
        temp<-temp[order(temp$p_val_adj,decreasing=TRUE),]
        Genes<-c(Genes,(rownames(temp)[1:min(G,numberofGenes[j])]))
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
  G<-which.min(conditionNumbers)+min(49,numberofGenes-1) #G is optimal gene number
  #
  Genes<-c()
  j=1
  for (i in unique(id)){
    if(numberofGenes[j]>0){
      temp<-paste("cluster_lrTest.table.",i,sep="")
      temp<-as.name(temp)
      temp<-eval(parse(text = temp))
      temp<-temp[order(temp$p_val_adj,decreasing=TRUE),]
      Genes<-c(Genes,(rownames(temp)[1:min(G,numberofGenes[j])]))
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

##alternative differential expression method using MAST

#functions for DE

Mean.in.log2space=function(x,pseudo.count) {
  return(log2(mean(2^(x)-pseudo.count)+pseudo.count))
}

stat.log2=function(data.m, group.v, pseudo.count){
  #data.m=data.used.log2
  log2.mean.r <- aggregate(t(data.m), list(as.character(group.v)), function(x) Mean.in.log2space(x,pseudo.count))
  log2.mean.r <- t(log2.mean.r)
  colnames(log2.mean.r) <- paste("mean.group",log2.mean.r[1,], sep="")
  log2.mean.r = log2.mean.r[-1,]
  log2.mean.r = as.data.frame(log2.mean.r)
  log2.mean.r = varhandle::unfactor(log2.mean.r)  #from varhandle
  log2.mean.r[,1] = as.numeric(log2.mean.r[,1])
  log2.mean.r[,2] = as.numeric(log2.mean.r[,2])
  log2_foldchange = log2.mean.r$mean.group1-log2.mean.r$mean.group0
  results = data.frame(cbind(log2.mean.r$mean.group0,log2.mean.r$mean.group1,log2_foldchange))
  colnames(results) = c("log2.mean.group0","log2.mean.group1","log2_fc")
  rownames(results) = rownames(log2.mean.r)
  return(results)
}

v.auc = function(data.v,group.v) {
  prediction.use=prediction(data.v, group.v, 0:1)
  perf.use=performance(prediction.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  return(auc.use)
}
m.auc=function(data.m,group.v) {
  AUC=apply(data.m, 1, function(x) v.auc(x,group.v))
  AUC[is.na(AUC)]=0.5
  return(AUC)

}

#perform DE analysis using MAST
DEAnalysisMAST<-function(scdata,id,path){

  pseudo.count = 0.1
  data.used.log2   <- log2(scdata+pseudo.count)
  colnames(data.used.log2)<-make.unique(colnames(data.used.log2))
  diff.cutoff=0.5
  for (i in unique(id)){
    cells.symbol.list2     = colnames(data.used.log2)[which(id==i)]
    cells.coord.list2      = match(cells.symbol.list2, colnames(data.used.log2))
    cells.symbol.list1     = colnames(data.used.log2)[which(id != i)]
    cells.coord.list1      = match(cells.symbol.list1, colnames(data.used.log2))
    data.used.log2.ordered  = cbind(data.used.log2[,cells.coord.list1], data.used.log2[,cells.coord.list2])
    group.v <- c(rep(0,length(cells.coord.list1)), rep(1, length(cells.coord.list2)))
    #ouput
    log2.stat.result <- stat.log2(data.used.log2.ordered, group.v, pseudo.count)
    Auc <- m.auc(data.used.log2.ordered, group.v)
    bigtable <- data.frame(cbind(log2.stat.result, Auc))

    DE <- bigtable[bigtable$log2_fc >diff.cutoff,]
    dim(DE)
    if(dim(DE)[1]>1){
      data.1                 = data.used.log2[,cells.coord.list1]
      data.2                 = data.used.log2[,cells.coord.list2]
      genes.list = rownames(DE)
      log2fold_change        = cbind(genes.list, DE$log2_fc)
      colnames(log2fold_change) = c("gene.name", "log2fold_change")
      counts  = as.data.frame(cbind( data.1[genes.list,], data.2[genes.list,] ))
      groups  = c(rep("Cluster_Other", length(cells.coord.list1) ), rep(i, length(cells.coord.list2) ) )
      groups  = as.character(groups)
      data_for_MIST <- as.data.frame(cbind(rep(rownames(counts), dim(counts)[2]), melt(counts),rep(groups, each = dim(counts)[1]), rep(1, dim(counts)[1] * dim(counts)[2]) ))
      colnames(data_for_MIST) = c("Gene", "Subject.ID", "Et", "Population", "Number.of.Cells")
      vbeta = data_for_MIST
      vbeta.fa <-FromFlatDF(vbeta, idvars=c("Subject.ID"),
                            primerid='Gene', measurement='Et', ncells='Number.of.Cells',
                            geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'),
                            phenovars=c('Population'), id='vbeta all')
      vbeta.1 <- subset(vbeta.fa,Number.of.Cells==1)
      # .3 MAST
      head(colData(vbeta.1))
      zlm.output <- zlm(~ Population, vbeta.1, method='bayesglm', ebayes=TRUE)
      show(zlm.output)
      coefAndCI <- summary(zlm.output, logFC=TRUE)
      zlm.lr <- lrTest(zlm.output, 'Population')
      zlm.lr_pvalue <- melt(zlm.lr[,,'Pr(>Chisq)'])
      zlm.lr_pvalue <- zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == 'hurdle'),]



      lrTest.table <-  merge(zlm.lr_pvalue, DE, by.x = "primerid", by.y = "row.names")
      colnames(lrTest.table) <- c("Gene", "test.type", "p_value", paste("log2.mean.", "Cluster_Other", sep=""), paste("log2.mean.",i,sep=""), "log2fold_change", "Auc")
      cluster_lrTest.table <- lrTest.table[rev(order(lrTest.table$Auc)),]

      #. 4 save results
      write.csv(cluster_lrTest.table, file=paste(path,"/",i,"_lrTest.csv", sep=""))
      save(cluster_lrTest.table, file=paste(path,"/",i,"_MIST.RData", sep=""))
    }
  }
}



#build signature matrix using genes identified by DEAnalysisMAST()
buildSignatureMatrixMAST<-function(scdata,id,path,diff.cutoff=0.5,pval.cutoff=0.01){
  #compute differentially expressed genes for each cell type
  DEAnalysisMAST(scdata,id,path)

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

  #need to reduce number of genes
  #for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
  #for each, iterate and choose matrix with lowest condition number
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

data <- readRDS("SO.Final.B.FinalMatrix2.nonNormalised.dgCMatrix.rds")

data <- as.data.frame(readRDS("SO.Final.B.FinalMatrix2.nonNormalised.dgCMatrix.rds"))
dim(data)
meta <- readRDS("SO.Final.B.FinalMatrix2.MetaData.rds")
dim(meta)

colnames(meta)

labels <- as.vector(meta$ClusterIdent)


labels <- as.vector(meta$CellType) # Grab celltype labels from metadata celltype column
Signature<-buildSignatureMatrixMAST(scdata=data,
                                    id=labels,
                                    path="DWLS.results",
                                    diff.cutoff=0.5,
                                    pval.cutoff=0.01)

save(Signature,file="DWLS.RData")

# 2 step process of making your mouse gene symbols (to avoding a formatting error)
Name <- data$mmusculus_homolog_associated_gene_name
rownames(CPM) <- Name
#

CPM.trim <- CPM[rownames(Signature),] # filters your counts for genes maintained in the Signature Matrix

samples<- colnames(CPM.trim) # grabbing a vector of your subject names

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
