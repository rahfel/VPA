
######## Scripts to run differential expression analysis ############

load("RData")
colnames(mat_merged)

mat_merged=mat_merged[,-c(10,15)]
GAERS.val2=GAERS.val[-4]
NEC.val2=NEC.val[-2]

Table2=read.csv("S.D.Genes.csv",sep=",",row.names=1)

Table2$Samples=rownames(Table2)
Group=Table2[,c("Samples","Col.","Group")]

samples.id=c(GAERS.val2,GAERS.control)
data.tmp2=mat_merged[,samples.id]

data.tmp=data.tmp2[ rowSums(data.tmp2 > 0) >= 1, ]
dim(data.tmp)

Group.E.GAERS.N.GAERS=rbind(Group[Group$Col. == "E-GAERS",],Group[Group$Col. == "N-GAERS",])
group=c(rep(1,length(GAERS.val2)),rep(0,length(GAERS.control)))
Sex=Group.E.GAERS.N.GAERS$Group
design=model.matrix(~group+Sex)


y <- DGEList(counts=data.tmp,group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef="group")
RES=topTags(lrt,n=nrow(data.tmp))
RES2=RES$table
length(rownames(RES2)[which(RES2$FDR < 0.05)])

######## Scripts to run differential splicing analysis ############

library(edgeR)


Sample.Info.DS=read.csv("S.D.Genes.Sample.Info.DS.csv",sep=",",row.names=1)
rownames(Sample.Info.DS)=Sample.Info.DS$Sample.ID
GTF <- rtracklayer::import('Rattus_norvegicus.Rnor_6.0.94.chr.gtf')
GTF_DF=as.data.frame(GTF)
head(GTF_DF)

N_SE=readRDS("N_SE.rds")

Matrix=N_SE$counts

Sample.Info.DS2=Sample.Info.DS[c(rownames(Sample.Info.DS)[which(Sample.Info.DS$Group == "E-GAERS")],rownames(Sample.Info.DS)[which(Sample.Info.DS$Group == "N-GAERS")]),]

Matrix.E.G.N.G=Matrix[,c(rownames(Sample.Info.DS)[which(Sample.Info.DS$Group == "E-GAERS")],rownames(Sample.Info.DS)[which(Sample.Info.DS$Group == "N-GAERS")])]

rownames(Matrix.E.G.N.G)
y.all <- DGEList(counts=Matrix.E.G.N.G, genes=N_SE$annotation)
dim(y.all)
head(y.all$genes)

Groups2=c(rep(1,7),rep(0,8))

Groups.Keep <- factor(Groups2)
keep <- filterByExpr(y.all, group=Groups.Keep)
table(keep)
y <- y.all[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples
plotMDS(y)

Gender=factor(Sample.Info.DS2$Gender)
design <- model.matrix(~ Groups.Keep+Gender)
design


y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

qlf <- glmQLFTest(fit, coef=2)
topTags(qlf)

is.de <- decideTests(qlf, p.value=0.05)
summary(is.de)


sp <- diffSpliceDGE(fit, coef=2, geneid="GeneID", exonid="Start")

#The top spliced genes under the Simes’ method

Gene.Spliced=topSpliceDGE(sp,test="gene",n=nrow(y))

length(rownames(Gene.Spliced)[which(Gene.Spliced$FDR < 0.05)])

write.csv(Gene.Spliced,"Gene.Spliced.QLF.E.GAERS.N.GAERS.2022.csv")


Simes.Spliced=topSpliceDGE(sp,test="Simes",n=nrow(y))

length(rownames(Simes.Spliced)[which(Simes.Spliced$FDR < 0.05)])

write.csv(Simes.Spliced,"Simes.Spliced.E.GAERS.N.GAERS.2022.csv")
