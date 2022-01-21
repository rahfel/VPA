

library(edgeR)
GTF <- rtracklayer::import('Rattus_norvegicus.Rnor_6.0.94.chr.gtf')
GTF_DF=as.data.frame(GTF)
head(GTF_DF)

library(edgeR)

N.BAMFiles=c("/Users/rafel137/Documents/Aus_Samples2/IGF106022.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106024.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106028.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106038.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106046.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106047.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106051.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106021.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106031.bam"
,"/Users/rafel137/Documents/Aus_Samples2/IGF106033.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106040.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106042.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106045.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106048.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106020.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106026.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106030.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106032.bam"
,"/Users/rafel137/Documents/Aus_Samples2/IGF106036.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106041.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106044.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106049.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106043.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106023.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106025.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106027.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106035.bam"
,"/Users/rafel137/Documents/Aus_Samples2/IGF106037.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106039.bam","/Users/rafel137/Documents/Aus_Samples2/IGF106050.bam")


#N_SE <- featureCounts(N.BAMFiles, annot.ext="/Users/rafel137/Documents/Aus_Samples2/Rattus_norvegicus.Rnor_6.0.94.chr.gtf",isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id",useMetaFeatures=FALSE, allowMultiOverlap=TRUE,isPairedEnd=T)
#saveRDS(N_SE,"N_SE.rds")


N_SE=readRDS("N_SE.rds")

y.all <- DGEList(counts=N_SE$counts, genes=N_SE$annotation)
dim(y.all)
head(y.all$genes)


Groups2=c(rep(1,14),rep(0,16))

y=y.all

y$samples

Groups.Keep <- factor(Groups2)
keep <- filterByExpr(y, group=Groups.Keep)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples
plotMDS(y)


design <- model.matrix(~ Groups.Keep)
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

#The top spliced genes under the Simes’ method are shown below
topSpliceDGE(sp, test="Simes", n=30)

topSpliceDGE(sp, test="gene", n=30)


#par(mfrow=c(1,2))
plotSpliceDGE(sp, geneid="ENSRNOG00000022804", genecol="GeneID")
plotSpliceDGE(sp, geneid="ENSRNOG00000018200", genecol="GeneID")


Gene.Spliced=topSpliceDGE(sp,test="gene",n=nrow(y))

length(rownames(Gene.Spliced)[which(Gene.Spliced$FDR < 0.05)])

write.csv(Gene.Spliced,"Gene.Spliced.QLF.E.N.csv")

Simes.Spliced=topSpliceDGE(sp,test="Simes",n=nrow(y))

length(rownames(Simes.Spliced)[which(Simes.Spliced$FDR < 0.05)])

write.csv(Simes.Spliced,"Simes.Spliced2.E.N.csv")




#Two testing methods at the gene-level are provided. The Simes’ method is likely to be more
#powerful when only a minority of the exons for a gene are differentially spliced. The F-tests
#are likely to be powerful for genes in which several exons are differentially spliced.

Gene.Spliced=topSpliceDGE(sp,test="gene",n=nrow(y))

length(rownames(Gene.Spliced)[which(Gene.Spliced$FDR < 0.05)])

write.csv(Gene.Spliced,"Gene.Spliced.E.N.csv")

Simes.Spliced=topSpliceDGE(sp,test="Simes",n=nrow(y))

length(rownames(Simes.Spliced)[which(Simes.Spliced$FDR < 0.05)])

write.csv(Simes.Spliced,"Simes.Spliced.E.N.csv")



################################################################################
#### Create a gene sets and entr HUMAN ID ### STEP1: GET ENTID from biomart (Rat to Human) ## STEP2: Create a csv first then convert it to txt and open VI and replace comma with blank

##cat AS.Genes.ENTID.csv | tr  ',' '\n' > AS.Genes.ENTID2.txt

#vi AS.Genes.ENTID2.txt + %s/,/ /g



PATH=$PATH:/Users/rafel137/Downloads/magma_v1/
Gene_annot_out=/Users/rafel137/Downloads/magma_v1/PD2019_W35.W10.Adult_Brain.magma.genes.raw
DEGenes_File=/Users/rafel137/Documents/Aus_Samples2/AS.Genes.ENTID2.txt

magma --gene-results $Gene_annot_out --set-annot $DEGenes_File 'col=1,2' --out PD2019.Adult_Brain.AS.Genes.ENTID2_AllSpliced_ENTID_DEgenes

C=read.table("FATSL2018.Adult_Brain.AS.Genes.ENTID2_AllSpliced_ENTID_DEgenes.gsa.out",sep="",header=T)
C$FDR=p.adjust(C$P, method = "BH", n = length(C$P))
C
rownames(C)=C$VARIABLE
write.table(C,"FATSL2018.Adult_Brain.B_C_TWENTYFIVE_AllRegulated_ENTID_DEgenes.gsa.out.txt")


########

Gene.Spliced.GeneID=Gene.Spliced$GeneID[which(Gene.Spliced$FDR < 0.05)]
saveRDS(Gene.Spliced.GeneID,"Gene.Spliced.GeneID.rds")

setwd("/Users/rafel137/Documents/Aus_Samples2")

Gene.Spliced.GeneID=readRDS("Gene.Spliced.GeneID.rds")
library(ggfortify)
library(DESeq2)
library("factoextra")
library(ggfortify)
library(DESeq2)
library(ggfortify)
library(DESeq2)
library("factoextra")
mat_merged2=readRDS("mat_merged.rds")
colnames(mat_merged2)
mat_merged=mat_merged2[,-c(10,15)]
colnames(mat_merged)
dim(mat_merged)

rlogcounts <- rlog(mat_merged)
library(reshape2)
load("RData")


GAERS.val2=GAERS.val[-4]
NEC.val2=NEC.val[-2]
A=list(GAERS.val2,GAERS.control,NEC.val2,NEC.control)
names(A)=c("E-GAERS","N-GAERS","E-NEC","N-NEC")
groups=melt(A)
colnames(groups)=c("SampleID","Groups")
groups2 = groups[with(groups, order(SampleID)), ]
rownames(groups2) = 1:nrow(groups2)
groups2
groups2B=as.factor(groups2$Groups)

# run PCA
pcDat <- prcomp(t(rlogcounts))
# plot PCA
autoplot(pcDat)
fviz_pca_ind(pcDat)

fviz_pca_ind(pcDat,
             col.ind = groups2B, # color by groups
             palette = c("black","grey","lightpink","red4"),
             addEllipses = F, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
             )


###

DataA=readRDS("mat_merged.rds")
#use = rowSums(DataA > 1) >=16
Data = DataA[ ,-c(10,15) ]


cormM = cor(log2(Data+1), method = "spearman")
heatmap.MIN.Val = max(0, (trunc(min(cormM) * 10) - 1) / 9)
heatmap.MAX.Val = 1
heatmap.RANGE.Val = (heatmap.MAX.Val - heatmap.MIN.Val) / 100
grouplabelcol = rainbow(ncol(Data), s=1, alpha = 1)
names(grouplabelcol)=c(colnames(Data))
anno_colors <- list(pheno = grouplabelcol)



groups2$Exposure=groups2$Groups
groups2$Exposure=gsub("N-GAERS","N",groups2$Exposure)
groups2$Exposure=gsub("E-GAERS","E",groups2$Exposure)
groups2$Exposure=gsub("N-NEC","N",groups2$Exposure)
groups2$Exposure=gsub("E-NEC","E",groups2$Exposure)


ColAnnot=data.frame(Exposure=groups2$Exposure)
rownames(ColAnnot)=groups2$SampleID
RowAnnot=data.frame(Group=groups2$Group)
rownames(RowAnnot)=groups2$SampleID


anno_colors2 = list(Exposure = c("E"="black", "N"="white"),Groups = c("N-GAERS"="lightgrey","E-GAERS"="darkgrey","N-NEC"="lightpink","E-GAERS"="salmon"))
pheatmap(cormM,
         breaks = seq(heatmap.MIN.Val, heatmap.MAX.Val, heatmap.RANGE.Val),
         fontsize = 15,
         fontsize_row=12,
         fontsize_col=12,
         annotation_colors = anno_colors2,
         annotation_col    =ColAnnot,
         annotation_row    =RowAnnot,
         cluster_cols = T,
         cluster_rows = T)

#

Gene.Spliced.GeneID2=c("ENSRNOG00000014104","ENSRNOG00000031997","ENSRNOG00000018220","ENSRNOG00000054058","ENSRNOG00000020989","ENSRNOG00000001557","ENSRNOG00000003081",
"ENSRNOG00000005900","ENSRNOG00000018992","ENSRNOG00000018602","ENSRNOG00000000525","ENSRNOG00000002544","ENSRNOG00000020906","ENSRNOG00000001439",
"ENSRNOG00000038572","ENSRNOG00000038572","ENSRNOG00000003657","ENSRNOG00000050071","ENSRNOG00000014287","ENSRNOG00000018184","ENSRNOG00000013583",
"ENSRNOG00000010589","ENSRNOG00000010838","ENSRNOG00000012579","ENSRNOG00000018778","ENSRNOG00000006617")


Data = DataA[ Gene.Spliced.GeneID2,-c(10,15) ]
Data2=as.data.frame(t(Data))

heatmap(Data)

#
aggdata <-aggregate(mtcars, by=list(cyl,vs),
  FUN=mean, na.rm=TRUE)

anno_colors2 = list(Exposure = c("E"="black", "N"="white"),Groups = c("N-GAERS"="lightgrey","E-GAERS"="darkgrey","N-NEC"="lightpink","E-GAERS"="salmon"))
pheatmap(log2(Data+1),
         fontsize = 15,
         fontsize_row=12,
         fontsize_col=12,
         annotation_colors = anno_colors2,
         annotation_col    =ColAnnot,
         cluster_cols = T,
         cluster_rows = T)

######

setwd("/Users/rafel137/Documents/Aus_Samples2")

Gene.Spliced.GeneID=readRDS("Gene.Spliced.GeneID.rds")
library(ggfortify)
library(DESeq2)
library("factoextra")
library(ggfortify)
library(DESeq2)
library(ggfortify)
library(DESeq2)
library("factoextra")
mat_merged2=readRDS("mat_merged.rds")
colnames(mat_merged2)
mat_merged=mat_merged2[,-c(10,15)]
colnames(mat_merged)
dim(mat_merged)
load("RData")
library(gsubfn)
GAERS.val2=GAERS.val[-4]
NEC.val2=NEC.val[-2]

#ENSRNOG00000022804	Foxp4
#ENSRNOG00000014104  Myo5b
#ENSRNOG00000050071	Cdc45
#ENSRNOG00000018602	Camta1
#ENSRNOG00000026649	Dnmt3a
#
#ENSRNOG00000014811 Cyhr1
#

Names=c("ENSRNOG00000022804","ENSRNOG00000014104","ENSRNOG00000018200","ENSRNOG00000018602","ENSRNOG00000014811","ENSRNOG00000053122")     # create list of names



mat_merged2=readRDS("mat_merged.rds")
colnames(mat_merged2)
mat_merged=mat_merged2[Names,-c(10,15)]

NEC=mat_merged[Names,NEC.val2]
NEC.N=mat_merged[Names,NEC.control]
GAERS=mat_merged[Names,GAERS.val2]
GAERS.N=mat_merged[Names,GAERS.control]

A=list(NEC,GAERS,NEC.N,GAERS.N)
names(A)=c("E","E","N","N")

A2=melt(A)

#A2$Var1=gsub("ENSRNOG00000022804","Foxp4",A2$Var1) 0.498018511190823
#A2$Var1=gsub("ENSRNOG00000014104","Myo5b",A2$Var1) 7.54024888896506E-09
#A2$Var1=gsub("ENSRNOG00000018200","Gad2",A2$Var1) 0.310036247705586
#A2$Var1=gsub("ENSRNOG00000018602","Camta1",A2$Var1) 0.000583148229070548
#A2$Var1=gsub("ENSRNOG00000014811","Cyhr1",A2$Var1) 0.150978253628741
#A2$Var1=gsub("ENSRNOG00000053122","SCN1A",A2$Var1) 0.150978253628741

ENSRNOG00000053122

A2$Var1=gsub("ENSRNOG00000022804","Foxp4",A2$Var1)
A2$Var1=gsub("ENSRNOG00000014104","Myo5b",A2$Var1)
A2$Var1=gsub("ENSRNOG00000018200","Gad2",A2$Var1)
A2$Var1=gsub("ENSRNOG00000018602","Camta1",A2$Var1)
A2$Var1=gsub("ENSRNOG00000014811","Cyhr1",A2$Var1)
A2$Var1=gsub("ENSRNOG00000053122","SCN1A",A2$Var1)


A2$Var2=NULL

Fills <- c("#293241","#ee6c4d")
DECON.PLOT=ggplot(A2, aes(x = Var1, y = value)) + geom_boxplot(aes(fill = L1)) + scale_fill_manual(values = Fills) + theme_bw()


#



DeconResults=A2 %>%
  group_by(Var1) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Var1, values_from = value) %>%
  select(-row)
dim(DeconResults)


Wilcoxon.Results <- lapply(2:7, function(x) pairwise.wilcox.test(DeconResults[[x]], DeconResults$L1))
names(Wilcoxon.Results) <- names(DeconResults)[2:7]
Wilcoxon.Results
Wilcoxon.Results.PVALUE=sapply(Wilcoxon.Results, function(x) {
  pValue <- x$p.value
  n <- outer(rownames(pValue), colnames(pValue), paste, sep='v')
  pValue <- as.vector(pValue)
  names(pValue) <- n
  pValue
})




Wilcoxon.Results.B=melt(Wilcoxon.Results.PVALUE)
Wilcoxon.Results.B$FDR=p.adjust(Wilcoxon.Results.B$value, method = "BH", n = length(Wilcoxon.Results.B$value))

Wilcoxon.Results.B$log_fdr=-log10(Wilcoxon.Results.B$FDR)

Wilcoxon.Results.B$Genes=rownames( Wilcoxon.Results.B)

Wilcoxon.Results.B$Genes=gsub(".NvE","",Wilcoxon.Results.B$Genes)

Wilcoxon.Results.PLOT=ggplot(Wilcoxon.Results.B,aes(x = Genes, y = log_fdr, fill = log_fdr)) + geom_tile(colour = "black") +
labs(x = "Genes", y = "", fill = "-log"[1][0]~"(FDR)") + scale_fill_gradient(low = "#dddddd", high = "#777777", na.value = "white") + geom_text(aes(label = case_when(log_fdr > 1.30103 ~ "**"))) +
theme(legend.key.size = unit(0.5, "cm"), legend.position = "right", plot.margin = margin(t = 0, r = 0.5, b = 0.5, l = 0.5, "cm"))



FigList=list(DECON.PLOT,Wilcoxon.Results.PLOT)

cowplot::plot_grid(plotlist = FigList,
                           ncol = 1,
                           align = "v", axis = "lr",
                           rel_heights = c(1, 0.65))










gggg


setwd("/Users/rafel137/Documents/Aus_Samples2")

Gene.Spliced.GeneID=readRDS("Gene.Spliced.GeneID.rds")
library(ggfortify)
library(DESeq2)
library("factoextra")
library(ggfortify)
library(DESeq2)
library(ggfortify)
library(DESeq2)
library("factoextra")
mat_merged2=readRDS("mat_merged.rds")
colnames(mat_merged2)
mat_merged=mat_merged2[,-c(10,15)]
colnames(mat_merged)
dim(mat_merged)
load("RData")
library(gsubfn)
GAERS.val2=GAERS.val[-4]
NEC.val2=NEC.val[-2]

#ENSRNOG00000022804	Foxp4
#ENSRNOG00000014104  Myo5b
#ENSRNOG00000050071	Cdc45
#ENSRNOG00000018602	Camta1
#ENSRNOG00000026649	Dnmt3a
#
#ENSRNOG00000014811 Cyhr1
#

DEGs=read.csv("TB2021.E.NEC.vs.N.NEC.csv",sep=",",row.names=1)
M30=read.csv("M30.txt",header =F)

DEGs2=DEGs[intersect(rownames(DEGs),M30$V1),]
dim(DEGs2)
DEGs2$Log10.Pval=-log10(DEGs2$PValue)
plot(DEGs2$logFC,DEGs2$Log10.Pval)


M30_DEGs <- DEGs2 %>%
             filter(FDR < 0.05)

DEGs2 %>%
  ggplot(aes(x=logFC,y=Log10.Pval)) +
  geom_point(alpha=0.3) +
  geom_point(data=M30_DEGs,
             aes(x=logFC,y=Log10.Pval),
             color='red',
             size=1)+ xlab("Log2 fold change") + ylab("-log10(P.Value)")





pheatmap(test)


mat_merged2=readRDS("mat_merged.rds")
colnames(mat_merged2)
Rownames=intersect(rownames(mat_merged2),M30$V1)
length(Rownames)

NEC=c(NEC.val2,NEC.control,GAERS.val2,GAERS.control)

NEC.MX=as.matrix(mat_merged2[M30$V1,NEC])

annotation_col = data.frame(
    Group = factor(c(rep(c("E.NEC"),7),rep(c("N.NEC"),8),rep(c("E.GAERS"),7),rep(c("N.GAERS"),8))))
rownames(annotation_col) = NEC
annotation_col


Group = c(E.NEC = "blue", N.NEC = "lightblue",E.GAERS="red",N.GAERS="lightpink")

library(pheatmap)

INT.Genes=c("ENSRNOG00000054508","ENSRNOG00000053122","ENSRNOG00000055401","ENSRNOG00000060599","ENSRNOG00000011568",
"ENSRNOG00000003680","ENSRNOG00000007377","ENSRNOG00000005206","ENSRNOG00000005509",
"ENSRNOG00000024650","ENSRNOG00000058539","ENSRNOG00000011568","ENSRNOG00000039949")

NEC.MX2=NEC.MX[INT.Genes,]

rownames(NEC.MX2)=gsub("ENSRNOG00000054508","Foxp2",rownames(NEC.MX2))
rownames(NEC.MX2)=gsub("ENSRNOG00000053122","Scn1a",rownames(NEC.MX2))
rownames(NEC.MX2)=gsub("ENSRNOG00000055401","Kcnc1",rownames(NEC.MX2))
rownames(NEC.MX2)=gsub("ENSRNOG00000060599","Gabrb3",rownames(NEC.MX2))
rownames(NEC.MX2)=gsub("ENSRNOG00000011568","Rspo3",rownames(NEC.MX2))
rownames(NEC.MX2)=gsub("ENSRNOG00000003680","Gabrb2",rownames(NEC.MX2))
rownames(NEC.MX2)=gsub("ENSRNOG00000007377","Slit3",rownames(NEC.MX2))
rownames(NEC.MX2)=gsub("ENSRNOG00000005206","Kcnq3",rownames(NEC.MX2))
rownames(NEC.MX2)=gsub("ENSRNOG00000005509","Pak7",rownames(NEC.MX2))
rownames(NEC.MX2)=gsub("ENSRNOG00000024650","Ckap2",rownames(NEC.MX2))
rownames(NEC.MX2)=gsub("ENSRNOG00000058539","Ccnb1",rownames(NEC.MX2))
rownames(NEC.MX2)=gsub("ENSRNOG00000011568","Rspo3",rownames(NEC.MX2))
rownames(NEC.MX2)=gsub("ENSRNOG00000039949","Rad54b",rownames(NEC.MX2))


pheatmap(log2(NEC.MX2+1), annotation_col = annotation_col)


ENSRNOG00000054508	Foxp2
ENSRNOG00000053122	Scn1a
ENSRNOG00000055401	Kcnc1
ENSRNOG00000060599	Gabrb3
ENSRNOG00000011568	Rspo3
ENSRNOG00000003680	Gabrb2
ENSRNOG00000007377	Slit3
ENSRNOG00000005206	Kcnq3
ENSRNOG00000005509	Pak7

ENSRNOG00000024650	Ckap2
ENSRNOG00000058539	Ccnb1
ENSRNOG00000011568	Rspo3
ENSRNOG00000039949	Rad54b



mat_merged2=readRDS("mat_merged.rds")
colnames(mat_merged2)
mat_merged=mat_merged2[,-c(10,15)]

NEC=mat_merged[INT.Genes,NEC.val2]
NEC.N=mat_merged[INT.Genes,NEC.control]
GAERS=mat_merged[INT.Genes,GAERS.val2]
GAERS.N=mat_merged[INT.Genes,GAERS.control]

A=list(NEC,NEC.N,GAERS,GAERS.N)
names(A)=c("E.NEC","N.NEC","E.GAERS","N.GAERS")

A2=melt(A)

A2$Var1=gsub("ENSRNOG00000054508","Foxp2",A2$Var1)
A2$Var1=gsub("ENSRNOG00000053122","Scn1a",A2$Var1)
A2$Var1=gsub("ENSRNOG00000055401","Kcnc1",A2$Var1)
A2$Var1=gsub("ENSRNOG00000060599","Gabrb3",A2$Var1)
A2$Var1=gsub("ENSRNOG00000011568","Rspo3",A2$Var1)
A2$Var1=gsub("ENSRNOG00000003680","Gabrb2",A2$Var1)
A2$Var1=gsub("ENSRNOG00000007377","Slit3",A2$Var1)
A2$Var1=gsub("ENSRNOG00000005206","Kcnq3",A2$Var1)
A2$Var1=gsub("ENSRNOG00000005509","Pak7",A2$Var1)
A2$Var1=gsub("ENSRNOG00000024650","Ckap2",A2$Var1)
A2$Var1=gsub("ENSRNOG00000058539","Ccnb1",A2$Var1)
A2$Var1=gsub("ENSRNOG00000011568","Rspo3",A2$Var1)
A2$Var1=gsub("ENSRNOG00000039949","Rad54b",A2$Var1)



A2$Var2=NULL

Fills <- c("blue","lightblue","red","lightpink")
ggplot(A2, aes(x = Var1, y = log2(value))) + geom_boxplot(aes(fill = L1)) + scale_fill_manual(values = Fills) + theme_bw() +
 scale_x_discrete(name ="Selected.M30genes") + scale_y_discrete(name ="log2(Expression.Values)")


#



DeconResults=A2 %>%
  group_by(Var1) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Var1, values_from = value) %>%
  select(-row)
dim(DeconResults)


Wilcoxon.Results <- lapply(2:7, function(x) pairwise.wilcox.test(DeconResults[[x]], DeconResults$L1))
names(Wilcoxon.Results) <- names(DeconResults)[2:7]
Wilcoxon.Results
Wilcoxon.Results.PVALUE=sapply(Wilcoxon.Results, function(x) {
  pValue <- x$p.value
  n <- outer(rownames(pValue), colnames(pValue), paste, sep='v')
  pValue <- as.vector(pValue)
  names(pValue) <- n
  pValue
})




Wilcoxon.Results.B=melt(Wilcoxon.Results.PVALUE)
Wilcoxon.Results.B$FDR=p.adjust(Wilcoxon.Results.B$value, method = "BH", n = length(Wilcoxon.Results.B$value))

Wilcoxon.Results.B$log_fdr=-log10(Wilcoxon.Results.B$FDR)

Wilcoxon.Results.B$Genes=rownames( Wilcoxon.Results.B)

Wilcoxon.Results.B$Genes=gsub(".NvE","",Wilcoxon.Results.B$Genes)

Wilcoxon.Results.PLOT=ggplot(Wilcoxon.Results.B,aes(x = Genes, y = log_fdr, fill = log_fdr)) + geom_tile(colour = "black") +
labs(x = "Genes", y = "", fill = "-log"[1][0]~"(FDR)") + scale_fill_gradient(low = "#dddddd", high = "#777777", na.value = "white") + geom_text(aes(label = case_when(log_fdr > 1.30103 ~ "**"))) +
theme(legend.key.size = unit(0.5, "cm"), legend.position = "right", plot.margin = margin(t = 0, r = 0.5, b = 0.5, l = 0.5, "cm"))



FigList=list(DECON.PLOT,Wilcoxon.Results.PLOT)

cowplot::plot_grid(plotlist = FigList,
                           ncol = 1,
                           align = "v", axis = "lr",
                           rel_heights = c(1, 0.65))


######################
###############################################################################
#################

N_SE=readRDS("N_SE.rds")

y.all <- DGEList(counts=N_SE$counts, genes=N_SE$annotation)
dim(y.all)
head(y.all$genes)


Groups2=c(rep(1,14),rep(0,16))

y=y.all

y$samples

Groups.Keep <- factor(Groups2)
keep <- filterByExpr(y, group=Groups.Keep)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples
plotMDS(y)


design <- model.matrix(~ Groups.Keep)
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

#The top spliced genes under the Simes’ method are shown below
topSpliceDGE(sp, test="Simes", n=30)

topSpliceDGE(sp, test="gene", n=30)


#par(mfrow=c(1,2))
plotSpliceDGE(sp, geneid="ENSRNOG00000022804", genecol="GeneID")
plotSpliceDGE(sp, geneid="ENSRNOG00000018200", genecol="GeneID")


Gene.Spliced=topSpliceDGE(sp,test="gene",n=nrow(y))

length(rownames(Gene.Spliced)[which(Gene.Spliced$FDR < 0.05)])

write.csv(Gene.Spliced,"Gene.Spliced.QLF.E.N.csv")

Simes.Spliced=topSpliceDGE(sp,test="Simes",n=nrow(y))

length(rownames(Simes.Spliced)[which(Simes.Spliced$FDR < 0.05)])

write.csv(Simes.Spliced,"Simes.Spliced.QLF.E.N.csv")


###################################################################################################################
#################

N_SE=readRDS("N_SE.rds")

y.all <- DGEList(counts=N_SE$counts, genes=N_SE$annotation)
dim(y.all)
head(y.all$genes)


Groups2=c(rep(1,14),rep(0,16))

y=y.all

y$samples

Groups.Keep <- factor(Groups2)
keep <- filterByExpr(y, group=Groups.Keep)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples

plotMDS(y, col=c(rep(c(1),7),rep(c(2),7),rep(c(3),8),rep(c(4),8)))


##legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)


design <- model.matrix(~ Groups.Keep)
design


y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)

fit <- glmFit(y, design, robust=TRUE)

lrt <- glmLRT(fit, coef=2)
topTags(lrt)

is.de <- decideTests(lrt, p.value=0.05)
summary(is.de)

plotMD(lrt)

abline(h=c(-1, 1), col="blue")

sp <- diffSpliceDGE(fit, coef=2, geneid="GeneID", exonid="Start")

#The top spliced genes under the Simes’ method are shown below
topSpliceDGE(sp, test="Simes", n=30)

topSpliceDGE(sp, test="gene", n=30)


#par(mfrow=c(1,2))
plotSpliceDGE(sp, geneid="ENSRNOG00000022804", genecol="GeneID")
plotSpliceDGE(sp, geneid="ENSRNOG00000018200", genecol="GeneID")


Gene.Spliced=topSpliceDGE(sp,test="gene",n=nrow(y))

length(rownames(Gene.Spliced)[which(Gene.Spliced$FDR < 0.05)])

write.csv(Gene.Spliced,"Gene.Spliced.LRT.E.N.csv")

Simes.Spliced=topSpliceDGE(sp,test="Simes",n=nrow(y))

length(rownames(Simes.Spliced)[which(Simes.Spliced$FDR < 0.05)])

write.csv(Simes.Spliced,"Simes.Spliced.LRT.E.N.csv")





## LIMMA ##################################
#####################################

N_SE=readRDS("N_SE.rds")

dge <- DGEList(counts=N_SE$counts, genes=N_SE$annotation)



##where counts is a matrix of exon-level counts, and GeneID identifies which gene each exon belongs
##to. Then filter and normalize:

A <- rowSums(dge$counts)
dge <- dge[A>10,, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

##Then apply the voom transformation and fit a linear model:
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
##Now we can test for differential splicing associated with any coefficient in the linear model. First
##run the diffSplice function:
ex <- diffSplice(fit, geneid="GeneID")
#Then
topSplice(ex, coef=2, test="simes")
#will find genes that show evidence of differential splicing associated with the second coefficient in
#the linear model. The output is similar that from the limma topTable function. More detail can be
#obtained by
topSplice(ex, coef=2, test="t")
#which will show individual exons that are enriched or depleted relative to other exons in the same
#gene. To display the pattern of exons in the top genes:
plotSplice(ex)

Simes.Spliced= topSplice(ex, coef=2, test="simes",n=nrow(ex))

length(rownames(Simes.Spliced)[which(Simes.Spliced$FDR < 0.05)])

write.csv(Simes.Spliced,"Simes.Spliced2.E.N.LIMMA.csv")


t.Spliced= topSplice(ex, coef=2, test="t",n=nrow(ex))

length(rownames(t.Spliced)[which(t.Spliced$FDR < 0.05)])

write.csv(t.Spliced,"t.Spliced2.E.N.LIMMA.csv")



###################################################################################################################
#################

N_SE=readRDS("N_SE.rds")

y.all <- DGEList(counts=N_SE$counts, genes=N_SE$annotation)
dim(y.all)
head(y.all$genes)


Groups2=c(rep(1,14),rep(0,16))

y=y.all

y$samples

Groups.Keep <- factor(Groups2)
keep <- filterByExpr(y, group=Groups.Keep)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples

plotMDS(y, col=c(rep(c(1),7),rep(c(2),7),rep(c(3),8),rep(c(4),8)))


##legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)


design <- model.matrix(~ Groups.Keep)
design


y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)

fit <- glmFit(y, design, robust=TRUE)

lrt <- glmLRT(fit, coef=2)
topTags(lrt)

is.de <- decideTests(lrt, p.value=0.05)
summary(is.de)

plotMD(lrt)

abline(h=c(-1, 1), col="blue")

sp <- diffSpliceDGE(fit, coef=2, geneid="GeneID", exonid="Start")

#The top spliced genes under the Simes’ method are shown below
topSpliceDGE(sp, test="Simes", n=30)

topSpliceDGE(sp, test="gene", n=30)


#par(mfrow=c(1,2))
plotSpliceDGE(sp, geneid="ENSRNOG00000022804", genecol="GeneID")
plotSpliceDGE(sp, geneid="ENSRNOG00000018200", genecol="GeneID")


Gene.Spliced=topSpliceDGE(sp,test="gene",n=nrow(y))

length(rownames(Gene.Spliced)[which(Gene.Spliced$FDR < 0.05)])

write.csv(Gene.Spliced,"Gene.Spliced.LRT.E.N.csv")

Simes.Spliced=topSpliceDGE(sp,test="Simes",n=nrow(y))

length(rownames(Simes.Spliced)[which(Simes.Spliced$FDR < 0.05)])

write.csv(Simes.Spliced,"Simes.Spliced.LRT.E.N.csv")

UPPPP+SPLICE

ENSRNOG00000050071=CDC45
ENSRNOG00000038572=Ncapg
ENSRNOG00000003081=Ribc1
ENSRNOG00000003657=Pkmyt1
ENSRNOG00000001439=Srrm3
ENSRNOG00000014287=Stk11

SRRM3 is a neuronal regulator of alternative splicing required
for motor coordination
SRRM3-SRRM4 switch the GABA effect from excitatory to
inhibitory by upregulating Kcc2= ENSRNOG00000018111==DOWN-REGULATED


LRT
ENSRNOG00000020906=Pola2
ENSRNOG00000008829=Sorbs3



DOWN+Spliced

ENSRNOG00000018778
ENSRNOG00000002544
ENSRNOG00000056246
ENSRNOG00000018992
ENSRNOG00000031997
ENSRNOG00000014104
Gene stable ID	Gene name	Human gene name	Human gene stable ID
ENSRNOG00000002544	Kifap3	KIFAP3	ENSG00000075945
ENSRNOG00000014104	Myo5b	MYO5B	ENSG00000167306
ENSRNOG00000018778	Cadm1	CADM1	ENSG00000182985
ENSRNOG00000018992	Dpysl3	DPYSL3	ENSG00000113657
ENSRNOG00000031997	Ttll7	TTLL7	ENSG00000137941
ENSRNOG00000056246	Gls	GLS

LRT

ENSRNOG00000000774
ENSRNOG00000016516
ENSRNOG00000006997
ENSRNOG00000016204
ENSRNOG00000056038
ENSRNOG00000000774	Gabbr1	GABBR1	ENSG00000204681
ENSRNOG00000006997	App	APP	ENSG00000142192
ENSRNOG00000016204	Lrrc24	LRRC24	ENSG00000254402
ENSRNOG00000016516	Mbp	MBP	ENSG00000197971
ENSRNOG00000056038	Ehbp1l1	EHBP1L1	ENSG00000173442

BOTH

ENSRNOG00000006617
ENSRNOG00000013583
ENSRNOG00000010838
ENSRNOG00000054058
ENSRNOG00000018220
ENSRNOG00000018184
ENSRNOG00000000525
ENSRNOG00000010589
ENSRNOG00000018602
ENSRNOG00000013202
ENSRNOG00000001557
ENSRNOG00000012579
ENSRNOG00000020989
ENSRNOG00000005900

ENSRNOG00000000525	Pi16	PI16	ENSG00000164530
ENSRNOG00000001557	Cxadr
ENSRNOG00000005900	Slc2a6	SLC2A6	ENSG00000160326
ENSRNOG00000006617	AABR07060487.1
ENSRNOG00000010589	Zfp622	ZNF622	ENSG00000173545
ENSRNOG00000010838	Araf	ARAF	ENSG00000078061
ENSRNOG00000012579	Tlcd1	TLCD1	ENSG00000160606
ENSRNOG00000013202	Akap7	AKAP7	ENSG00000118507
ENSRNOG00000013583	Tbc1d8	TBC1D8	ENSG00000204634
ENSRNOG00000018184	Tpm1	TPM1	ENSG00000140416
ENSRNOG00000018220	Pde4dip	PDE4DIP	ENSG00000178104
ENSRNOG00000018602	Camta1	CAMTA1	ENSG00000171735
ENSRNOG00000020989	Tm7sf2	TM7SF2	ENSG00000149809
ENSRNOG00000054058	Osbpl1a	OSBPL1A	ENSG00000141447

AKAP7
In their study, Jones et al. determined that knockout of the mouse AKAP7 gene results in redistribution of PKA-RIIβ out of MF axons and terminals, loss of cAMP-induced MF-LTP and specific deficits in behavior including non-cued spatial and contextual pattern recognition (Jones et al., 2016). These findings are particularly interesting as this is the first functional report of a role for pre-synaptic PKA scaffolding by an AKAP.
##The calmodulin-binding transcription activator CAMTA1 is required for long-term memory formation in mice
Foxp4

ENSRNOG00000000387	Slc25a16	DNA2	ENSG00000138346
ENSRNOG00000000839	Nfkbil1	NFKBIL1	ENSG00000204498
ENSRNOG00000002412	Mapk7	MAPK7	ENSG00000166484
ENSRNOG00000002898	Nme7	NME7	ENSG00000143156
ENSRNOG00000005673	Runx1t1	RUNX1T1	ENSG00000079102
ENSRNOG00000007445	Asph
ENSRNOG00000011879	Nfat5	NFAT5	ENSG00000102908
ENSRNOG00000012644	Nup58	NUP58	ENSG00000139496
ENSRNOG00000013465	Tepp	TEPP	ENSG00000159648
ENSRNOG00000013884	Psd3	PSD3	ENSG00000156011
ENSRNOG00000014382	C2cd5	C2CD5	ENSG00000111731
ENSRNOG00000014436	Ppip5k1	PPIP5K1	ENSG00000168781
ENSRNOG00000014811	Cyhr1	CYHR1	ENSG00000187954
ENSRNOG00000016539	Mxd3	RAB24	ENSG00000169228
ENSRNOG00000017441	Tpm3	TPM3	ENSG00000143549
ENSRNOG00000017905	Map1lc3b	MAP1LC3B2	ENSG00000258102
ENSRNOG00000018200	Gad2	GAD2	ENSG00000136750
ENSRNOG00000018350	Shtn1	SHTN1	ENSG00000187164
ENSRNOG00000018863	Abcc10	ABCC10	ENSG00000124574
ENSRNOG00000019550	Slc11a2	SLC11A2	ENSG00000110911
ENSRNOG00000019799	Pcdhgc3	PCDHGC3	ENSG00000240184
ENSRNOG00000020500	Scamp3	CLK2	ENSG00000176444
ENSRNOG00000020578	Ceacam1	CEACAM8	ENSG00000124469
ENSRNOG00000021108	Slc22a12	SLC22A12	ENSG00000197891
ENSRNOG00000022804	Foxp4	FOXP4	ENSG00000137166
ENSRNOG00000023152	Tmem201	TMEM201	ENSG00000188807
ENSRNOG00000024482	Tnrc18	TNRC18	ENSG00000182095
ENSRNOG00000024501	Rgs3	RGS3	ENSG00000138835
ENSRNOG00000025000	Cep104	CEP104	ENSG00000116198
ENSRNOG00000025592	Tspan31	TSPAN31	ENSG00000135452
ENSRNOG00000026649	Dnmt3a	DNMT3A	ENSG00000119772
ENSRNOG00000028856	Pknox2	PKNOX2	ENSG00000165495
ENSRNOG00000030101	Traip	TRAIP	ENSG00000183763

#Dnmt3a regulates gene expression in inhibitory neurons by writing all mCH and some mCG,
#GAD2 Alternative Transcripts in the Human Prefrontal
##Elevated Foxp4 Expression Coincides with Neuronal Differentiation and


#####
setwd("/Users/rafel137/Documents/Aus_Samples2")

Gene.Spliced.GeneID=readRDS("Gene.Spliced.GeneID.rds")
library(ggfortify)
library(DESeq2)
library("factoextra")
library(ggfortify)
library(DESeq2)
library(ggfortify)
library(DESeq2)
library("factoextra")

load("RData")
library(gsubfn)
GAERS.val2=GAERS.val[-4]
NEC.val2=NEC.val[-2]

mat_merged2=readRDS("mat_merged.rds")
colnames(mat_merged2)
Rownames=rownames(mat_merged2)
length(Rownames)

NEC=c(NEC.val2,NEC.control,GAERS.val2,GAERS.control)

NEC.MX=as.matrix(mat_merged2[,NEC])

annotation_col = data.frame(
    Group = factor(c(rep(c("E.NEC"),7),rep(c("N.NEC"),8),rep(c("E.GAERS"),7),rep(c("N.GAERS"),8))))
rownames(annotation_col) = NEC
annotation_col


Group = c(E.NEC = "blue", N.NEC = "lightblue",E.GAERS="red",N.GAERS="lightpink")

library(pheatmap)

INT.Genes=c("ENSRNOG00000053122","ENSRNOG00000013847","ENSRNOG00000022804","ENSRNOG00000018200","ENSRNOG00000018602",
"ENSRNOG00000018778")

NEC.MX2=NEC.MX[INT.Genes,]


mat_merged2=readRDS("mat_merged.rds")
colnames(mat_merged2)
mat_merged=mat_merged2[,-c(10,15)]

NEC=mat_merged[INT.Genes,NEC.val2]
NEC.N=mat_merged[INT.Genes,NEC.control]
GAERS=mat_merged[INT.Genes,GAERS.val2]
GAERS.N=mat_merged[INT.Genes,GAERS.control]

A=list(NEC,NEC.N,GAERS,GAERS.N)
names(A)=c("E.NEC","N.NEC","E.GAERS","N.GAERS")

A2=melt(A)

A2$Var1=gsub("ENSRNOG00000053122","Scn1a",A2$Var1)
A2$Var1=gsub("ENSRNOG00000013847","Nova2",A2$Var1)


A2$Var1=gsub("ENSRNOG00000022804","Foxp4",A2$Var1)
A2$Var1=gsub("ENSRNOG00000018200","Gad2",A2$Var1)


A2$Var1=gsub("ENSRNOG00000018602","Camta1",A2$Var1)
A2$Var1=gsub("ENSRNOG00000018778","Cadm1",A2$Var1)




A2$Var2=NULL

Fills <- c("blue","lightblue","red","lightpink")
ggplot(A2, aes(x = Var1, y = log2(value))) + geom_boxplot(aes(fill = L1)) + scale_fill_manual(values = Fills) + theme_bw() +
 scale_x_discrete(name ="Selected.M30genes") + scale_y_discrete(name ="log2(Expression.Values)")


#



DeconResults=A2 %>%
  group_by(Var1) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Var1, values_from = value) %>%
  select(-row)
dim(DeconResults)


Wilcoxon.Results <- lapply(2:7, function(x) pairwise.wilcox.test(DeconResults[[x]], DeconResults$L1))
names(Wilcoxon.Results) <- names(DeconResults)[2:7]
Wilcoxon.Results
Wilcoxon.Results.PVALUE=sapply(Wilcoxon.Results, function(x) {
  pValue <- x$p.value
  n <- outer(rownames(pValue), colnames(pValue), paste, sep='v')
  pValue <- as.vector(pValue)
  names(pValue) <- n
  pValue
})




Wilcoxon.Results.B=melt(Wilcoxon.Results.PVALUE)
Wilcoxon.Results.B$FDR=p.adjust(Wilcoxon.Results.B$value, method = "BH", n = length(Wilcoxon.Results.B$value))

Wilcoxon.Results.B$log_fdr=-log10(Wilcoxon.Results.B$FDR)

Wilcoxon.Results.B$Genes=rownames( Wilcoxon.Results.B)

Wilcoxon.Results.B$Genes=gsub(".NvE","",Wilcoxon.Results.B$Genes)

Wilcoxon.Results.PLOT=ggplot(Wilcoxon.Results.B,aes(x = Genes, y = log_fdr, fill = log_fdr)) + geom_tile(colour = "black") +
labs(x = "Genes", y = "", fill = "-log"[1][0]~"(FDR)") + scale_fill_gradient(low = "#dddddd", high = "#777777", na.value = "white") + geom_text(aes(label = case_when(log_fdr > 1.30103 ~ "**"))) +
theme(legend.key.size = unit(0.5, "cm"), legend.position = "right", plot.margin = margin(t = 0, r = 0.5, b = 0.5, l = 0.5, "cm"))



FigList=list(DECON.PLOT,Wilcoxon.Results.PLOT)

cowplot::plot_grid(plotlist = FigList,
                           ncol = 1,
                           align = "v", axis = "lr",
                           rel_heights = c(1, 0.65))


######################
###############################################################################
#################

library(edgeR)
library(ggfortify)
library(DESeq2)
library("factoextra")
library(ggfortify)
library(DESeq2)
library(ggfortify)
library(DESeq2)
library("factoextra")
load("RData")
library(gsubfn)
GAERS.val2=GAERS.val[-4]
NEC.val2=NEC.val[-2]

NEC=c(GAERS.val2,NEC.val2,GAERS.control,NEC.control)

N_SE=readRDS("N_SE.rds")
N_SE2=N_SE$counts
colnames(N_SE2)=gsub(".bam","",colnames(N_SE2))
N_SE3=N_SE2[,NEC]
colnames(N_SE3)
dim(N_SE3)
y.all <- DGEList(counts=N_SE3, genes=N_SE$annotation)
dim(y.all)
head(y.all$genes)


Groups2=c(rep(1,14),rep(0,16))

y=y.all

y$samples

Groups.Keep <- factor(Groups2)
keep <- filterByExpr(y, group=Groups.Keep)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples
plotMDS(y)
plotMDS(y, dim.plot = c(2,3),col=c(rep(c(1),14),rep(c(2),16)))

plotMDS(y, dim.plot = c(1,2),col=c(rep(c(1),7),rep(c(2),7),rep(c(3),8),rep(c(4),8)))

design <- model.matrix(~ Groups.Keep)
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

#The top spliced genes under the Simes’ method are shown below
topSpliceDGE(sp, test="Simes", n=30)

topSpliceDGE(sp, test="gene", n=30)


#par(mfrow=c(1,2))
plotSpliceDGE(sp, geneid="ENSRNOG00000022804", genecol="GeneID")
plotSpliceDGE(sp, geneid="ENSRNOG00000018200", genecol="GeneID")


Gene.Spliced=topSpliceDGE(sp,test="gene",n=nrow(y))

length(rownames(Gene.Spliced)[which(Gene.Spliced$FDR < 0.05)])

write.csv(Gene.Spliced,"Gene.Spliced.QLF.All.E.All.N.csv.csv")

Simes.Spliced=topSpliceDGE(sp,test="Simes",n=nrow(y))

length(rownames(Simes.Spliced)[which(Simes.Spliced$FDR < 0.05)])

write.csv(Simes.Spliced,"Simes.Spliced2.All.E.All.N.csv")


###############################################################################

###############################################################################
#################

library(edgeR)
library(ggfortify)
library(DESeq2)
library("factoextra")
library(ggfortify)
library(DESeq2)
library(ggfortify)
library(DESeq2)
library("factoextra")
load("RData")
library(gsubfn)
GAERS.val2=GAERS.val[-4]
NEC.val2=NEC.val[-2]

NEC=c(NEC.val2,NEC.control)

N_SE=readRDS("N_SE.rds")
N_SE2=N_SE$counts
colnames(N_SE2)=gsub(".bam","",colnames(N_SE2))
N_SE3=N_SE2[,NEC]
colnames(N_SE3)
dim(N_SE3)
y.all <- DGEList(counts=N_SE3, genes=N_SE$annotation)
dim(y.all)
head(y.all$genes)


Groups2=c(rep(1,7),rep(0,8))

y=y.all

y$samples

Groups.Keep <- factor(Groups2)
keep <- filterByExpr(y, group=Groups.Keep)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples
plotMDS(y)
plotMDS(y, col=c(rep(c(1),7),rep(c(2),8)))


design <- model.matrix(~ Groups.Keep)
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

#The top spliced genes under the Simes’ method are shown below
topSpliceDGE(sp, test="Simes", n=30)

topSpliceDGE(sp, test="gene", n=30)


#par(mfrow=c(1,2))
plotSpliceDGE(sp, geneid="ENSRNOG00000022804", genecol="GeneID")
plotSpliceDGE(sp, geneid="ENSRNOG00000018200", genecol="GeneID")


Gene.Spliced=topSpliceDGE(sp,test="gene",n=nrow(y))

length(rownames(Gene.Spliced)[which(Gene.Spliced$FDR < 0.05)])

write.csv(Gene.Spliced,"Gene.Spliced.QLF.E.NEC.N.NEC.csv")

Simes.Spliced=topSpliceDGE(sp,test="Simes",n=nrow(y))

length(rownames(Simes.Spliced)[which(Simes.Spliced$FDR < 0.05)])

write.csv(Simes.Spliced,"Simes.Spliced2.E.NEC.N.NEC.csv")


###############################################################################
#################

library(edgeR)
library(ggfortify)
library(DESeq2)
library("factoextra")
library(ggfortify)
library(DESeq2)
library(ggfortify)
library(DESeq2)
library("factoextra")
load("RData")
library(gsubfn)
GAERS.val2=GAERS.val[-4]
NEC.val2=NEC.val[-2]

NEC=c(GAERS.val2,GAERS.control)

N_SE=readRDS("N_SE.rds")
N_SE2=N_SE$counts
colnames(N_SE2)=gsub(".bam","",colnames(N_SE2))
N_SE3=N_SE2[,NEC]
colnames(N_SE3)
dim(N_SE3)
y.all <- DGEList(counts=N_SE3, genes=N_SE$annotation)
dim(y.all)
head(y.all$genes)


Groups2=c(rep(1,7),rep(0,8))

y=y.all

y$samples

Groups.Keep <- factor(Groups2)
keep <- filterByExpr(y, group=Groups.Keep)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples
plotMDS(y)
plotMDS(y, col=c(rep(c(1),7),rep(c(2),8)))


design <- model.matrix(~ Groups.Keep)
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

#The top spliced genes under the Simes’ method are shown below
topSpliceDGE(sp, test="Simes", n=30)

topSpliceDGE(sp, test="gene", n=30)


#par(mfrow=c(1,2))
plotSpliceDGE(sp, geneid="ENSRNOG00000022804", genecol="GeneID")
plotSpliceDGE(sp, geneid="ENSRNOG00000018200", genecol="GeneID")


Gene.Spliced=topSpliceDGE(sp,test="gene",n=nrow(y))

length(rownames(Gene.Spliced)[which(Gene.Spliced$FDR < 0.05)])

write.csv(Gene.Spliced,"Gene.Spliced.QLF.E.GAERS.N.GAERS.csv")

Simes.Spliced=topSpliceDGE(sp,test="Simes",n=nrow(y))

length(rownames(Simes.Spliced)[which(Simes.Spliced$FDR < 0.05)])

write.csv(Simes.Spliced,"Simes.Spliced2.E.GAERS.N.GAERS.csv")

#################

library(edgeR)
library(ggfortify)
library(DESeq2)
library("factoextra")
library(ggfortify)
library(DESeq2)
library(ggfortify)
library(DESeq2)
library("factoextra")
load("RData")
library(gsubfn)
GAERS.val2=GAERS.val[-4]
NEC.val2=NEC.val[-2]

NEC=c(GAERS.control,NEC.control)

N_SE=readRDS("N_SE.rds")
N_SE2=N_SE$counts
colnames(N_SE2)=gsub(".bam","",colnames(N_SE2))
N_SE3=N_SE2[,NEC]
colnames(N_SE3)
dim(N_SE3)
y.all <- DGEList(counts=N_SE3, genes=N_SE$annotation)
dim(y.all)
head(y.all$genes)


Groups2=c(rep(1,8),rep(0,8))

y=y.all

y$samples

Groups.Keep <- factor(Groups2)
keep <- filterByExpr(y, group=Groups.Keep)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples
plotMDS(y)
plotMDS(y, col=c(rep(c(1),8),rep(c(2),8)))


design <- model.matrix(~ Groups.Keep)
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

#The top spliced genes under the Simes’ method are shown below
topSpliceDGE(sp, test="Simes", n=30)

topSpliceDGE(sp, test="gene", n=30)


#par(mfrow=c(1,2))
plotSpliceDGE(sp, geneid="ENSRNOG00000022804", genecol="GeneID")
plotSpliceDGE(sp, geneid="ENSRNOG00000018200", genecol="GeneID")


Gene.Spliced=topSpliceDGE(sp,test="gene",n=nrow(y))

length(rownames(Gene.Spliced)[which(Gene.Spliced$FDR < 0.05)])

write.csv(Gene.Spliced,"Gene.Spliced.QLF.N.GAERS.N.NEC.csv")

Simes.Spliced=topSpliceDGE(sp,test="Simes",n=nrow(y))

length(rownames(Simes.Spliced)[which(Simes.Spliced$FDR < 0.05)])

write.csv(Simes.Spliced,"Simes.Spliced2.N.GAERS.N.NEC.csv")

#################

library(edgeR)
library(ggfortify)
library(DESeq2)
library("factoextra")
library(ggfortify)
library(DESeq2)
library(ggfortify)
library(DESeq2)
library("factoextra")
load("RData")
library(gsubfn)
GAERS.val2=GAERS.val[-4]
NEC.val2=NEC.val[-2]

NEC=c(GAERS.val2,NEC.val2)

N_SE=readRDS("N_SE.rds")
N_SE2=N_SE$counts
colnames(N_SE2)=gsub(".bam","",colnames(N_SE2))
N_SE3=N_SE2[,NEC]
colnames(N_SE3)
dim(N_SE3)
y.all <- DGEList(counts=N_SE3, genes=N_SE$annotation)
dim(y.all)
head(y.all$genes)


Groups2=c(rep(1,7),rep(0,7))

y=y.all

y$samples

Groups.Keep <- factor(Groups2)
keep <- filterByExpr(y, group=Groups.Keep)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples
plotMDS(y)
plotMDS(y, col=c(rep(c(1),7),rep(c(2),7)))


design <- model.matrix(~ Groups.Keep)
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

#The top spliced genes under the Simes’ method are shown below
topSpliceDGE(sp, test="Simes", n=30)

topSpliceDGE(sp, test="gene", n=30)


#par(mfrow=c(1,2))
plotSpliceDGE(sp, geneid="ENSRNOG00000022804", genecol="GeneID")
plotSpliceDGE(sp, geneid="ENSRNOG00000018200", genecol="GeneID")


Gene.Spliced=topSpliceDGE(sp,test="gene",n=nrow(y))

length(rownames(Gene.Spliced)[which(Gene.Spliced$FDR < 0.05)])

write.csv(Gene.Spliced,"Gene.Spliced.QLF.E.GAERS.E.NEC.csv")

Simes.Spliced=topSpliceDGE(sp,test="Simes",n=nrow(y))

length(rownames(Simes.Spliced)[which(Simes.Spliced$FDR < 0.05)])

write.csv(Simes.Spliced,"Simes.Spliced2.E.GAERS.E.NEC.csv")


##########################################################################################################################################################################
##########################################################################################################################################################################
