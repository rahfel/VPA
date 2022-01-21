
#########################################################

library(ggfortify)
library(DESeq2)
library("factoextra")
library(ggfortify)
library(DESeq2)
setwd("/Users/rafel137/Documents/Aus_Samples2")
library(ggfortify)
library(DESeq2)
library("factoextra")
mat_merged2=readRDS("mat_merged.rds")
colnames(mat_merged2)

mat_merged=mat_merged2[,-c(10,15)]
colnames(mat_merged)

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
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
             )


#####################

## GAERS.val vs GAERS.control
samples.id=c(GAERS.val2,GAERS.control)
data.tmp2=mat_merged[,samples.id]

data.tmp=data.tmp2[ rowSums(data.tmp2 > 0) >= 1, ]
dim(data.tmp)
group=c(rep(1,length(GAERS.val2)),rep(0,length(GAERS.control)))
design=model.matrix(~group)

y <- DGEList(counts=data.tmp,group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef="group")
res.GAERS.val.GAERS.control=topTags(lrt,n=nrow(data.tmp))
AB.GAERS.val.GAERS.control=res.GAERS.val.GAERS.control$table
length(rownames(AB.GAERS.val.GAERS.control)[which(AB.GAERS.val.GAERS.control$FDR < 0.05)])
write.csv(AB.GAERS.val.GAERS.control,"TB2021.E.GAERS.vs.N.GAERS.csv")

#########################################################################################

====NEC.val vs NEC.control
samples.id=c(NEC.val2,NEC.control)
data.tmp2=mat_merged[,samples.id]

data.tmp=data.tmp2[ rowSums(data.tmp2 > 0) >= 1, ]
dim(data.tmp)
group=c(rep(1,length(NEC.val2)),rep(0,length(NEC.control)))
design=model.matrix(~group)

y <- DGEList(counts=data.tmp,group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef="group")
res.NEC.val.NEC.control=topTags(lrt,n=nrow(data.tmp))
AB.NEC.val.NEC.control=res.NEC.val.NEC.control$table
length(rownames(AB.NEC.val.NEC.control)[which(AB.NEC.val.NEC.control$FDR < 0.05)])
write.csv(AB.NEC.val.NEC.control,"TB2021.E.NEC.vs.N.NEC.csv")

#########################################################################################

# GAERS vs NEC
samples.id=c(GAERS.control,NEC.control)
data.tmp2=mat_merged[,samples.id]
data.tmp=data.tmp2[ rowSums(data.tmp2 > 0) >= 1, ]
dim(data.tmp)
group=c(rep(1,length(GAERS.control)),rep(0,length(NEC.control)))
design=model.matrix(~group)
y <- DGEList(counts=data.tmp,group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef="group")
res.GAERS.control.NEC.control=topTags(lrt,n=nrow(data.tmp))
AB.GAERS.control.NEC.control=res.GAERS.control.NEC.control$table
length(rownames(AB.GAERS.control.NEC.control)[which(AB.GAERS.control.NEC.control$FDR < 0.05)])
write.csv(AB.GAERS.control.NEC.control,"TB2021.N.GAERS.vs.N.NEC.csv")


#########################################################################################

# val
samples.id=c(GAERS.val2,NEC.val2)
data.tmp2=mat_merged[,samples.id]
data.tmp=data.tmp2[ rowSums(data.tmp2 > 0) >= 1, ]
dim(data.tmp)
group=c(rep(1,length(GAERS.val2)),rep(0,length(NEC.val2)))
design=model.matrix(~group)
y <- DGEList(counts=data.tmp,group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef="group")
res.GAERS.val.NEC.val=topTags(lrt,n=nrow(data.tmp))
AB.GAERS.val.NEC.val=res.GAERS.val.NEC.val$table
length(rownames(AB.GAERS.val.NEC.val)[which(AB.GAERS.val.NEC.val$FDR < 0.05)])
write.csv(AB.GAERS.val.NEC.val,"TB2021.E.GAERS.vs.E.NEC.csv")

#########################################################################################
#########################################################################################

# val
samples.id=c(GAERS.val2,NEC.control)
data.tmp2=mat_merged[,samples.id]
data.tmp=data.tmp2[ rowSums(data.tmp2 > 0) >= 1, ]
dim(data.tmp)
group=c(rep(1,length(GAERS.val2)),rep(0,length(NEC.control)))
design=model.matrix(~group)
y <- DGEList(counts=data.tmp,group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef="group")
res.GAERS.val.NEC.control=topTags(lrt,n=nrow(data.tmp))
AB.GAERS.val.NEC.control=res.GAERS.val.NEC.control$table
length(rownames(AB.GAERS.val.NEC.control)[which(AB.GAERS.val.NEC.control$FDR < 0.05)])
write.csv(AB.GAERS.val.NEC.control,"TB2021.E.GAERS.vs.N.NEC.csv")

#########################################################################################


samples.id=c(GAERS.val2,NEC.val2,GAERS.control,NEC.control)
data.tmp2=mat_merged[,samples.id]
data.tmp=data.tmp2[ rowSums(data.tmp2 > 0) >= 1, ]
dim(data.tmp)

group=c(rep(1,14),rep(0,16))
design=model.matrix(~group)

y <- DGEList(counts=data.tmp,group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef="group")
res.Val.Control=topTags(lrt,n=nrow(data.tmp))
AB.Val.Control=res.Val.Control$table
length(rownames(AB.Val.Control)[which(AB.Val.Control$FDR < 0.05)])
write.csv(AB.Val.Control,"TB2021.All.E.vs.All.N.csv")


#########################################################################################

########### change  TB2021 to TB2022 ####################################################### FINAL #######

library(UpSetR)
library(ggplot2)
library(reshape2)
library(dplyr)
setwd("/Users/rafel137/Documents/Aus_Samples2/")
rm(list=ls())
library(WebGestaltR)
tmp = list.files(pattern="TB2022.*.csv")
names(tmp)=tmp
Allfiles = lapply(tmp,read.csv,sep=",",row.names=1)
names(Allfiles)=gsub(".csv","",names(Allfiles))
names(Allfiles)

Function1=function(data,x){rownames(data)[which(data[[x]] < 0.05)]}
AllRegulated=lapply(Allfiles,Function1,x="FDR")
str(AllRegulated)

Function2=function(data,x,y){rownames(data)[which(data[[y]] > 0 & data[[x]] < 0.05)]}
UpRegulated=lapply(Allfiles, Function2,x="FDR",y="logFC")
str(UpRegulated)


Function3=function(data,x,y){rownames(data)[which(data[[y]] < 0 & data[[x]] < 0.05)]}
DownRegulated=lapply(Allfiles, Function3,x="FDR", y="logFC")
str(DownRegulated)



UpRegulated=list(UpRegulated)
DownRegulated=list(DownRegulated)
library(dplyr)
library(reshape2)
UP=melt(UpRegulated)
DOWN=melt(DownRegulated)
length(unique(UP$value))
length(unique(DOWN$value))

head(UP)
UP2=dcast(data = UP,formula = value~L2,value.var = "L1")
UP2[is.na(UP2)] <- 0
upset(UP2,mb.ratio = c(0.3, 0.7),text.scale= c(1.5))

head(DOWN)
DOWN2=dcast(data = DOWN,formula = value~L2,value.var = "L1")
DOWN2[is.na(DOWN2)] <- 0
upset(DOWN2,mb.ratio = c(0.3, 0.7),text.scale= c(1.5))


AllRegulated1=list(AllRegulated)
ALL=melt(AllRegulated1)
head(ALL)
ALL$L2=gsub("TB2022.","",ALL$L2)
table(ALL$L2)

ALL2=dcast(data = ALL,formula = value~L2,value.var = "L1")
ALL2[is.na(ALL2)] <- 0

upset(ALL2, sets = c("All.E.vs.All.N","E.NEC.vs.N.NEC","E.GAERS.vs.N.GAERS","E.GAERS.vs.E.NEC","N.GAERS.vs.N.NEC"), mb.ratio = c(0.55, 0.45),order.by = "degree", keep.order = TRUE)

upset(ALL2,text.scale= c(1.5),sets=c("N.GAERS.vs.N.NEC","E.GAERS.vs.E.NEC","All.E.vs.All.N","E.GAERS.vs.N.GAERS","E.NEC.vs.N.NEC"),sets.bar.color =c("#000000", "red4", "#cccccc","salmon","#8c8c8c"),set_size.show=F,keep.order = TRUE)

upset(ALL2,text.scale= c(1.5),sets=c("N.GAERS.vs.N.NEC","E.GAERS.vs.E.NEC","All.E.vs.All.N","E.GAERS.vs.N.GAERS","E.NEC.vs.N.NEC")
,sets.bar.color =c("#000000", "red4", "#cccccc","salmon","#8c8c8c"),set_size.show=F,keep.order = TRUE,show.numbers = "no", nintersects = 10)

ggplot(mtcars, aes(x = wt, y = qsec)) + geom_line(color = "blue", size = 2) + geom_point(color = "red", size = 3)



str(AllRegulated)


library(ComplexHeatmap)

tiff("/Users/rafel137/Documents/Aus_Samples2/Main2.tiff", units="in", width=5, height=3.5, res=300)

lt = list(All.E.vs.All.N = AllRegulated[[1]],
          E.NEC.vs.N.NEC = AllRegulated[[4]],
          E.GAERS.vs.N.GAERS = AllRegulated[[3]])

m = make_comb_mat(lt,mode = "intersect")
cs = comb_size(m)
UpSet(m)

UpSet(m, pt_size = unit(5, "mm"), lwd = 4,
      comb_col = c("brown", "black", "red")[comb_degree(m)],column_names_gp = gpar(fontsize = 15),comb_order = order(comb_degree(m), -cs),
      set_order=c("E.NEC.vs.N.NEC","E.GAERS.vs.N.GAERS","All.E.vs.All.N"),
      row_names_gp = gpar(fontsize = 10),top_annotation = HeatmapAnnotation(
          "Number of Genes" = anno_barplot(
              comb_size(m),
              border = FALSE,
              gp = gpar(fill = c("brown", "black", "red")[comb_degree(m)],fontsize = 16),
              height = unit(5, "cm"),
              axis_param = list(side = "left"),
              right_annotation = upset_right_annotation(m),
              add_numbers = TRUE,axis=T)))

###
dev.off()

###

comb_col = c("brown", "black", "red")[comb_degree(m)],column_names_gp = gpar(fontsize = 5),comb_order = order(comb_degree(m), -cs)
set_order=c("E.NEC.vs.N.NEC","E.GAERS.vs.N.GAERS","All.E.vs.All.N")
row_names_gp = gpar(fontsize = 10)


####


lt = list(All.E.vs.All.N = AllRegulated[[1]],
          E.NEC.vs.N.NEC = AllRegulated[[4]],
          E.GAERS.vs.N.GAERS = AllRegulated[[3]],
          E.GAERS.vs.E.NEC=AllRegulated[[2]],
          N.GAERS.vs.N.NEC=AllRegulated[[5]])
m = make_comb_mat(lt,mode = "intersect",top_n_sets = 5)
cs = comb_size(m)

tiff("/Users/rafel137/Documents/Aus_Samples2/Supp3.tiff", units="in", width=7, height=3.5, res=300)

UpSet(m, pt_size = unit(2, "mm"), lwd = 3,
    comb_col = c("brown", "black", "red","salmon","darkgrey")[comb_degree(m)],column_names_gp = gpar(fontsize = 15),comb_order = order(comb_degree(m), -cs),
set_order=c("E.NEC.vs.N.NEC","E.GAERS.vs.N.GAERS","All.E.vs.All.N","N.GAERS.vs.N.NEC","E.GAERS.vs.E.NEC"),
    row_names_gp = gpar(fontsize = 10),top_annotation = HeatmapAnnotation(
        "Number of Genes" = anno_barplot(
            comb_size(m),
            border = FALSE,
            gp = gpar(fill = c("brown", "black", "red")[comb_degree(m)],fontsize = 16),
            height = unit(5, "cm"),
            axis_param = list(side = "left"),
            right_annotation = upset_right_annotation(m),
            add_numbers = TRUE,axis=T)))

dev.off()




#######

library(dplyr)
library(tidyr)
rm(list=ls())
setwd("/Users/rafel137/Documents/Aus_Samples/LDSC_HMAGMA_RESULTS_2021/")
A=read.csv("All.Final.TB_FDR_LDSC2021.csv",sep=",",header=T,row.names="file")
A$X=NULL
A$Comparison=rownames(A)
A$Comparison=gsub("All.VPA.Exposed.VS.All.Non.Exposed","All-E.vs.All-N",A$Comparison)
A$Comparison=gsub("Non.Exposed.GAERS.vs.Non.Exposed.NEC","N-GAERS.vs.N-NEC",A$Comparison)
A$Comparison=gsub("VPA.Exposed.GAERS.vs.Non.Exposed.GAERS","E-GAERS.vs.N-GAERS",A$Comparison)
A$Comparison=gsub("VPA.Exposed.GAERS.vs.VPA.Exposed.NEC","E-GAERS.vs.E-NEC",A$Comparison)
A$Comparison=gsub("VPA.Exposed.NEC.vs.Non.Exposed.NEC","E-NEC.vs.N-NEC",A$Comparison)
A2=melt(A)
A2
A2$variable=gsub("IQSavageJansen2018.ALL","IQ",A2$variable)


A3=A2 %>% separate(variable, c("GWAS.SummaryStat","DEGeneSets"), sep = "([.?:])")
A3$DEGeneSets=gsub("ALL","All.Dysregulated",A3$DEGeneSets)
A3$DEGeneSets=gsub("DOWN","DownRegulated",A3$DEGeneSets)
A3$DEGeneSets=gsub("UP","UpRegulated",A3$DEGeneSets)
table(A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("ADHD2017","ADHD",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("ASD2017","ASD",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("BIP2019","BD",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("CDG2019","CDG",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("EPI2018","EPI",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("SCZ2020","SCZ",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("FATSL2019","WHR",A3$GWAS.SummaryStat)

A3$GWAS.SummaryStat=gsub("IQ2018","IQ",A3$GWAS.SummaryStat)
tiff("/Users/rafel137/Documents/Aus_Samples2/Main4.tiff", units="in", width=10, height=7, res=300)

P=ggplot(A3,aes(x = Comparison, y = -log10(value), fill = Comparison)) +
    geom_col(alpha = 0.8, width = 0.85) +
    scale_y_continuous(expand = c(0, 0.5)) +
    coord_flip() +
    facet_grid(ordered(GWAS.SummaryStat, levels = c("ADHD", "ASD", "BD","SCZ","IQ","EPI","CDG","WHR")) ~ DEGeneSets) +
    labs(y = "-log10(FDR)") + geom_hline(yintercept=-log10(0.05),size=0.25) +
    theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
    plot.title = element_text(size = 15, face = "bold"),
    strip.text.y = element_text(angle = 270, face = "bold"),
    strip.placement = "outside",
    axis.title.x = element_text(margin = margin(t = 1, b = 1, unit = "cm")),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 4),
    legend.position = "none",
    panel.grid.major.y = element_blank(),) + theme_bw()
P+ scale_fill_manual(values=c("All-E.vs.All-N"="black","E-NEC.vs.N-NEC"="red4","N-GAERS.vs.N-NEC"="#cccccc","E-GAERS.vs.E-NEC"="salmon","E-GAERS.vs.N-GAERS"="#8c8c8c"))
ordered(GWAS.SummaryStat, levels = c("ADHD", "ASD", "BD","SCZ","IQ","EPI","CDG","WHR"))
dev.off()

#######


################################################################################




setwd("/Users/rafel137/Documents/Aus_Samples2")
library(dplyr)
library(ggplot2)
library(reshape2)
library(gplots)
library(tidyverse)



UP.BP=read.csv("UpRegulated.0.9.BP.csv",sep=",")
UP.CC=read.csv("UpRegulated.0.9.CC.csv",sep=",")
UP.MF=read.csv("UpRegulated.0.9.MF.csv",sep=",")

DOWN.BP=read.csv("DownRegulated.0.9.BP.csv",sep=",")
DOWN.CC=read.csv("DownRegulated.0.9.CC.csv",sep=",")
DOWN.MF=read.csv("DownRegulated.0.9.MF.csv",sep=",")

ALL.BP=read.csv("AllRegulated.0.9.BP.csv",sep=",")
ALL.CC=read.csv("AllRegulated.0.9.CC.csv",sep=",")
ALL.MF=read.csv("AllRegulated.0.9.MF.csv",sep=",")



UP=Reduce(rbind,list(UP.BP,UP.CC,UP.MF))
UP2=UP[rownames(UP)[which(UP$overlap > 10)],]
UP3 <- na.omit(UP2)

UP3$gene_list=gsub("TB2021_E_GAERS_vs_N_GAERS","E-GAERS.vs.N-GAERS",UP3$gene_list)
UP3$gene_list=gsub("TB2021_N_GAERS_vs_N_NEC","N-GAERS.vs.N-NEC",UP3$gene_list)


UP3.BP.All.E=UP3[UP3$gene_list == "All-E.vs.All-N" & UP3$go_type == "BP" ,]


UP3.BP.All.E2=UP3.BP.All.E %>% group_by(parent_term) %>% summarise(fdr = min(fdr)) %>% mutate(parent_term = fct_reorder(parent_term, -fdr))
UP3.BP.All.E2$Log10=-log10(as.numeric(UP3.BP.All.E2$fdr))

UP3.BP.All.E2$Log10[which(!is.finite(UP3.BP.All.E2$Log10))] <- max(UP3.BP.All.E2$Log10[is.finite(UP3.BP.All.E2$Log10)])


DOWN=Reduce(rbind,list(DOWN.BP,DOWN.MF))
DOWN2=DOWN[rownames(DOWN)[which(DOWN$overlap > 10)],]
DOWN3 <- na.omit(DOWN2)

DOWN3$gene_list=gsub("TB2021_E_GAERS_vs_N_GAERS","E-GAERS.vs.N-GAERS",DOWN3$gene_list)
DOWN3$gene_list=gsub("TB2021_N_GAERS_vs_N_NEC","N-GAERS.vs.N-NEC",DOWN3$gene_list)

DOWN3.BP.All.E=DOWN3[DOWN3$gene_list == "All-E.vs.All-N" & DOWN3$go_type == "BP" ,]


DOWN3.BP.All.E2=DOWN3.BP.All.E %>% group_by(parent_term) %>% summarise(fdr = min(fdr)) %>% mutate(parent_term = fct_reorder(parent_term, -fdr))
DOWN3.BP.All.E2$Log10=-log10(as.numeric(DOWN3.BP.All.E2$fdr))
DOWN3.BP.All.E2$Log10[which(!is.finite(DOWN3.BP.All.E2$Log10))] <- max(DOWN3.BP.All.E2$Log10[is.finite(DOWN3.BP.All.E2$Log10)])
DOWN3.BP.All.E2$Log10=-1*(DOWN3.BP.All.E2$Log10)
par(mar = c(5, 5, 5, 5))
#barplot(height=sort(DOWN3.BP.All.E2$Log10), names=DOWN3.BP.All.E2$parent_term,col="grey",horiz=T, las=1,xlim = c(0, 10))

DOWN3.BP.All.E3 <- Reduce(rbind,by(DOWN3.BP.All.E2, DOWN3.BP.All.E2["Log10"], head, n=15))[1:15,]


A=rbind(UP3.BP.All.E2,DOWN3.BP.All.E3)
A2 <- A[order(A$parent_term, -A$Log10),]

tiff("/Users/rafel137/Documents/Aus_Samples2/Main2.tiff", units="in", width=7, height=6, res=300)

par(mar = c(5, 23, 4, 4))
barplot(height=A2$Log10, names=A2$parent_term,horiz=T,las = 1,cex.names=0.7,xlim=c(-15,10),xlab="-Log10(FDR)",main="GO Parental terms",
col =ifelse(A2$Log10>0,"salmon","lightblue"))

dev.off()


#########################################################

###


library(ggfortify)
library(DESeq2)
library("factoextra")
library(ggfortify)
library(DESeq2)
setwd("/Users/rafel137/Documents/Aus_Samples")
library(ggfortify)
library(DESeq2)
library("factoextra")
mat_merged=readRDS("mat_merged.rds")
mat_merged2=mat_merged[,-c(10,15)]
colnames(mat_merged2)
rlogcounts <- rlog(mat_merged2)
library(reshape2)
load("RData")
A=list(GAERS.val,GAERS.control,NEC.val,NEC.control)
names(A)=c("E-GAERS","N-GAERS","E-NEC","N-NEC")
groups=melt(A)
colnames(groups)=c("SampleID","Groups")
groups2 = groups[with(groups, order(SampleID)), ]
rownames(groups2) = 1:nrow(groups2)
groups2
groups3=groups2[-c(10,15),]
groups2B=as.factor(groups3$Groups)

# run PCA
pcDat <- prcomp(t(rlogcounts))
# plot PCA
autoplot(pcDat)
fviz_pca_ind(pcDat)

tiff("/Users/rafel137/Documents/Aus_Samples2/Supp2.tiff", units="in", width=7, height=4, res=300)


fviz_pca_ind(pcDat,
             col.ind = groups2B, # color by groups
             palette = c("#DC3977","#045275","#F0746E","#089099"),
             addEllipses = FALSE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE,label="none"
)

dev.off()
###
############################\

setwd("/Users/rafel137/Documents/Aus_Samples/DECON_MOUSE/DWLS.results")

rm(list=ls())
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)



N.NEC=melt(read.csv("CellFractions_PFCref_DWLS.N.NEC.MOUSEEMB18.csv",sep=",",header=1,row.names=1))
N.NEC$Group=c(rep("N.NEC",nrow(N.NEC)))
colnames(N.NEC)=c("Celltypes","Estimated.Celltype.Fraction","Group")


E.NEC=melt(read.csv("CellFractions_PFCref_DWLS.E.NEC.MOUSEEMB18.csv",sep=",",header=1,row.names=1))
E.NEC$Group=c(rep("E.NEC",nrow(E.NEC)))
colnames(E.NEC)=c("Celltypes","Estimated.Celltype.Fraction","Group")

E.GAERS=melt(read.csv("CellFractions_PFCref_DWLS.E.GAERS.MOUSEEMB18.csv",sep=",",header=1,row.names=1))
E.GAERS$Group=c(rep("E.GAERS",nrow(E.GAERS)))
colnames(E.GAERS)=c("Celltypes","Estimated.Celltype.Fraction","Group")

N.GAERS=melt(read.csv("CellFractions_PFCref_DWLS.N.GAERS.MOUSEEMB18.csv",sep=",",header=1,row.names=1))
N.GAERS$Group=c(rep("N.GAERS",nrow(N.GAERS)))
colnames(N.GAERS)=c("Celltypes","Estimated.Celltype.Fraction","Group")


A=Reduce(rbind,list(N.NEC,N.GAERS,E.NEC,E.GAERS))


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

#Wilcoxon.Results.B=melt(Wilcoxon.Results.PVALUE[Wilcoxon.Results.names,])

Wilcoxon.Results.B=melt(Wilcoxon.Results.PVALUE)
Wilcoxon.Results.B$FDR=p.adjust(Wilcoxon.Results.B$value, method = "BH", n = length(Wilcoxon.Results.B$value))

colnames(Wilcoxon.Results.B)=c("Comparison","Celltypes","pvalue","FDR")
Wilcoxon.Results.B$log_fdr=-log10(Wilcoxon.Results.B$FDR)

Rownames=c(rownames(Wilcoxon.Results.B)[which(Wilcoxon.Results.B$Comparison == "N.NECvN.GAERS")],rownames(Wilcoxon.Results.B)[which(Wilcoxon.Results.B$Comparison == "N.NECvE.NEC")],
rownames(Wilcoxon.Results.B)[which(Wilcoxon.Results.B$Comparison == "N.NECvE.GAERS")])

Wilcoxon.Results.B2=Wilcoxon.Results.B[intersect(rownames(Wilcoxon.Results.B),Rownames),]


Wilcoxon.Results.PLOT=ggplot(Wilcoxon.Results.B2,aes(x = Celltypes, y = fct_rev(Comparison), fill = log_fdr)) + geom_tile(colour = "black") +
labs(x = "Celltype", y = "", fill = "-log"[1][0]~"(FDR)") + scale_fill_gradient(low = "white", high = "grey", na.value = "white") + geom_text(aes(label = case_when(log_fdr > 1.5 ~ "*" , log_fdr > 2 ~ "**"))) +
theme(legend.key.size = unit(0.5, "cm"), legend.position = "right", plot.margin = margin(t = 0, r = 0.5, b = 0.5, l = 0.5, "cm"))


FigList=list(DECON.PLOT,Wilcoxon.Results.PLOT)

cowplot::plot_grid(plotlist = FigList,
                           ncol = 1,
                           align = "v", axis = "lr",
                           rel_heights = c(1, 0.35))

tiff("/Users/rafel137/Documents/Aus_Samples2/Supp1.tiff", units="in", width=9.25, height=4, res=300)


cowplot::plot_grid(plotlist = FigList,
                           ncol = 1,
                           align = "v", axis = "lr",
                           rel_heights = c(1, 0.35))
dev.off()


tiff("/Users/rafel137/Documents/Aus_Samples2/Supp1b.tiff", units="in", width=9.25, height=4, res=300)

ggplot(A, aes(x = Celltypes, y = Estimated.Celltype.Fraction)) + geom_boxplot(aes(fill = Group)) + scale_fill_manual(values = Fills) + theme_bw()

dev.off()
##############################################################################
