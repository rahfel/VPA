setwd("/Users/rafel137/Documents/Aus_Samples_2022")

library(ggfortify)
library(DESeq2)
library("factoextra")
library(ggfortify)
library(DESeq2)
library(ggfortify)
library(DESeq2)
library("factoextra")
library(reshape2)
library(dplyr)
rm(list=ls())

##### Main Figure 2 ##########

rm(list=ls())


tmp = list.files(pattern="TB2022.*Gender.csv")
names(tmp)=tmp
Allfiles = lapply(tmp,read.csv,sep=",",row.names=1)
names(Allfiles)=gsub(".Gender.csv","",names(Allfiles))
names(Allfiles)
names(Allfiles)=gsub("TB2022.","",names(Allfiles))
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

upset(ALL2,text.scale= c(1.5),sets=c("N.GAERS.vs.N.NEC","E.GAERS.vs.E.NEC","All.E.vs.All.N","E.GAERS.vs.N.GAERS","E.NEC.vs.N.NEC"),sets.bar.color =c("#000000", "red4", "#cccccc","salmon","#8c8c8c"),set_size.show=F,keep.order = TRUE)

upset(ALL2,text.scale= c(1.5),sets=c("N.GAERS.vs.N.NEC","E.GAERS.vs.E.NEC","All.E.vs.All.N","E.GAERS.vs.N.GAERS","E.NEC.vs.N.NEC")
,sets.bar.color =c("#000000", "red4", "#cccccc","salmon","#8c8c8c"),set_size.show=F,keep.order = TRUE,show.numbers = "no", nintersects = 10)



str(AllRegulated)


library(ComplexHeatmap)


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

########################################

##### Main Figure 3 ##########



setwd("/Users/rafel137/Documents/Aus_Samples_2022")
library(dplyr)
library(ggplot2)
library(reshape2)
library(gplots)
library(tidyverse)



UP.BP=read.csv("UpRegulated.0.9.BP.2022.csv",sep=",")

DOWN.BP=read.csv("DownRegulated.0.9.BP.2022.csv",sep=",")


UP.BP=na.omit(UP.BP[rownames(UP.BP)[which(UP.BP$overlap > 10)],])

UP.BP$gene_list=gsub("TB2021_E_GAERS_vs_N_GAERS","E-GAERS.vs.N-GAERS",UP.BP$gene_list)
UP.BP$gene_list=gsub("TB2021_N_GAERS_vs_N_NEC","N-GAERS.vs.N-NEC",UP.BP$gene_list)


UP.BP2=UP.BP %>% group_by(parent_term) %>% summarise(fdr = min(fdr)) %>% mutate(parent_term = fct_reorder(parent_term, -fdr))
UP.BP2$Log10=-log10(as.numeric(UP.BP2$fdr))

UP.BP2$Log10[which(!is.finite(UP.BP2$Log10))] <- max(UP.BP2$Log10[is.finite(UP.BP2$Log10)])


DOWN.BP=na.omit(DOWN.BP[rownames(DOWN.BP)[which(DOWN.BP$overlap > 10)],])

DOWN.BP$gene_list=gsub("TB2021_E_GAERS_vs_N_GAERS","E-GAERS.vs.N-GAERS",DOWN.BP$gene_list)
DOWN.BP$gene_list=gsub("TB2021_N_GAERS_vs_N_NEC","N-GAERS.vs.N-NEC",DOWN.BP$gene_list)


DOWN.BP2=DOWN.BP %>% group_by(parent_term) %>% summarise(fdr = min(fdr)) %>% mutate(parent_term = fct_reorder(parent_term, -fdr))
DOWN.BP2$Log10=-log10(as.numeric(DOWN.BP2$fdr))

DOWN.BP2$Log10[which(!is.finite(DOWN.BP2$Log10))] <- max(DOWN.BP2$Log10[is.finite(DOWN.BP2$Log10)])
DOWN.BP2$Log10=-1*(DOWN.BP2$Log10)
par(mar = c(5, 5, 5, 5))

DOWN.BP3 <- Reduce(rbind,by(DOWN.BP2, DOWN.BP2["Log10"], head, n=15))[1:15,]


A=rbind(UP.BP2,DOWN.BP3)
A2 <- A[order(A$parent_term, -A$Log10),]

par(mar = c(5, 23, 4, 4))
barplot(height=A2$Log10, names=A2$parent_term,horiz=T,las = 1,cex.names=0.7,xlim=c(-15,10),xlab="-Log10(FDR)",main="GO: Biological processes",
col =ifelse(A2$Log10>0,"salmon","lightblue"))



#########

##### Main Figure 4 ##########


library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
rm(list=ls())
setwd("/Users/rafel137/Documents/Aus_Samples_2022/LDSC_2022/")
A=read.csv("All.Final.TB_FDR_LDSC2022.csv",sep=",",header=T,row.names="file")
A$X=NULL
A$Comparison=rownames(A)
A2=melt(A)
A2

A3=A2 %>% separate(variable, c("GWAS.SummaryStat","DEGeneSets"), sep = "([.?:])")
A3$DEGeneSets=gsub("ALL","All.Dysregulated",A3$DEGeneSets)
A3$DEGeneSets=gsub("DOWN","DownRegulated",A3$DEGeneSets)
A3$DEGeneSets=gsub("UP","UpRegulated",A3$DEGeneSets)
table(A3$GWAS.SummaryStat)


ggplot(A3,aes(x = Comparison, y = -log10(value), fill = Comparison)) +
    geom_col(alpha = 0.8, width = 0.85) +
    scale_y_continuous(expand = c(0, 0.5)) +
    coord_flip() +
    facet_grid(ordered(GWAS.SummaryStat, levels = c("ADHD", "ASD", "BD","SCZ","IQ","EPI","CDG","WHR")) ~ DEGeneSets) +
    labs(y = "-log10(FDR)") + geom_hline(yintercept=-log10(0.05),size=0.25) +
    theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 4),
    legend.position = "none",
    panel.grid.major.y = element_blank(),) + theme_bw() + scale_fill_manual(values=c("All.E.vs.All.N"="black","E.NEC.vs.N.NEC"="red4","N.GAERS.vs.N.NEC"="#cccccc","E.GAERS.vs.E.NEC"="salmon","E.GAERS.vs.N.GAERS"="#8c8c8c"))



##### Supplementary Figure 1 ################

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


Wilcoxon.Results.PLOT=ggplot(Wilcoxon.Results.B2,aes(x = Celltypes, y = fct_rev(Comparison), fill = log_fdr)) + geom_tile(colour = "black") +
labs(x = "Celltype", y = "", fill = "-log"[1][0]~"(FDR)") + scale_fill_gradient(low = "white", high = "grey", na.value = "white") + geom_text(aes(label = case_when(log_fdr > 1.5 ~ "*" , log_fdr > 2 ~ "**"))) +
theme(legend.key.size = unit(0.5, "cm"), legend.position = "right", plot.margin = margin(t = 0, r = 0.5, b = 0.5, l = 0.5, "cm"))


FigList=list(DECON.PLOT,Wilcoxon.Results.PLOT)

cowplot::plot_grid(plotlist = FigList,
                           ncol = 1,
                           align = "v", axis = "lr",
                           rel_heights = c(1, 0.35))


ggplot(A, aes(x = Celltypes, y = Estimated.Celltype.Fraction)) + geom_boxplot(aes(fill = Group)) + scale_fill_manual(values = Fills) + theme_bw()




##### Supplementary Figure 2 ################

load("RData")
colnames(mat_merged)

mat_merged=mat_merged[,-c(10,15)]
colnames(mat_merged)

rlogcounts <- rlog(mat_merged)



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
             palette = c("#DC3977","#045275","#F0746E","#089099"),
             addEllipses = FALSE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = F,label="none"
)
