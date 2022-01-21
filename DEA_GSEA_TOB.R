
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

tiff("test.tiff", units="in", width=5, height=3.5, res=300)

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

comb_col = c("brown", "black", "red")[comb_degree(m)],column_names_gp = gpar(fontsize = 5),comb_order = order(comb_degree(m), -cs),
set_order=c("E.NEC.vs.N.NEC","E.GAERS.vs.N.GAERS","All.E.vs.All.N"),
row_names_gp = gpar(fontsize = 10)


####


lt = list(All.E.vs.All.N = AllRegulated[[1]],
          E.NEC.vs.N.NEC = AllRegulated[[4]],
          E.GAERS.vs.N.GAERS = AllRegulated[[3]],
          E.GAERS.vs.E.NEC=AllRegulated[[2]],
          N.GAERS.vs.N.NEC=AllRegulated[[5]])
m = make_comb_mat(lt,mode = "intersect",top_n_sets = 5)
cs = comb_size(m)

tiff("test2.tiff", units="in", width=7, height=3.5, res=300)

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









upset(ALL2,
      sets = c("E.NEC.vs.N.NEC", "E.GAERS.vs.N.GAERS", "All.E.vs.All.N"),
      order.by="degree", matrix.color="black", point.size=5,
      sets.bar.color=c("#8c8c8c","salmon","#cccccc"),show.numbers = "no")

UpSet(ALL2, pt_size = unit(5, "mm"), lwd = 3,
    comb_col = c("red", "blue", "black")[comb_degree(m)])


xxxx
library('UpSetR')
movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"),
                    header=T, sep=";" )
require(ggplot2); require(plyr); require(gridExtra); require(grid);
## Loading required package: ggplot2
## Loading required package: plyr
## Loading required package: gridExtra
## Loading required package: grid
upset(movies,
      sets = c("Action", "Comedy", "Drama"),
      order.by="degree", matrix.color="blue", point.size=5,
      sets.bar.color=c("maroon","blue","orange"))




xxxx
## Initial inputs on the first example
movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"),
                    header=T, sep=";" )
## comma -> semicolon
data = movies; sets = c("Action", "Comedy", "Drama");
      order.by="degree"; matrix.color="blue"; point.size=5;
      sets.bar.color=c("maroon","blue","orange")



xxxxxx



## Piece of code we introduced
for(i in 1:3) {
      j <- which(Matrix_layout$y == i & Matrix_layout$value == 1)
      if(length(j) > 0) Matrix_layout$color[j] <- c("maroon","blue","orange")[i]
  }

xxxx

## Modified internal upset() code

startend <- UpSetR:::FindStartEnd(data)
  first.col <- startend[1]
  last.col <- startend[2]

  if(color.pal == 1){
    palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2",
                 "#7F7F7F", "#BCBD22", "#17BECF")
  } else{
    palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                 "#CC79A7")
  }

  if(is.null(intersections) == F){
    Set_names <- unique((unlist(intersections)))
    Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
    New_data <- UpSetR:::Wanted(data, Sets_to_remove)
    Num_of_set <- UpSetR:::Number_of_sets(Set_names)
    if(keep.order == F){
      Set_names <- UpSetR:::order_sets(New_data, Set_names)
    }
    All_Freqs <- UpSetR:::specific_intersections(data, first.col, last.col, intersections, order.by, group.by, decreasing,
                                        cutoff, main.bar.color, Set_names)
  } else if(is.null(intersections) == T){
    Set_names <- sets
    if(is.null(Set_names) == T || length(Set_names) == 0 ){
      Set_names <- UpSetR:::FindMostFreq(data, first.col, last.col, nsets)
    }
    Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
    New_data <- UpSetR:::Wanted(data, Sets_to_remove)
    Num_of_set <- UpSetR:::Number_of_sets(Set_names)
    if(keep.order == F){
    Set_names <- UpSetR:::order_sets(New_data, Set_names)
    }
    All_Freqs <- UpSetR:::Counter(New_data, Num_of_set, first.col, Set_names, nintersects, main.bar.color,
                         order.by, group.by, cutoff, empty.intersections, decreasing)
  }
  Matrix_setup <- UpSetR:::Create_matrix(All_Freqs)
  labels <- UpSetR:::Make_labels(Matrix_setup)
  #Chose NA to represent NULL case as result of NA being inserted when at least one contained both x and y
  #i.e. if one custom plot had both x and y, and others had only x, the y's for the other plots were NA
  #if I decided to make the NULL case (all x and no y, or vice versa), there would have been alot more if/else statements
  #NA can be indexed so that we still get the non NA y aesthetics on correct plot. NULL cant be indexed.
  att.x <- c(); att.y <- c();
  if(is.null(attribute.plots) == F){
    for(i in seq_along(attribute.plots$plots)){
      if(length(attribute.plots$plots[[i]]$x) != 0){
        att.x[i] <- attribute.plots$plots[[i]]$x
      }
      else if(length(attribute.plots$plots[[i]]$x) == 0){
        att.x[i] <- NA
      }
      if(length(attribute.plots$plots[[i]]$y) != 0){
        att.y[i] <- attribute.plots$plots[[i]]$y
      }
      else if(length(attribute.plots$plots[[i]]$y) == 0){
        att.y[i] <- NA
      }
    }
  }

  BoxPlots <- NULL
  if(is.null(boxplot.summary) == F){
    BoxData <- UpSetR:::IntersectionBoxPlot(All_Freqs, New_data, first.col, Set_names)
    BoxPlots <- list()
    for(i in seq_along(boxplot.summary)){
      BoxPlots[[i]] <- UpSetR:::BoxPlotsPlot(BoxData, boxplot.summary[i], att.color)
    }
  }

  customAttDat <- NULL
  customQBar <- NULL
  Intersection <- NULL
  Element <- NULL
  legend <- NULL
  EBar_data <- NULL
  if(is.null(queries) == F){
    custom.queries <- UpSetR:::SeperateQueries(queries, 2, palette)
    customDat <- UpSetR:::customQueries(New_data, custom.queries, Set_names)
    legend <- UpSetR:::GuideGenerator(queries, palette)
    legend <- UpSetR:::Make_legend(legend)
    if(is.null(att.x) == F && is.null(customDat) == F){
      customAttDat <- UpSetR:::CustomAttData(customDat, Set_names)
    }
    customQBar <- UpSetR:::customQueriesBar(customDat, Set_names, All_Freqs, custom.queries)
  }
  if(is.null(queries) == F){
    Intersection <- UpSetR:::SeperateQueries(queries, 1, palette)
    Matrix_col <- UpSetR:::intersects(QuerieInterData, Intersection, New_data, first.col, Num_of_set,
                             All_Freqs, expression, Set_names, palette)
    Element <- UpSetR:::SeperateQueries(queries, 1, palette)
    EBar_data <-UpSetR:::ElemBarDat(Element, New_data, first.col, expression, Set_names,palette, All_Freqs)
  } else{
    Matrix_col <- NULL
  }

  Matrix_layout <- UpSetR:::Create_layout(Matrix_setup, matrix.color, Matrix_col, matrix.dot.alpha)



xxxxx


Matrix_layout
##    y x value  color alpha Intersection
## 1  1 1     1   blue   1.0         1yes
## 2  2 1     1   blue   1.0         1yes
## 3  3 1     1   blue   1.0         1yes
## 4  1 2     0 gray83   0.5          4No
## 5  2 2     1   blue   1.0         2yes
## 6  3 2     1   blue   1.0         2yes
## 7  1 3     1   blue   1.0         3yes
## 8  2 3     0 gray83   0.5          8No
## 9  3 3     1   blue   1.0         3yes
## 10 1 4     1   blue   1.0         4yes
## 11 2 4     1   blue   1.0         4yes
## 12 3 4     0 gray83   0.5         12No
## 13 1 5     0 gray83   0.5         13No
## 14 2 5     0 gray83   0.5         14No
## 15 3 5     1   blue   1.0         5yes
## 16 1 6     0 gray83   0.5         16No
## 17 2 6     1   blue   1.0         6yes
## 18 3 6     0 gray83   0.5         18No
## 19 1 7     1   blue   1.0         7yes
## 20 2 7     0 gray83   0.5         20No
## 21 3 7     0 gray83   0.5         21No







for (i in 1:length(AllRegulated)) {
  write.table(AllRegulated[i], file=paste0(names(AllRegulated)[i], ".AllRegulated.txt"),sep="\t",row.names=F)
}

for (i in 1:length(UpRegulated)) {
  write.table(UpRegulated[i], file=paste0(names(UpRegulated)[i], ".UpRegulated.txt"),sep="\t",row.names=F)
}


for (i in 1:length(DownRegulated)) {
  write.table(DownRegulated[i], file=paste0(names(DownRegulated)[i], ".DownRegulated.txt"),sep="\t",row.names=F)
}

Function4=function(data){rownames(data)}
AllGenes=lapply(Allfiles, Function4)
str(AllGenes)


for (i in 1:length(AllGenes)) {
  write.table(AllGenes[i], file=paste0(names(AllGenes)[i], ".AllRegulated.AllGenes.txt"),sep="\t",row.names=F)
}
for (i in 1:length(AllGenes)) {
  write.table(AllGenes[i], file=paste0(names(AllGenes)[i], ".UpRegulated.AllGenes.txt"),sep="\t",row.names=F)
}
for (i in 1:length(AllGenes)) {
  write.table(AllGenes[i], file=paste0(names(AllGenes)[i], ".DownRegulated.AllGenes.txt"),sep="\t",row.names=F)
}
str(AllGenes)


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
ALL2=dcast(data = ALL,formula = value~L2,value.var = "L1")
ALL2[is.na(ALL2)] <- 0
upset(ALL2,mb.ratio = c(0.3, 0.7),text.scale= c(1.5),sets.bar.color =c("#000000", "red4", "#cccccc","salmon","#8c8c8c"))

upset(ALL2,text.scale= c(1.5),sets.bar.color =c("#000000", "red4", "#cccccc","salmon","#8c8c8c"),set_size.show=F
)




rm(list=ls())
tmp = list.files(pattern="TB2021.*.txt")
Allfiles = lapply(tmp,read.delim,sep="\t")
filenames=list.files(pattern="TB2021.*.txt", full.names=TRUE)
filenames
names(Allfiles) = gsub(".*/(.*)\\..*", "\\1", filenames)
names(Allfiles)
Allfiles2 <- Allfiles[ grepl("*Regulated$", names(Allfiles)) ]
names(Allfiles2)


Inpath="/Users/rafel137/Documents/Aus_Samples2/"

for (i in 1:length(Allfiles2)){
WebGestaltR(enrichMethod="ORA", organism="rnorvegicus",
 enrichDatabase="geneontology_Biological_Process", interestGeneFile=paste0(Inpath,names(Allfiles2)[i],".txt"),
 interestGeneType="ensembl_gene_id", referenceGeneFile=paste0(Inpath,names(Allfiles2)[i],".AllGenes.txt"),
referenceGeneType="ensembl_gene_id", isOutput=TRUE,
 outputDirectory=Inpath, projectName=paste0("BP.",names(Allfiles2)[i]))}


 for (i in 1:length(Allfiles2)){
 WebGestaltR(enrichMethod="ORA", organism="rnorvegicus",
  enrichDatabase="geneontology_Molecular_Function", interestGeneFile=paste0(Inpath,names(Allfiles2)[i],".txt"),
  interestGeneType="ensembl_gene_id", referenceGeneFile=paste0(Inpath,names(Allfiles2)[i],".AllGenes.txt"),
 referenceGeneType="ensembl_gene_id", isOutput=TRUE,
  outputDirectory=Inpath, projectName=paste0("MF.",names(Allfiles2)[i]))}

  for (i in 1:length(Allfiles2)){
  WebGestaltR(enrichMethod="ORA", organism="rnorvegicus",
   enrichDatabase="geneontology_Cellular_Component", interestGeneFile=paste0(Inpath,names(Allfiles2)[i],".txt"),
   interestGeneType="ensembl_gene_id", referenceGeneFile=paste0(Inpath,names(Allfiles2)[i],".AllGenes.txt"),
  referenceGeneType="ensembl_gene_id", isOutput=TRUE,
   outputDirectory=Inpath, projectName=paste0("CC.",names(Allfiles2)[i]))}

##############################


rm(list=ls())
library(gplots)
library(reshape2)
library(ggplot2)
library(rutils)
dirs = list.dirs()
dirs = grep("_BP_TB2021_", dirs, value = TRUE)
tmp = list.files(path=dirs,pattern = "^enrichment_results_(.*)txt$",full.names = TRUE)
Allfiles = lapply(tmp,read.delim,sep="\t")
str(Allfiles)
filenames=list.files(path=dirs,pattern = "^enrichment_results_(.*)txt$",full.names = TRUE)
names(Allfiles) = gsub(".*/(.*)\\..*", "\\1", filenames)
names(Allfiles)
Allfiles2=lapply(Allfiles,function (x) {x <- x[,c("FDR","geneSet","description")]})
names(Allfiles2) = gsub("enrichment_results_BP_", "", names(Allfiles2))
names(Allfiles2)
ALL <- Allfiles2[ grepl("*_AllRegulated$", names(Allfiles2)) ]
str(ALL)
ALL2=melt(ALL)
ALL2$L1=gsub("_*_AllRegulated","",ALL2$L1)
head(ALL2)
ALL2$variable=gsub("FDR","BP",ALL2$variable)
colnames(ALL2)=c("go_id","go_term","go_type","fdr","gene_list")

Allfiles3=lapply(Allfiles,function (x) {x <- x[,c("overlap","geneSet","description")]})
names(Allfiles3) = gsub("enrichment_results_BP_", "", names(Allfiles3))
names(Allfiles3)
ALL.B <- Allfiles3[ grepl("*_AllRegulated$", names(Allfiles3)) ]
str(ALL.B)
ALL2B=melt(ALL.B)
ALL2B$L1=gsub("_*_AllRegulated","",ALL2B$L1)
head(ALL2B)
colnames(ALL2B)=c("go_id","go_term","variable","overlap","gene_list")
head(ALL2B)

ALL2$overlap=ALL2B$overlap

library(tibble)
ALL3=as_tibble(ALL2)
head(ALL3)
 table(ALL3$gene_list)



   ALL3$gene_list=gsub("TB2021_All_E_vs_All_N","All-E.vs.All-N",ALL3$gene_list)
   ALL3$gene_list=gsub("TB2021_E_GAERS_vs_E_NEC","E-GAERS.vs.E-NEC",ALL3$gene_list)
   ALL3$gene_list=gsub("TB2021_E_GAERS_vs_N_GAERS ","E-GAERS.vs.N-GAERS",ALL3$gene_list)
   ALL3$gene_list=gsub("TB2021_E_NEC_vs_N_NEC","E-NEC.vs.N-NEC",ALL3$gene_list)
   ALL3$gene_list=gsub("TB2021_N_GAERS_vs_N_NEC","N-GAERS.vs.N-NEC",ALL3$gene_list)


unique(ALL3$go_term)


ALL4=go_reduce(
   pathway_df = ALL3,
   threshold = 0.9,
   scores = NULL,
   measure = "Wang")


rowSums(table(ALL4$gene_list,ALL4$parent_term))

summary(ALL4$parent_sim_score)

ALL4 %>% dplyr::count(parent_term)

unique(ALL4$parent_term)

ALL4$Direction=c(rep("ALL",length(nrow(ALL4))))

write.csv(ALL4,"AllRegulated.0.9.BP.csv")


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
UP3.BP.All.E2=UP3.BP.All.E %>% group_by(parent_term) %>% summarise(fdr = min(fdr)) %>% mutate(parent_term = fct_reorder(parent_term, -fdr)) %>%


ggplot(UP3.BP.All.E2,mapping = aes(x = -log10(fdr), y = parent_term)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) + theme_bw()

UP3.BP.All.E2=UP3.BP.All.E %>% group_by(parent_term) %>% summarise(fdr = min(fdr)) %>% mutate(parent_term = fct_reorder(parent_term, -fdr))
UP3.BP.All.E2$Log10=-log10(as.numeric(UP3.BP.All.E2$fdr))

UP3.BP.All.E2$Log10[which(!is.finite(UP3.BP.All.E2$Log10))] <- max(UP3.BP.All.E2$Log10[is.finite(UP3.BP.All.E2$Log10)])
par(mar = c(5, 25, 5, 5))
barplot(height=sort(UP3.BP.All.E2$Log10), names=UP3.BP.All.E2$parent_term,col="grey",horiz=T, las=1,xlim = c(0, 10))




DOWN=Reduce(rbind,list(DOWN.BP,DOWN.MF))
DOWN2=DOWN[rownames(DOWN)[which(DOWN$overlap > 10)],]
DOWN3 <- na.omit(DOWN2)

DOWN3$gene_list=gsub("TB2021_E_GAERS_vs_N_GAERS","E-GAERS.vs.N-GAERS",DOWN3$gene_list)
DOWN3$gene_list=gsub("TB2021_N_GAERS_vs_N_NEC","N-GAERS.vs.N-NEC",DOWN3$gene_list)

DOWN3.BP.All.E=DOWN3[DOWN3$gene_list == "All-E.vs.All-N" & DOWN3$go_type == "BP" ,]
DOWN3.BP.All.E2=DOWN3.BP.All.E %>% group_by(parent_term) %>% summarise(fdr = min(fdr)) %>% mutate(parent_term = fct_reorder(parent_term, -fdr)) %>%


ggplot(DOWN3.BP.All.E2,mapping = aes(x = -log10(fdr), y = parent_term)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) + theme_bw()

DOWN3.BP.All.E2=DOWN3.BP.All.E %>% group_by(parent_term) %>% summarise(fdr = min(fdr)) %>% mutate(parent_term = fct_reorder(parent_term, -fdr))
DOWN3.BP.All.E2$Log10=-log10(as.numeric(DOWN3.BP.All.E2$fdr))

DOWN3.BP.All.E2$Log10[which(!is.finite(DOWN3.BP.All.E2$Log10))] <- max(DOWN3.BP.All.E2$Log10[is.finite(DOWN3.BP.All.E2$Log10)])
DOWN3.BP.All.E2$Log10=-1*(DOWN3.BP.All.E2$Log10)
par(mar = c(5, 25, 5, 5))
barplot(height=sort(DOWN3.BP.All.E2$Log10), names=DOWN3.BP.All.E2$parent_term,col="grey",horiz=T, las=1,xlim = c(0, 10))


barplot(height=sort(A$Log10), names=A$parent_term,col="grey",horiz=T, las=1,xlim = c(-10, 10))




P=ggplot(UP3,aes(x = desc(-log10(fdr)), y = parent_term, fill = gene_list)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    facet_grid( ~ go_type) + theme_bw()

P+scale_fill_manual(values=c("All-E.vs.All-N"="black","E-NEC.vs.N-NEC"="red4","N-GAERS.vs.N-NEC"="#cccccc","E-GAERS.vs.E-NEC"="salmon","E-GAERS.vs.N-GAERS"="#8c8c8c"))



UP3.BP.All.E=UP3[UP3$gene_list == "All-E.vs.All-N" & UP3$go_type == "BP" ,]
UP3.BP.All.E %>% group_by(parent_term) %>% summarise(fdr = min(fdr)) %>%
mutate(parent_term = fct_reorder(parent_term, -fdr)) %>%
ggplot(UP3.BP.All.E,mapping = aes(x = -log10(fdr), y = parent_term)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) + theme_bw()

barplot(height=UP3.BP.All.E$fdr, names=UP3.BP.All.E$parent_term,col="#69b3a2",horiz=T, las=1)




UP3.CC=UP3[UP3$go_type == "CC",]
ggplot(UP3.CC,aes(x = -log10(fdr), y = parent_term, fill = gene_list)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) + theme_bw() + scale_fill_manual(values=c("All-E.vs.All-N"="black","E-NEC.vs.N-NEC"="red4","N-GAERS.vs.N-NEC"="#cccccc","E-GAERS.vs.E-NEC"="salmon","E-GAERS.vs.N-GAERS"="#8c8c8c"))


UP3.MF=UP3[UP3$go_type == "MF",]
ggplot(UP3.MF,aes(x = -log10(fdr), y = parent_term, fill = gene_list)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) + theme_bw() + scale_fill_manual(values=c("All-E.vs.All-N"="black","E-NEC.vs.N-NEC"="red4","N-GAERS.vs.N-NEC"="#cccccc","E-GAERS.vs.E-NEC"="salmon","E-GAERS.vs.N-GAERS"="#8c8c8c"))








ALL=Reduce(rbind,list(ALL.BP,ALL.MF))
ALL2=ALL[rownames(ALL)[which(ALL$overlap > 20)],]
ALL3 <- na.omit(ALL2)

ALL3$gene_list=gsub("TB2021_E_GAERS_vs_N_GAERS","E-GAERS.vs.N-GAERS",ALL3$gene_list)
ALL3$gene_list=gsub("TB2021_N_GAERS_vs_N_NEC","N-GAERS.vs.N-NEC",ALL3$gene_list)



P=ggplot(ALL3,aes(x = -log10(fdr), y = parent_term, fill = gene_list)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    facet_grid( ~ go_type) + theme_bw()

P+scale_fill_manual(values=c("All-E.vs.All-N"="black","E-NEC.vs.N-NEC"="red4","N-GAERS.vs.N-NEC"="#cccccc","E-GAERS.vs.E-NEC"="salmon","E-GAERS.vs.N-GAERS"="#8c8c8c"))



DOWN=Reduce(rbind,list(DOWN.BP,DOWN.MF))
DOWN2=DOWN[rownames(DOWN)[which(DOWN$overlap > 10)],]
DOWN3 <- na.omit(DOWN2)

DOWN3$gene_list=gsub("TB2021_E_GAERS_vs_N_GAERS","E-GAERS.vs.N-GAERS",DOWN3$gene_list)
DOWN3$gene_list=gsub("TB2021_N_GAERS_vs_N_NEC","N-GAERS.vs.N-NEC",DOWN3$gene_list)

DOWN3.BP.All.E=DOWN3[DOWN3$gene_list == "All-E.vs.All-N" & DOWN3$go_type == "BP" ,]
DOWN3.BP.All.E %>% group_by(parent_term) %>% summarise(fdr = min(fdr)) %>%
mutate(parent_term = fct_reorder(parent_term, -fdr)) %>%
ggplot(DOWN3.BP.All.E,mapping = aes(x = -log10(fdr), y = parent_term)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) + theme_bw()


P=ggplot(DOWN3,aes(x = -log10(fdr), y = parent_term, fill = gene_list)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    facet_grid( ~ go_type) + theme_bw()

P+scale_fill_manual(values=c("All-E.vs.All-N"="black","E-NEC.vs.N-NEC"="red4","N-GAERS.vs.N-NEC"="#cccccc","E-GAERS.vs.E-NEC"="salmon","E-GAERS.vs.N-GAERS"="#8c8c8c"))



DOWN3.BP.All.E=DOWN3[DOWN3$gene_list == "All-E.vs.All-N" & DOWN3$go_type == "BP" ,]
DOWN3.BP.All.E %>% group_by(parent_term) %>% summarise(fdr = min(fdr)) %>%
mutate(parent_term = fct_reorder(parent_term, -fdr)) %>%
ggplot(DOWN3.BP.All.E,mapping = aes(x = -log10(fdr), y = parent_term)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) + theme_bw()

barplot(height=DOWN3.BP.All.E$fdr, names=DOWN3.BP.All.E$parent_term,col="#69b3a2",horiz=T, las=1)



#########################################################################################
#########################################################################################





#########################################################################################
#########################################################################################
#########################################################################################

library(ggplot2)
library(reshape2)
library(dplyr)
setwd("/Users/rafel137/Documents/Aus_Samples2/")
rm(list=ls())
library(WebGestaltR)

Rat_6_Human_Description=read.csv("Rat_6_Human_Description.csv",sep=",",header=T)

tmp = list.files(pattern="TB2021.*.csv")
names(tmp)=tmp
Allfiles = lapply(tmp,read.csv,sep=",",row.names=1)
names(Allfiles)=gsub(".csv","",names(Allfiles))
names(Allfiles)

for (i in 1:length(Allfiles)) {
Allfiles[[i]]$Gene.Name= Rat_6_Human_Description$Gene.name[match(paste(rownames(Allfiles[[i]])), paste(Rat_6_Human_Description$Gene.stable.ID))]
Allfiles[[i]]$Gene.description= Rat_6_Human_Description$Gene.description[match(paste(rownames(Allfiles[[i]])), paste(Rat_6_Human_Description$Gene.stable.ID))]
write.csv(Allfiles[[i]], file=paste0(names(Allfiles)[[i]], "2.csv"),sep=",")
}

#########################################################################################

#########################################################################################

library(ggplot2)
library(reshape2)
library(dplyr)
setwd("/Users/rafel137/Documents/Aus_Samples2/")
rm(list=ls())
library(WebGestaltR)

Rat_6_Human_Description=read.csv("Rat_6_Human_Description.csv",sep=",",header=T)

tmp = list.files(pattern="Gene.Spliced.QLF.*.csv")
names(tmp)=tmp
Allfiles = lapply(tmp,read.csv,sep=",",row.names=1)
names(Allfiles)=gsub(".csv","",names(Allfiles))
names(Allfiles)
Allfiles=Allfiles[-c(4,6)]
names(Allfiles)


for (i in 1:length(Allfiles)) {
Allfiles[[i]]$Gene.Name= Rat_6_Human_Description$Gene.name[match(paste(Allfiles[[i]]$GeneID), paste(Rat_6_Human_Description$Gene.stable.ID))]
Allfiles[[i]]$Gene.description= Rat_6_Human_Description$Gene.description[match(paste(Allfiles[[i]]$GeneID), paste(Rat_6_Human_Description$Gene.stable.ID))]
write.csv(Allfiles[[i]], file=paste0(names(Allfiles)[[i]], ".edgeR.csv"),sep=",")
}

#########################################################################################


#########################################################################################
#########################################################################################
################################################################################


setwd("/Users/rafel137/Documents/Aus_Samples/HMAGMA")
rm(list=ls())
tmp = list.files(pattern = "*.Adult_Brain.*.gsa.out",full.names = TRUE)
tmp

Allfiles = lapply(tmp,read.table,sep="",header=T)
str(Allfiles)
filenames=list.files(pattern="*.Adult_Brain.*.gsa.out", full.names=TRUE)
filenames

names(Allfiles) = gsub(".*/(.*)\\..*", "\\1", filenames)
names(Allfiles)

for( i in 1:length(Allfiles)){
    Allfiles[[i]]$FDR <- p.adjust(Allfiles[[i]]$P, method = "BH", n = length(Allfiles[[i]]$P))
}

Allfiles2=lapply(Allfiles,function (x) {x <- x[,c("FULL_NAME","FDR")]})
names(Allfiles2) = names(Allfiles)

Allfiles2a=melt(Allfiles2)
Allfiles2a$FULL_NAME =gsub("TB.","",Allfiles2a$FULL_NAME)


head(Allfiles2a)

library(dplyr)
Allfiles2a$Regulated=Allfiles2a$FULL_NAME
#Allfiles2a2=Allfiles2a %>% separate(L1, c("A", "Regulated2"))  %>% separate(FULL_NAME, c("GeneSet","B"))
 Allfiles2a2= Allfiles2a
head(Allfiles2a2)
Allfiles2a2$L1 =gsub("_ENTID_DEgenes.gsa","",Allfiles2a2$L1)
Allfiles2a2$L1 =gsub("Adult_Brain.TB_","",Allfiles2a2$L1)
Allfiles2a2$L1 =gsub("Adult_Brain.TB_","",Allfiles2a2$L1)
Allfiles2a2$L1 =gsub("_",".",Allfiles2a2$L1)

Allfiles2a2$FULL_NAME=gsub(".AllRegulated","",Allfiles2a2$FULL_NAME)
Allfiles2a2$FULL_NAME=gsub(".UpRegulated","",Allfiles2a2$FULL_NAME)
Allfiles2a2$FULL_NAME=gsub(".DownRegulated","",Allfiles2a2$FULL_NAME)


Allfiles2a2a=dcast(Allfiles2a2,FULL_NAME~L1)
head(Allfiles2a2a)
write.csv(Allfiles2a2a,"TB_FDR_HMAGMA2021.BIP.csv")

#######################


################################################################################ DO NOT USE THIS BECAUSE IT IS FDR corrected for ALL THE TEST CARRIED OUT ###############
setwd("/Users/rafel137/Documents/Aus_Samples/HMAGMA/")
library(dplyr)
library(tidyr)
library(reshape2)
rm(list=ls())
A=read.csv("TB_FDR_HMAGMA2021.BIP.csv",sep=",",header=T,row.names=1)

A$Comparison=A$FULL_NAME

A$Comparison=gsub("All.VPA.Exposed.VS.All.Non.Exposed","All-E.vs.All-N",A$Comparison)
A$Comparison=gsub("Non.Exposed.GAERS.vs.Non.Exposed.NEC","N-GAERS.vs.N-NEC",A$Comparison)
A$Comparison=gsub("VPA.Exposed.GAERS.vs.Non.Exposed.GAERS","E-GAERS.vs.N-GAERS",A$Comparison)
A$Comparison=gsub("VPA.Exposed.GAERS.vs.VPA.Exposed.NEC","E-GAERS.vs.E-NEC",A$Comparison)
A$Comparison=gsub("VPA.Exposed.NEC.vs.Non.Exposed.NEC","E-NEC.vs.N-NEC",A$Comparison)
A2=melt(A)
A2$variable=gsub("all.epilepsy","EPI2018",A2$variable)

A3=A2 %>% separate(variable, c("GWAS.SummaryStat","DEGeneSets"), sep = "([.])")
A3$DEGeneSets=gsub("ALL","All.dysregulated",A3$DEGeneSets)
A3$DEGeneSets=gsub("DOWN","DownRegulated",A3$DEGeneSets)
A3$DEGeneSets=gsub("UP","UpRegulated",A3$DEGeneSets)
head(A3)
A3$GWAS.SummaryStat=gsub("ADHD2017","ADHD",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("ASD2017","ASD",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("BIP2019","BD",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("CDG2019","CDG",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("EPI2018","EPI",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("FATSL2018","WHR",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("IQJansen2018","IQ",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("PGCBI2021","PGC",A3$GWAS.SummaryStat)
A3$GWAS.SummaryStat=gsub("SCZ2020","SCZ",A3$GWAS.SummaryStat)
head(A3)

P=ggplot(A3,aes(x = Comparison, y = -log10(value), fill = Comparison)) +
    geom_col(alpha = 0.8, width = 0.85) +
    scale_y_continuous(expand = c(0, 0.01)) +
    coord_flip() +
    facet_grid(ordered(GWAS.SummaryStat, levels = c("ADHD", "ASD", "BD","SCZ","IQ","EPI","CDG","WHR","PGC")) ~ DEGeneSets) +
    labs(y = "-log10(FDR)") + geom_hline(yintercept=-log10(0.05),size=0.25) +
    theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
    plot.title = element_text(size = 15, face = "bold"),
    strip.text.y = element_text(angle = 270, face = "bold"),
    strip.placement = "outside",
    axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 4),
    legend.position = "none",
    panel.grid.major.y = element_blank(),) + theme_bw()
P+ scale_fill_manual(values=c("All-E.vs.All-N"="black","E-NEC.vs.N-NEC"="red4","N-GAERS.vs.N-NEC"="#cccccc","E-GAERS.vs.E-NEC"="salmon","E-GAERS.vs.N-GAERS"="#8c8c8c"))
ordered(GWAS.SummaryStat, levels = c("ADHD", "ASD", "BD","SCZ","IQ","EPI","CDG","WHR","PGC"))

################################################################################


####################################################

library(data.table)
library(ggplot2)
library(cowplot)
rm(list=ls())
setwd("/Users/rafel137/Documents/Aus_Samples/LDSC_RESULTS2/LDSC_Results_IQ")
FOLDERS=list.files()
count=0
for(FOLD in FOLDERS){
	All.files=list.files(path=FOLD,pattern=".results")
	for(file in All.files){
		data=fread(sprintf("%s/%s",FOLD,file))
		tmp=data[1,]
		tmp$FOLD=FOLD
		tmp$file=file
		if(count==0){
			Results.Table=tmp
		}else{
			Results.Table=rbind(Results.Table,tmp)
		}
		count=count+1
	}
}
Results.Table$ct=gsub("_\\w[^_]*.results","",Results.Table$file,perl=TRUE)
Results.Table$ct=gsub(".results","",Results.Table$file,perl=TRUE)
Results.Table$ct=gsub(paste(FOLDERS,collapse="|"),"",Results.Table$ct)
Results.Table$ct=gsub("_$","",Results.Table$ct)
ggplot(Results.Table) + geom_bar(aes(x=ct,y=Enrichment_p,fill=FOLD),stat="identity",position="dodge") + scale_y_log10() + coord_flip() + theme_cowplot() + geom_hline(yintercept=0.05/dim(Results.Table)[1])

Results.Table$FDR= p.adjust(Results.Table$Enrichment_p, method = "BH", n = length(Results.Table$Enrichment_p))

Results.Table$file=gsub("CDG.","",Results.Table$file)
Results.Table$file=gsub(".results","",Results.Table$file)
Results.Table$file=gsub("TB.","",Results.Table$file)

Results.Table=as.data.frame(Results.Table)
rownames(Results.Table)=Results.Table$file
head(Results.Table)

write.table(Results.Table,"PGCBIP2021.DOWN.Results.Table.txt",sep="\t",row.names=T,quote=FALSE)

##########################################################################################################

#######################


library(data.table)
library(ggplot2)
library(cowplot)
rm(list=ls())
setwd("/Users/rafel137/Documents/Aus_Samples/LDSC_ALL")
tmp = list.files(pattern = "*.txt",full.names = TRUE)

Allfiles = lapply(tmp,read.table,sep="\t")
str(Allfiles)
filenames=list.files(pattern="*.txt", full.names=TRUE)
names(Allfiles) = gsub(".*/(.*)\\..*", "\\1", filenames)
names(Allfiles)



Allfiles2=lapply(Allfiles,function (x) {x <- x[,c("file","FDR")]})
names(Allfiles2) = names(Allfiles)

Allfiles2a=melt(Allfiles2)
head(Allfiles2a)
Allfiles2a$file =gsub(".AllRegulated.results","",Allfiles2a$file)
Allfiles2a$file =gsub(".AllRegulateds","",Allfiles2a$file)
Allfiles2a$file =gsub(".AllRegulated","",Allfiles2a$file)
Allfiles2a$L1 =gsub(".Results.Table","",Allfiles2a$L1)
head(Allfiles2a)

Allfiles2a$file=gsub("TB.","", Allfiles2a$file)
FINAL=dcast(Allfiles2a,file~L1)
dim(FINAL)
head(FINAL)
write.csv(FINAL,"TB_FDR_ALL_LDSC2021.BIP.csv")

#########################################################################################
#########################################################################################

######################################################################



library(data.table)
library(ggplot2)
library(cowplot)
rm(list=ls())
setwd("/Users/rafel137/Documents/Aus_Samples/LDSC_DOWN")
tmp = list.files(pattern = "*.txt",full.names = TRUE)

Allfiles = lapply(tmp,read.table,sep="\t")
str(Allfiles)
filenames=list.files(pattern="*.txt", full.names=TRUE)
names(Allfiles) = gsub(".*/(.*)\\..*", "\\1", filenames)
names(Allfiles)



Allfiles2=lapply(Allfiles,function (x) {x <- x[,c("file","FDR")]})
names(Allfiles2) = names(Allfiles)

Allfiles2a=melt(Allfiles2)
head(Allfiles2a)
Allfiles2a$file =gsub(".DownRegulated.results","",Allfiles2a$file)
Allfiles2a$file =gsub(".DownRegulateds","",Allfiles2a$file)
Allfiles2a$file =gsub(".DownRegulated","",Allfiles2a$file)
Allfiles2a$L1 =gsub(".Results.Table","",Allfiles2a$L1)
head(Allfiles2a)

Allfiles2a$file=gsub("TB.","", Allfiles2a$file)
FINAL=dcast(Allfiles2a,file~L1)
dim(FINAL)
head(FINAL)
write.csv(FINAL,"TB_FDR_DOWN_LDSC2021.BIP.csv")


############ MOVE ALL CSV FILES TO  setwd("/Users/rafel137/Documents/Aus_Samples/LDSC_HMAGMA_RESULTS_2021")


setwd("/Users/rafel137/Documents/Aus_Samples/LDSC_HMAGMA_RESULTS_2021/")

tmp = list.files(pattern = "*.BIP.csv",full.names = TRUE)

Allfiles = lapply(tmp,read.csv,sep=",",row.names=1)
str(Allfiles)
filenames=list.files(pattern="*.BIP.csv", full.names=TRUE)
names(Allfiles) = gsub(".*/(.*)\\..*", "\\1", filenames)
names(Allfiles)
Allfiles2=melt(Allfiles)

Allfiles2a =dcast(Allfiles2,file~variable)
Allfiles2a

write.csv(Allfiles2a,"All.Final.TB_FDR_LDSC2021.BIP.csv")



############ MOVE ALL CSV FILES TO  setwd("/Users/rafel137/Documents/Aus_Samples/LDSC_HMAGMA_RESULTS_2021")

#### Open the two files and make sure the names are the same TB_FDR_HMAGMA2021.csv and All.Final.TB_FDR_LDSC2021.csv
################################################################################

#########################################################################################
########################################################################################
#########################################################################################
#########################################################################################

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

dev.copy2pdf(file="LDSC.ALL.pdf")

################################################################################

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
labs(x = "Celltype", y = "", fill = "-log"[1][0]~"(FDR)") + scale_fill_gradient(low = "#dddddd", high = "#777777", na.value = "white") + geom_text(aes(label = case_when(log_fdr > 1 ~ "*" , log_fdr > 1.30103 ~ "**"))) +
theme(legend.key.size = unit(0.5, "cm"), legend.position = "right", plot.margin = margin(t = 0, r = 0.5, b = 0.5, l = 0.5, "cm"))


FigList=list(DECON.PLOT,Wilcoxon.Results.PLOT)

cowplot::plot_grid(plotlist = FigList,
                           ncol = 1,
                           align = "v", axis = "lr",
                           rel_heights = c(1, 0.65))
##############################################################################


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


par(mar = c(5, 23, 4, 4))
barplot(height=A2$Log10, names=A2$parent_term,horiz=T,las = 1,cex.names=0.7,xlim=c(-15,10),xlab="-Log10(FDR)",main="GO Parental terms",
col =ifelse(A2$Log10>0,"salmon","lightblue"))













DOWN3.BP.All.E=DOWN3[DOWN3$gene_list == "All-E.vs.All-N" & DOWN3$go_type == "BP" ,]
DOWN3.BP.All.E %>% group_by(parent_term) %>% summarise(fdr = min(fdr)) %>%
mutate(parent_term = fct_reorder(parent_term, -fdr)) %>%
ggplot(DOWN3.BP.All.E,mapping = aes(x = -log10(fdr), y = parent_term)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) + theme_bw()

barplot(height=DOWN3.BP.All.E$fdr, names=DOWN3.BP.All.E$parent_term,col="#69b3a2",horiz=T, las=1)




A=rbind(UP3.BP.All.E2,DOWN3.BP.All.E2)
