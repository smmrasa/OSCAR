#making normalized count from HTseq counts
ss<-read.csv("SampleSheet.csv") #a csv file contain file nema (HTseq Count) the sample name and groups 

sampleTable <- data.frame(sampleName = ss$SampleName,
                          fileName = ss$FileName,
                          condition = ss$Group)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = "./HTseqCount",
                                  design= ~ condition)
dds <- DESeq(dds)
count<-counts(dds,normalized=TRUE)
write.csv(count,"./Counts_norm_DESeq2_AllGenes.csv")

#####################################################

#Figure 4e Person correlation 
count<-read.csv("Counts_norm_DESeq2_AllGenes.csv",row.names = 1,check.names = FALSE)
df<-count[apply(count,1,IQR)>0.1,] #selecting the genes with IQR>0.1
corM<-matrix(rep(0,ncol(df)^2),ncol=ncol(df))
colnames(corM)<-colnames(df)
rownames(corM)<-colnames(df)
for(i in 1:ncol(df)){
  for(j in 1:ncol(df)){
    corM[i,j]<-cor(df[,i],df[,j])
  }
}
pdf("correlation.pdf")
pheatmap(corM)
dev.off()

####################################################################

#Figure 4F
countPCA <- count[apply(count > 10, 1, sum)>2.99 & apply(log2(count+0.1), 1, IQR) >= 1.5, ]
t.log<-log(countPCA)
df<-t.log[ , ! apply( t.log , 2 , function(x) all(is.na(x)) ) ]
df <- df[!is.infinite(rowSums(df)),]
pca<-prcomp(t(df),center = TRUE,scale. = TRUE)
percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
percentage <- c(paste( "PC1(", as.character(percentage)[1] ,"%)", sep=""),paste( "PC2(", as.character(percentage)[2] ,"%)", sep="")) 
p1<-cbind(as.data.frame(pca$x[,1:2]),Group=ss$Group)
rownames(p1)<-ss$SampleName

pdf("PCAplot_NormalizedFilteredCounts.pdf")

p<-ggplot(p1,aes(x=p1[,1],y=p1[,2],col=Group))+
  geom_point(aes(shape=Group))+
  xlab(percentage[1])+
  ylab(percentage[2])

print(p)

dev.off()


#################################################

#Figure 5A, 5E, S5D
#pipeline for gene set analysis

if("ggplot2" %in% rownames(installed.packages()) == FALSE){
  install.packages("ggplot2",repos='http://cran.us.r-project.org')}
if("grid" %in% rownames(installed.packages()) == FALSE){
  install.packages("grid",repos='http://cran.us.r-project.org')}
if("gridExtra" %in% rownames(installed.packages()) == FALSE){
  install.packages("grid",repos='http://cran.us.r-project.org')}

library(ggplot2)
library(gridExtra)
library(grid)
dir.create("./csv")
l<-list.files()
cat("Please select the Sample Sheet CSV file (The first column should be the sample`s names and the second column should be the group`s names) \n" )
print(l)
ls<-l[as.numeric(readline(prompt="Enter a number: "))]
cat("Please select the Count csv file (the first column should be the gene name following by one column for every sample) \n" )
print(l)
lc<-l[as.numeric(readline(prompt="Enter a number: " ))]

cat("Please select the Gene Set CSV file (the first column should be the gene set name, the second column should be the gene symbol and the thisrd column should be the 
    direction (up or down regulated genes with 1 or -1 numbers)) \n" )
print(l)
lgs<-l[as.numeric(readline(prompt="Enter a number: "))]



cat("Show outliers in Boxplot(y/n)")
outl<-readline(prompt="y or n: ")
if(outl=="n"){outl<-NA} else {outl<-16}
cat("Fixed scales in Boxplot(y/n)")
fixs<-readline(prompt="y or n: ")
mins<-NA
maxs<-NA
if(fixs=="y"){
  
  mins<-as.numeric(readline(prompt="Minimum: "))
  maxs<-as.numeric(readline(prompt="Maximum: "))
}
cat("Cutoff for the minimum number of reads in 20 % of samples:")
minfilter<-readline(prompt="Enter a number: ")

ss<-read.csv(ls)
colnames(ss)<-c("SampleName","Group")

counts<-read.csv(lc,row.names = 1,check.names = FALSE)

if(!all(ss$SampleName %in% colnames(counts))){stop("Sample names does not match with the Counts!!")}
counts<-counts[,colnames(counts) %in% ss$SampleName]
counts<-counts[,match(ss$SampleName,colnames(counts))]
gn<-read.csv(lgs,stringsAsFactors = FALSE)
colnames(gn)<-c("setName", "gene","Direction")

if(!  all(gn$gene %in% rownames(counts))){
  cat("The following genes are not in the count list: \n")
  print(gn[! gn$gene %in% rownames(counts),])
}
gn<-gn[gn$gene %in% rownames(counts) & !is.na(gn$gene),]
pdf("Results.pdf")
totalres<-data.frame()
for(gs in unique(gn$setName)){
  gsel<-gn[gn$setName==gs,]
  csel<-counts[rownames(counts) %in% gsel$gene,]
  csel<-csel[apply(csel > minfilter, 1, sum)>ncol(csel)*0.2,]
  write.csv(csel,paste("./csv/",gs,"_Counts.csv",sep=""))
  if(nrow(csel)==0){next}
  gsel<-gsel[gsel$gene %in% rownames(csel),]
  gsel<-gsel[match(rownames(csel),gsel$gene),]
  csel<-t(scale(t(csel)))
  write.csv(csel,paste("./csv/",gs,"_Counts_Scaled.csv",sep=""))
  csel<-csel*gsel$Direction
  write.csv(csel,paste("./csv/",gs,"_Counts_Scaled_MultipliedByDirection.csv",sep=""))
  csel<-as.data.frame(csel)
  res<-data.frame()
  for(gr in unique(ss$Group)){
    sam<-ss[ss$Group==gr,]
    if(sum(colnames(csel) %in% sam$SampleName)<2){
      m<-csel[,colnames(csel) %in% sam$SampleName]
    } else {
      m<-apply(csel[,colnames(csel) %in% sam$SampleName],1,mean)
    }
    res<-rbind(res,data.frame(MeanZscore=m,group=gr))
  }
  res$group<-factor(res$group,levels=unique(res$group))
  p<-ggplot(res,aes(y=MeanZscore,x=group))+
    geom_boxplot(outlier.shape =outl )+
    theme(text=element_text(size=15))+
    scale_y_continuous(limits=c(mins,maxs))+
    ggtitle(gs)
  print(p)
  stat<-data.frame()
  for(i in 1:(length(levels(res$group))-1)){
    g1<-levels(res$group)[i]
    
    for(j in (i+1):length(levels(res$group))){
      g2<-levels(res$group)[j]
      pvalueg<-wilcox.test(res[res$group==g2,"MeanZscore"],res[res$group==g1,"MeanZscore"],alternative="greater",paired = TRUE)$p.value
      pvalueles<-wilcox.test(res[res$group==g2,"MeanZscore"],res[res$group==g1,"MeanZscore"],alternative="less",paired = TRUE)$p.value
      pvalue2side<-wilcox.test(res[res$group==g2,"MeanZscore"],res[res$group==g1,"MeanZscore"],alternative="two.sided",paired = TRUE)$p.value
      
      stat<-rbind(stat,data.frame(Group1=g1,Group2=g2,Wilxox_pvalue_greater=pvalueg,Wilxox_pvalue_less=pvalueles,Wilxox_pvalue_2side=pvalue2side))  
    }
  }
  grid.newpage()
  if(nrow(stat)>25){
    np<-seq(-1,nrow(stat),25)
    if(np[length(np)]!=nrow(stat)){np<-c(np,nrow(stat))}
    for(page in 2:length(np)){
      grid.newpage()
      grid.table(stat[(np[page-1]+1):np[page],], rows = NULL,theme=ttheme_default(base_size = 10))  
    }
  }else{
    grid.newpage()
    grid.table(stat, rows = NULL,theme=ttheme_default(base_size = 10))  
  }
  
  
  totalres<-rbind(totalres,cbind(res,GeneSet=gs))
}
totalres$GeneSet<-factor(totalres$GeneSet,levels=unique(totalres$GeneSet))
p<-ggplot(totalres,aes(x=GeneSet,y=MeanZscore,fill=group))+
  geom_boxplot(outlier.shape =outl )+
  scale_y_continuous(limits=c(mins,maxs))+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
print(p)
dev.off()

####################################################################################
#Figure 5g idetification of the cell cycle 

fpkm<-read.csv("Gene_exp_FPKM.csv",check.names = FALSE) # gene expression in FPKM (Fragment per kilo base pair transcript per million sequenced reads)

library(SingleCellExperiment)
fpkm<-as.matrix(fpkm)
class(fpkm) <- "numeric"
sce <- SingleCellExperiment(list(counts=fpkm))
set.seed(100)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
library(scran)
assignments <- cyclone(sce, mm.pairs)

ph<-data.frame(Sample=colnames(fpkm),phase=assignments$phases,assignments$normalized.scores)
ph<-ph[order(ph$phase,ph$Sample),]
write.csv(ph,"CellCyclePhase.csv")

########################################################


