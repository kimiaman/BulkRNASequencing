#Kimia Mansouri
#BULKRNASeq
#FirstStep:Set your working directory

#Reading the data#optional#
sample1Data=data.frame()
sample1File=readLines("counts/ERR188044.count")
sample1File
for (i in sample1File){
  sample1Data=rbind(sample1Data,i)
}

#Read data
sample1=read.table("counts/ERR188044.count")
sample1=sample1[-grep("__",sample1[,1]),]
dim(sample1)


Counts=read.table("counts/ERR188044.count")[,1]
for (i in dir("counts/",full.names = T)){
  Counts=cbind(Counts,read.table(i)[,2])
}
Counts=Counts[-grep("__",Counts[,1]),]
dim(Counts)
View(Counts)
colnames(Counts)=c("GeneIDs",gsub(".count","",dir("counts/")))
View(Counts)

#Set your group files data
GroupFiles <- read.csv("geuvadis_phenodata.csv")
GroupFiles
GroupFiles[which(GroupFiles$population=="YRI"),1]
GroupFiles[which(GroupFiles$population=="GBR"),1]

Counts=Counts[,c("GeneIDs",GroupFiles[which(GroupFiles$population=="YRI"),1],GroupFiles[which(GroupFiles$population=="GBR"),1])]
Counts=data.frame(Counts)
for (i in 2:ncol(Counts)){
  Counts[,i]=as.integer(Counts[,i])
}
typeof(Counts)
sum(Counts$ERR188044)
86*1000000/sum(Counts$ERR188044)
29*1000000/sum(Counts$ERR188273)
View(Counts)

CountsCPM=Counts
for (i in 2:ncol(CountsCPM)){
  CountsCPM[,i]=CountsCPM[,i]*1000000/sum(Counts[,i])
}
View(CountsCPM)
colSums(Counts[,-1])

#notes about data normalization
#CPM
#RPKM, FPKM
#100000   1000     10000  1000
#1000     ?        100    ?
#10                10 


Counts*1000/GeneLength
 
GTF=readLines("chrX.gtf")
length(GTF)
GTF=GTF[grep("\texon\t",GTF)]
length(GTF)

ChrXGeneLength=c()

for (i in Counts$GeneIDs){
  CurrentGene=GTF[grep(i,GTF)]
  GeneLength=0
  for (j in CurrentGene){
    pos=as.integer(unlist(strsplit(j,split = "\t"))[4:5])
    GeneLength=GeneLength+pos[2]-pos[1]
  }
  ChrXGeneLength=c(ChrXGeneLength,GeneLength)
}
names(ChrXGeneLength)=Counts$GeneIDs
ChrXGeneLength
ChrXGeneLength["NM_033031"]

CountsRPK=Counts
for (i in 2:ncol(CountsRPK)){
  CountsRPK[,i]=CountsRPK[,i]*1000/ChrXGeneLength
}
View(CountsRPK)

CountsRPKM=Counts
for (i in 2:ncol(CountsRPKM)){
  CountsRPKM[,i]=CountsRPKM[,i]*1000000/sum(CountsRPKM[,i])
  CountsRPKM[,i]=CountsRPKM[,i]*1000/ChrXGeneLength
}
View(CountsRPKM)

CountsTPM=Counts
for (i in 2:ncol(CountsTPM)){
  CountsTPM[,i]=CountsTPM[,i]*1000/ChrXGeneLength
  CountsTPM[,i]=CountsTPM[,i]*1000000/sum(CountsTPM[,i])
}

TPM

DESEQ2

#DataNormalization
RPKM=function(Counts,GeneLength){
  CountsRPKM=Counts
  for (i in 2:ncol(CountsRPKM)){
    CountsRPKM[,i]=CountsRPKM[,i]*1000000/sum(CountsRPKM[,i])
    CountsRPKM[,i]=CountsRPKM[,i]*1000/GeneLength
  }
  return(CountsRPKM)
}
RPKM(data.frame(Genes=c("S","M","E","N"),Rep1=c(1000000,2000000,500000,0),Rep2=c(1200000,2500000,800000,0),Rep3=c(3000000,6000000,1500000,100000)),c(2000,4000,1000,10000))
RPKM(data.frame(Genes=c("S","M","E"),Rep1=c(6,8,16),Rep2=c(20,17,63),Rep3=c(30,40,36)),c(500,1000,1500))
TPM=function(Counts,GeneLength){
  CountsTPM=Counts
  for (i in 2:ncol(CountsTPM)){
    CountsTPM[,i]=CountsTPM[,i]*1000/GeneLength
    CountsTPM[,i]=CountsTPM[,i]*1000000/sum(CountsTPM[,i])
  }
  return(CountsTPM)
}
TPM(data.frame(Genes=c("S","M","E","N"),Rep1=c(1000000,2000000,500000,0),Rep2=c(1200000,2500000,800000,0),Rep3=c(3000000,6000000,1500000,100000)),c(2000,4000,1000,10000))

Counts=data.frame(Sample1=c(0,2,33),Sample2=c(10,6,55),Sample3=c(4,12,200))
rownames(Counts)=c("Gene1","Gene2","Gene3")
Counts
# 1) Loge
log(4)
log2(4)
log10(4)
Counts1=log(Counts)

# 2) Average Each row
Counts2=rowMeans(Counts1)

Counts2=c()
for (i in 1:nrow(Counts1)) {
  Counts2=c(Counts2,mean(as.numeric(Counts1[i,])))
}
Counts2

rowSums(Counts1)
rowMeans(Counts1)
colSums(Counts1)
colMeans(Counts1)

# 3) Filter out genes with inf
Counts1=Counts1[-grep("-Inf",Counts1),]
Counts3=Counts2[-which(Counts2==-Inf)]
Counts3

# 4) Subtract Avg log from log value of round 1
Counts4=Counts1-Counts3
Counts4

# 5 ) Mean of Each Sample
Counts5=colMeans(Counts4)
Counts5

# 6 ) Creating Scaling Factors by exponenting Counts5
ScalingFactor=exp(Counts5)
ScalingFactor

# 7) Divide Counts by Scaling factors
CountsNorm=matrix(0,nrow(Counts),ncol(Counts))
for (i in 1:ncol(Counts)){
  CountsNorm[,i]=Counts[,i]/ScalingFactor[i]
}
CountsNorm


BiocManager::install("DESeq2", lib = "R/win-library/")
BiocManager::install(c("RCurl", "S4Vectors"))
install.packages("cachem")
install.packages("colorspace")

install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("DESeq2", force = TRUE)
install.packages(c("bitops", "openssl"))
install.packages("RCurl", repos = "http://cran.us.r-project.org")
install.packages("RCurl", type = "binary")
install.packages("cli")
install.packages("fansi")
install.packages("utf8")



###################second part using Deseq2


library(DESeq2)

install.packages("BiocManager")
BiocManager::install("DESeq2")
require(DESeq2)

Counts=read.table("counts/ERR188044.count")[,2]
for (i in dir("counts/",full.names = T)[-1]){
  Counts=cbind(Counts,read.table(i)[,2])
}
rownames(Counts)=read.table("counts/ERR188044.count")[,1]
colnames(Counts)=gsub(".count","",dir("counts/"))
Counts=Counts[-grep("__",rownames(Counts)),]
dim(Counts)
GroupFiles=read.csv("geuvadis_phenodata.csv")
CountsMatrix=as.matrix(Counts)
Class=GroupFiles$population
Class=data.frame(Population=as.factor(Class))
# Class=data.frame(Type=as.factor(c("A","A","B","A","B","A","A","A","B","B")))
DESeq2Obj=DESeqDataSetFromMatrix(CountsMatrix,Class,~Population)
CountsNorm=DESeq(DESeq2Obj)
CountsNorm$sizeFactor
CountsNorm$Population

# Class=data.frame(Type=as.factor(c("A","A","B","A","B","A","A","A","B","B")))
DESeq2Obj=DESeqDataSetFromMatrix(CountsMatrix,Class,~Population)
NormalizedCounts=data.frame(counts(estimateSizeFactors(DESeq2Obj),normalized=T))
NormalizedCounts$Genes=rownames(NormalizedCounts)
rownames(NormalizedCounts)=NULL
dim(NormalizedCounts)
NormalizedCounts=NormalizedCounts[,c(13,1:12)]
write.table(NormalizedCounts,"NormalizedCounts.txt",row.names = F,quote = F)

DEGresult=data.frame(results(CountsNorm))
DEGresult$padj=p.adjust(DEGresult$pvalue,method = "BH")
DEGresult=DEGresult[order(DEGresult$padj),]
plotCounts(CountsNorm,gene = "NR_131238",intgroup = "Population")
plotCounts(CountsNorm,gene = "NM_013444",intgroup = "Population")
plotMA(CountsNorm)
dim(DEGresult)
View(DEGresult)
DEGresult$Genes=rownames(DEGresult)
DEGresult=DEGresult[,c(7,1:6)]
rownames(DEGresult)=NULL
write.table(DEGresult,"DEGresult.txt",row.names = F,quote = F)

library(ggplot2)
DEGresult$pvalue
ggplot(DEGresult[which(DEGresult$pvalue<1),],aes(x=log2FoldChange,y=-log10(pvalue)))+geom_point()+geom_vline(xintercept = c(-1,1))

DEGresult$Color="black"
DEGresult$Color[which(DEGresult$log2FoldChange>=1 & DEGresult$pvalue<=0.05)]="red"
DEGresult$Color[which(DEGresult$log2FoldChange<=-1 & DEGresult$pvalue<=0.05)]="blue"
View(DEGresult)
ggplot(DEGresult[which(DEGresult$pvalue<1),],aes(x=log2FoldChange,y=-log10(pvalue),colour=Color))+geom_point()+
  geom_vline(xintercept = c(-1,1))+theme_classic()+scale_color_manual(labels=c("nonDEG","DownRegulatedGenes","UpRegulatedGenes"),values = c("green","blue","red"))

View(Counts)

10^(-1)
10^(-2)
10^(-3)

View(DEGresult)


DEGresult=DEGresult[which(DEGresult$pvalue<0.05),]
DEGresult[which(DEGresult$log2FoldChange>=1),c(2,5)]
DEGresult[which(DEGresult$log2FoldChange<=-1),c(2,5)]


FoldChange
Log2FoldChange>1 <-1

Counts
# 1) Normalize Counts
for (i in 1:ncol(Counts)){
  Counts[,i]=Counts[,i]/CountsNorm$sizeFactor[i]
}
View(Counts)

# 2) Mean of each Group
YRI=rowMeans(Counts[,which(CountsNorm$Population=="YRI")])
GBR=rowMeans(Counts[,which(CountsNorm$Population=="GBR")])

# 3) FoldChange
FoldChange=YRI/GBR
FoldChange[1:10]

# 4) Log2FoldChange
log2FoldChange=log2(FoldChange)
log2FoldChange[1:10]

# 5) t.test
Pvalue=c()
YRIindex=which(CountsNorm$Population=="YRI")
GBRindex=which(CountsNorm$Population=="GBR")
for (i in 1:nrow(Counts)){
  Pvalue=c(Pvalue,t.test(Counts[i,YRIindex],Counts[i,GBRindex])$p.value)
}
Pvalue

# 6) DEG
results=data.frame(MeanYRI=YRI,MeanGBR=GBR,FoldChange,log2FoldChange,Pvalue)
View(results)
results=results[order(results$Pvalue),]
View(results)
results=results[which(results$Pvalue<0.05),]
View(results[order(results$log2FoldChange),])
View(results[order(results$log2FoldChange,decreasing = T),])

UpRegulateGenes=which(log2FoldChange>=1 & Pvalue<0.05)
DownRegulateGenes=which(log2FoldChange<=-1 & Pvalue < 0.05)
DEG=which(log2FoldChange>=1 | log2FoldChange<=-1)
length(UpRegulateGenes)
length(DownRegulateGenes)
length(DEG)
View(Counts[UpRegulateGenes,])

































