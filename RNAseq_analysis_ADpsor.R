### using DESeq2 to get expression matrix
library(DESeq2)

### or other locations storing the gene expression files output by HTSeq
directory <- "../HTSeq/"
sampleFiles <- list.files(directory)


### annotate the samples:
### e.g. sampleTable <- data.frame(sampleName=nameVariable,fileName=sampleFiles,Skin=skinVariable,patient=patientVariable,Sex=sexVariable)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=directory,design=~Skin)
rawddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=directory,design=~Skin)

### remove genes with on average <= 1 read per sample
### from 60,991 genes to 31,364 genes

ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq))>=dim(ddsHTSeq)[2],]


### reorder level of skin type (control first)

ddsHTSeq$Skin <- factor(ddsHTSeq$Skin,levels= c("CO_healthy","PSO_non-lesional","PSO_lesional","AD_non-lesional","AD_lesional"))


### normalization
ddsHTSeq <- DESeq(ddsHTSeq)


### differential expression

tempcomparison <- matrix(character(),7,2)

### psoriasis
tempcomparison[1,] <- c("CO_healthy","PSO_non-lesional")
tempcomparison[2,] <- c("CO_healthy","PSO_lesional")
tempcomparison[3,] <- c("PSO_non-lesional","PSO_lesional")

### AD
tempcomparison[4,] <- c("CO_healthy","AD_non-lesional")
tempcomparison[5,] <- c("CO_healthy","AD_lesional")

tempcomparison[6,] <- c("AD_non-lesional","AD_lesional")


### psoriasis vs AD
tempcomparison[7,] <- c("PSO_lesional","AD_lesional")



res <- list()

for (i in 1:dim(tempcomparison)[1]){
  print(i)
  
  tempdds <- ddsHTSeq[,(ddsHTSeq$Skin==tempcomparison[i,1]) | (ddsHTSeq$Skin==tempcomparison[i,2])]
  
  tempdds$Skin <- droplevels(tempdds$Skin)
  tempdds$patient <- droplevels(tempdds$patient)
  tempdds$Sex <- droplevels(tempdds$Sex)
  
  design(tempdds) <- ~ Sex+Skin
  
  tempdds <- DESeq(tempdds)
  res[[paste(tempcomparison[i,],collapse="_VS_")]] <- results(tempdds,addMLE=T,contrast=c("Skin",tempcomparison[i,2] ,tempcomparison[i,1]))
}


### DEGs: FDR <= 0.05 and |log2FC|>=1
res_DEGs <- list()

for (i in 1:length(res)){
  tempa <- res[[i]]
  temp1 <- tempa[!is.na(tempa$padj) & (tempa$padj<=0.05) & ((tempa$log2FoldChange>=1) | (tempa$log2FoldChange <= -1)),]
  temp2 <- tempa[!is.na(tempa$padj) & (tempa$padj<=0.05) & ((tempa$log2FoldChange >=1)),]
  temp3 <- tempa[!is.na(tempa$padj) & (tempa$padj<=0.05) & ((tempa$log2FoldChange <= -1)),]
  res_DEGs[[names(res)[i]]] <- temp1
  print(c(names(res_DEGs)[i],dim(temp1)[1],dim(temp2)[1],dim(temp3)[1]))
}




### using voom + limma 
temphtseqCount <- counts(ddsHTSeq,normalized=F)
library(edgeR)
DGEList <- DGEList(counts=temphtseqCount,genes=rownames(temphtseqCount))

# apply scale normalization
DGEList_scale <- calcNormFactors(DGEList)

tempdesign  <- model.matrix(~factor(ddsHTSeq$Skin)+factor(ddsHTSeq$Sex))

# use voom to convert the read counts to log2-cpm, with associated weights, ready for linear
DGEList_scale_voom <- voom(DGEList_scale,tempdesign,plot=F)


ddsHTSeq$MergedSkin <- ddsHTSeq$Skin



### differential expression 

tempcomparison <- matrix(character(),8,2)

### psoriasis
tempcomparison[1,] <- c("CO_healthy","PSO_non-lesional")
tempcomparison[2,] <- c("CO_healthy","PSO_lesional")
tempcomparison[3,] <- c("PSO_non-lesional","PSO_lesional")


### AD
tempcomparison[4,] <- c("CO_healthy","AD_non-lesional")
tempcomparison[5,] <- c("CO_healthy","AD_lesional")
tempcomparison[6,] <- c("AD_non-lesional","AD_lesional")


### psoriasis vs AD
tempcomparison[7,] <- c("PSO_non-lesional","AD_non-lesional")
tempcomparison[8,] <- c("PSO_lesional","AD_lesional")



library(limma)
limmavoom.res <- list()

for (i in 1:dim(tempcomparison)[1]){
  print(i)
  
  s="MergedSkin";
  
  templogical <- (ddsHTSeq[[s]]==tempcomparison[i,1]) | (ddsHTSeq[[s]]==tempcomparison[i,2])
  
  tempdesign  <- model.matrix(~factor(ddsHTSeq[[s]][templogical],levels=c(tempcomparison[i,1],tempcomparison[i,2]))+factor(ddsHTSeq$Sex[templogical]))
  
  
  # use voom() to convert the read counts to log2-cpm, with associated weights, ready for linear
  tempDGEList_scale <- DGEList_scale[,templogical]
  
  tempDGEList_scale_voom <- voom(tempDGEList_scale,tempdesign,plot=F)
  
  tempfit.voom <- lmFit(tempDGEList_scale_voom,tempdesign)
  tempfit.voom <- eBayes(tempfit.voom)
  limmavoom.res[[paste(tempcomparison[i,],collapse="_VS_")]] <- topTable(tempfit.voom,coef=2,number=dim(tempDGEList_scale_voom)[1],sort= "none")
  
}

### DEGs: FDR <= 0.05 and |log2FC|>=1
limmavoom.res_DEGs <- list()

for (i in 1:length(limmavoom.res)){
  tempa <- limmavoom.res[[i]]
  temp1 <- tempa[!is.na(tempa$adj.P.Val) & (tempa$adj.P.Val<=0.05) & ((tempa$logFC>=1) | (tempa$logFC <= -1)),]
  temp2 <- tempa[!is.na(tempa$adj.P.Val) & (tempa$adj.P.Val<=0.05) & ((tempa$logFC >=1)),]
  temp3 <- tempa[!is.na(tempa$adj.P.Val) & (tempa$adj.P.Val<=0.05) & ((tempa$logFC <= -1)),]
  limmavoom.res_DEGs[[names(limmavoom.res)[i]]] <- temp1
  print(c(names(limmavoom.res_DEGs)[i],dim(temp1)[1],dim(temp2)[1],dim(temp3)[1]))
  
}


### PC plot and Pheatmap


### inverse normalization on log values of genes
DGEList_scale_voom_invnorm <- t(apply(DGEList_scale_voom[[3]],1,function(x){ qnorm((rank(x)-(3/8))/(length(x)-2*(3/8)+1))}))
DGEList_scale_voom_invnorm.pca <- prcomp(t(DGEList_scale_voom_invnorm),scale=T)

library(ggplot2)
tempplot <- data.frame(Condition= as.character(ddsHTSeq$MergedSkin),PC1= DGEList_scale_voom_invnorm.pca[[5]][,1],PC2= DGEList_scale_voom_invnorm.pca[[5]][,2], PC3= DGEList_scale_voom_invnorm.pca[[5]][,3], PC4= DGEList_scale_voom_invnorm.pca[[5]][,4])

pdf("DGEList_scale_voom_invnorm.pca.pdf")
ggplot(tempplot,aes(x=PC1,y=PC2))+geom_point(aes(colour=Condition),size=5)
ggplot(tempplot,aes(x=PC1,y=PC3))+geom_point(aes(colour=Condition),size=5)
ggplot(tempplot,aes(x=PC2,y=PC3))+geom_point(aes(colour=Condition),size=5)
dev.off()


library(scatterplot3d)
pdf("DGEList_scale_voom_invnorm.pca_3dscatter.pdf",height=6,width=8)
scatterplot3d(tempplot[,2:4],color= rainbow(7)[as.numeric(tempplot[,1])],pch=16, type="p",cex.symbols=2)

legend("topleft", legend = levels(tempplot[,1]), col =  rainbow(7),pch = c(16), inset = -0.08, xpd = TRUE)
dev.off()

### Pheatmap
library(RColorBrewer)
library(pheatmap)

tempdf <- data.frame(Skin=ddsHTSeq$MergedSkin)
rownames(tempdf) <- rownames(colData(ddsHTSeq))
pdf("Pheatmap.pdf",width=8,height=10)
### all genes
pheatmap(DGEList_scale_voom_invnorm,show_rownames=F,show_colnames=F,annotation_col=tempdf)
dev.off()
###  color = colorRampPalette(c("navy", "white", "red"))(50)


### subplots including correlation plots

tempplot <- data.frame(Condition= as.character(ddsHTSeq$MergedSkin),PC1= DGEList_scale_voom_invnorm.pca[[5]][,1],PC2= DGEList_scale_voom_invnorm.pca[[5]][,2], PC3= DGEList_scale_voom_invnorm.pca[[5]][,3], PC4= DGEList_scale_voom_invnorm.pca[[5]][,4])

library(scatterplot3d)
pdf("limmavoom.res_correlation_3dscatter.pdf",width=8,height=8)
### FC vs p-value plot for NN vs psoriatic skin samples
plot(limmavoom.res[[2]][,2],-log10(limmavoom.res[[2]][,5]),xlab="log2FC",ylab="-log10(p)",pch=20,main="control vs psoriatic skin",cex.lab=1.3,xlim=c(-7,7),ylim=c(0,30))

### points(limmavoom.res_DEGs[[2]][,2],-log10(limmavoom.res_DEGs[[2]][,5]),col= "red",pch=20)


### FC vs p-value plot for NN vs AD skin samples
plot(limmavoom.res[[5]][,2],-log10(limmavoom.res[[5]][,5]),xlab="log2FC",ylab="-log10(p)",pch=20, main="control vs AD skin",cex.lab=1.3, xlim=c(-7,7),ylim=c(0,30))

### points(limmavoom.res_DEGs[[5]][,2],-log10(limmavoom.res_DEGs[[5]][,5]),col= "red",pch=20)



### plot FC vs FC
temporder <- rbind(c(2,5),c(1,4))

for (i in 1:dim(temporder)[1]){
  
  plot(limmavoom.res[[temporder[i,1]]][,2],limmavoom.res[[temporder[i,2]]][,2],xlim=c(-8,8),ylim=c(-8,8),pch=19,xlab= paste(names(limmavoom.res)[temporder[i,1]],"(log2FC)"),ylab= paste(names(limmavoom.res)[temporder[i,2]],"(log2FC)"),cex.lab=1.3,cex=1.3,col="grey")
  
  ### DE genes
  temp1 <- rownames(limmavoom.res_DEGs[[temporder[i,1]]])
  temp2 <- rownames(limmavoom.res_DEGs[[temporder[i,2]]])
  
  if (length(temp1)!=0){
    points(limmavoom.res[[temporder[i,1]]][temp1,2],limmavoom.res[[temporder[i,2]]][temp1,2],pch=19,col="orange")
  }
  
  if (length(temp2)!=0){
    points(limmavoom.res[[temporder[i,1]]][temp2,2],limmavoom.res[[temporder[i,2]]][temp2,2],pch=19,col="blue")
  }
  
  if (length(intersect(temp1,temp2))!=0){
    points(limmavoom.res[[temporder[i,1]]][intersect(temp1,temp2),2],limmavoom.res[[temporder[i,2]]][intersect(temp1,temp2),2] ,pch=19,col="red")
  }
  
  lines(c(-10,10),c(-10,10),lty=2)
  abline(h=1,v=1,lty=2,col="red")
  abline(h=-1,v=-1,lty=2,col="red")
  
  ### text(-3,4,paste("cor=",round(cor.test(limmavoom.res[[temporder[i,1]]][,2],limmavoom.res[[temporder[i,2]]][,2],method="spearman")[[4]],digits=2),sep=""),cex=1.5)
}





### 3D scatter plot for inverse normalized explimmavoom.ression data
scatterplot3d(tempplot[,2:4],color= rainbow(7)[as.numeric(tempplot[,1])],pch=16, type="p",cex.symbols=2)

legend("bottomright", legend = levels(tempplot[,1]), col =  rainbow(7),pch = c(16), inset = -0.08, xpd = TRUE,cex=1.1)


dev.off()



### venn diagram, 

library(VennDiagram)

venn.diagram(list("Normal vs Pso LS"= rownames(limmavoom.res_DEGs[[2]]), "Normal vs AD LS"= rownames(limmavoom.res_DEGs[[5]]), "Pso LS vs AD LS"= rownames(limmavoom.res_DEGs[[11]]) ),fill = c("red", "royalblue","yellow"), alpha = c(0.4, 0.4,0.4), cex = 1.5,cat.fontface = 4,lty =2, fontfamily =3,    filename = "Venndiagram.DEGs_lesional.tiff",main="Overlap between DEGs (lesional)");


venn.diagram(list("Normal vs Pso uninvolved"= rownames(limmavoom.res_DEGs[[1]]), "Normal vs AD uninvolved"= rownames(limmavoom.res_DEGs[[4]]), "Pso uninvolved vs AD uninvolved"= rownames(limmavoom.res_DEGs[[10]]) ),fill = c("red", "royalblue","yellow"), alpha = c(0.4, 0.4,0.4), cex = 1.5,cat.fontface = 4,lty =2, fontfamily =3,    filename = "Venndiagram.DEGs_uninvolved.tiff",main="Overlap between DEGs (uninvolved)");



### compute the euclidean distance between all PCs for samples within each disease
tempDist <- matrix(numeric(),0,2)
tempskins <- as.character(unique(colData(ddsHTSeq)$MergedSkin))
tempskins <- tempskins[ (tempskins== "AD_lesional") | (tempskins== "PSO_lesional")]

for (i in tempskins){
  tempR <- colnames(ddsHTSeq)[ddsHTSeq$MergedSkin==i]
  tempD <- DGEList_scale_voom_invnorm.pca[[5]][tempR,1:3]
  
  for (j in 1:(length(tempR)-1)){
    for (k in (j+1):(length(tempR))){
      
      a <- as.numeric(tempD[tempR[j],])
      b <- as.numeric(tempD[tempR[k],])
      
      tempDist <- rbind(tempDist,c(i,sum((a-b)^2)))
      
    }
  }
}


library(ggplot2)

pdf("PCDist_w_skin.pdf",width=8,height=3)
ggplot(data.frame(Skin=tempDist[,1],Dist_log=log10(as.numeric(tempDist[,2]))) ,aes(Dist_log,fill=Skin)
)+geom_density(alpha=0.5,adjust=2)+xlim(0,6)
dev.off()



### compare between our RNA-seq vs previous MA studies which have control vs AD comparison

### obtain summary statistics from GEO for the three AD microarray studies
### for fair comparison, only using genes appearing in both platform
MA_GEO <- list()
MA_GEO[["AD_2009"]] <- as.matrix(read.table("JACI_AD_2009.txt",header=T,sep="\t"))
MA_GEO[["AD_2011"]] <- as.matrix(read.table("JACI_AD_2011.txt",header=T,sep="\t"))
MA_GEO[["AD_2012"]] <- as.matrix(read.table("JACI_AD_201.txt",header=T,sep="\t"))

### only need to correct FC for the first 2 studies
for (i in 1:3){
  MA_GEO[[i]][,6] <- as.numeric(MA_GEO[[i]][,6])*-1
}

### only intersect genes with RNA-seq--> 16,783 genes
tempgenes <- unique(intersect(intersect(intersect(MA_GEO[[1]][,7],MA_GEO[[2]][,7]),MA_GEO[[3]][,7]),rownames(limmavoom.res[[5]])))

MA_GEO_DEGs <- list()
MA_GEO_DEGs$up <- list()
MA_GEO_DEGs$down <- list()
for (i in 1:3){
  MA_GEO_DEGs$up[[names(MA_GEO)[i]]] <- intersect(tempgenes,unique(MA_GEO[[i]][(as.numeric(MA_GEO[[i]][,2])<=0.05) & (as.numeric(MA_GEO[[i]][,6])>=1),7]))
  MA_GEO_DEGs$down[[names(MA_GEO)[i]]] <- intersect(tempgenes,unique(MA_GEO[[i]][(as.numeric(MA_GEO[[i]][,2])<=0.05) & (as.numeric(MA_GEO[[i]][,6])<=-1),7]))
}


MA_GEO_FC <- matrix(matrix(),length(tempgenes),4)
rownames(MA_GEO_FC) <- tempgenes
colnames(MA_GEO_FC) <- c("ThisStudy",names(MA_GEO))
MA_GEO_FC[,1] <- limmavoom.res[[5]][tempgenes,2]
MA_GEO_FC[,2] <- as.numeric(MA_GEO[[1]][match(tempgenes,MA_GEO[[1]][,7]),6])
MA_GEO_FC[,3] <- as.numeric(MA_GEO[[2]][match(tempgenes,MA_GEO[[2]][,7]),6])
MA_GEO_FC[,4] <- as.numeric(MA_GEO[[3]][match(tempgenes,MA_GEO[[3]][,7]),6])

summary(limmavoom.res[[5]][tempgenes,3])



pdf("MA_GEO_FC.pdf",height=12,width=4)
par(mfrow=c(3,1))
plot(MA_GEO_FC[,1],MA_GEO_FC[,2],xlim=c(-7,7),ylim=c(-7,7),col="grey",cex=1,pch=20,xlab= "log2(FC) in AD skin (this study)",ylab= "log2(FC) in AD skin (Guttman, et al. 2009)")
points(MA_GEO_FC[temp1,1],MA_GEO_FC[temp1,2],col="red",cex=0.8,pch=20)

plot(MA_GEO_FC[,1],MA_GEO_FC[,3],xlim=c(-7,7),ylim=c(-7,7),col="grey",cex=1,pch=20,xlab= "log2(FC) in AD skin (this study)",ylab= "log2(FC) in AD skin (Suarez, et al. 2011)")
points(MA_GEO_FC[temp1,1],MA_GEO_FC[temp1,3],col="red",cex=0.8,pch=20)

plot(MA_GEO_FC[,1],MA_GEO_FC[,4],xlim=c(-7,7),ylim=c(-7,7),col="grey",cex=1,pch=20,xlab= "log2(FC) in AD skin (this study)",ylab= "log2(FC) in AD skin (Gittler, et al. 2012)")
points(MA_GEO_FC[temp1,1],MA_GEO_FC[temp1,4],col="red",cex=0.8,pch=20)
dev.off()

_______________________________________________________________________________

### histogram plots to indicate the cytokine signatures for different AD/Pso/control/etc samples
### use Cytokine-upregulated genes
templist <- list()
for (i in names(cytokineinduced)){
  templist[[i]] <- rownames(cytokineinduced[[i]][cytokineinduced[[i]][,2]>0,])
}

Cytokine.induced_Samplewise.FC <- matrix(character(),153,11)
rownames(Cytokine.induced_Samplewise.FC) <- colnames(ddsHTSeq)
Cytokine.induced_Samplewise.FC[,1] <- as.character(ddsHTSeq$MergedSkin)
Cytokine.induced_Samplewise.FC  <- Cytokine.induced_Samplewise.FC[Cytokine.induced_Samplewise.FC[,1]!="AD_acute lesion",]
colnames(Cytokine.induced_Samplewise.FC) <- c("Skin",names(cytokineinduced))

for (i in names(cytokineinduced)){
  tempdata <- 2^DGEList_scale_voom[[3]][templist[[i]],rownames(Cytokine.induced_Samplewise.FC)]
  
  ### take the median of the control samples
  tempdata_normal <- apply(tempdata[,Cytokine.induced_Samplewise.FC[,1]=="CO_healthy"],1,function(x){summary(x)[3]})
  tempdata <- tempdata/tempdata_normal
  
  ### take the gene 75% quantile (as we would not assume 50% of cytokine up-regulated genes also up in lesional skin)
  tempdata <- apply(tempdata,2,function(x){summary(x)[5]})
  Cytokine.induced_Samplewise.FC[,i] <- tempdata
}

### plot histogram
tempplot <- data.frame(Skin=Cytokine.induced_Samplewise.FC[,1] ,
                       IL4_induced= as.numeric(Cytokine.induced_Samplewise.FC[,2]), 
                       IL13_induced= as.numeric(Cytokine.induced_Samplewise.FC[,3]), 
                       IL17A_induced= as.numeric(Cytokine.induced_Samplewise.FC[,4]), 
                       TNFa_induced= as.numeric(Cytokine.induced_Samplewise.FC[,6]), 
                       IFNa_induced= as.numeric(Cytokine.induced_Samplewise.FC[,7]), 
                       IFNg_induced= as.numeric(Cytokine.induced_Samplewise.FC[,8]),
                       IL36A_induced= as.numeric(Cytokine.induced_Samplewise.FC[,9]), 
                       IL36B_induced= as.numeric(Cytokine.induced_Samplewise.FC[,10]), 
                       IL36G_induced= as.numeric(Cytokine.induced_Samplewise.FC[,11]))
tempplot$Skin  <- factor(tempplot$Skin,levels= c("CO_healthy","PSO_non-lesional","PSO_lesional","AD_non-lesional","AD_lesional"))

tempplot <- tempplot[tempplot[,1]=="PSO_non-lesional" | (tempplot[,1]=="PSO_lesional") | (tempplot[,1]=="CO_healthy"),]

library(ggplot2)

pdf("Cytokine.induced_Samplewise.FC_histogram.pdf",width=10,height=5)
ggplot(tempplot,aes(IL17A_induced,fill=Skin))+geom_density(alpha=0.2,adjust=1,bw=3)+xlim(c(0,50))
ggplot(tempplot,aes(IL4_induced,fill=Skin))+geom_density(alpha=0.2,adjust=0.1,bw=1)+xlim(c(0.5,3))
ggplot(tempplot,aes(IL13_induced,fill=Skin))+geom_density(alpha=0.2,adjust=0.1,bw=1)+xlim(c(0.5,3))
ggplot(tempplot,aes(TNFa_induced,fill=Skin))+geom_density(alpha=0.2,adjust=0.2,bw=1.5)+xlim(c(0.5,5))
ggplot(tempplot,aes(IFNa_induced,fill=Skin))+geom_density(alpha=0.2,adjust=0.5,bw=2)+xlim(c(0.5,15))
ggplot(tempplot,aes(IFNg_induced,fill=Skin))+geom_density(alpha=0.2,adjust=0.5,bw=2)+xlim(c(0.5,7))
ggplot(tempplot,aes(IL36A_induced,fill=Skin))+geom_density(alpha=0.2,adjust=0.5,bw=1)+xlim(c(0.5,7)) 
### CORRELATIONS
plot(tempplot[tempplot[,1]=="PSO_lesional",4],tempplot[tempplot[,1]=="PSO_lesional",5],xlab="IL17A induced",ylab="TNFa induced",pch=19,main="Psoriatic skin samples")
plot(tempplot[tempplot[,1]=="PSO_lesional",4],tempplot[tempplot[,1]=="PSO_lesional",6],xlab="IL17A induced",ylab="IFNa induced",pch=19,main="Psoriatic skin samples")
plot(tempplot[tempplot[,1]=="PSO_lesional",4],tempplot[tempplot[,1]=="PSO_lesional",7],xlab="IL17A induced",ylab="IFNg induced",pch=19,main="Psoriatic skin samples")
plot(tempplot[tempplot[,1]=="PSO_lesional",4],tempplot[tempplot[,1]=="PSO_lesional",2],xlab="IL17A induced",ylab="IL4 induced",pch=19,main="Psoriatic skin samples")
dev.off()



