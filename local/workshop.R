# install Bioconductor to install libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")
BiocManager::install(c("edgeR","RColorBrewer","org.Mm.eg.db","AnnotationDbi","svglite","EnhancedVolcano"))




# load all libraries
library(edgeR)
library(gplots)
library(RColorBrewer)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(EnhancedVolcano)

#import read counts####

ND1<-read.delim("NTC_D1.csv",header=F)
ND2<-read.delim("NTC_D2.csv",header=F)
ND3<-read.delim("NTC_D3.csv",header=F)
NG1<-read.delim("NTC_G1.csv",header=F)
NG2<-read.delim("NTC_G2.csv",header=F)
NG3<-read.delim("NTC_G3.csv",header=F)
DD1<-read.delim("DKO_D1.csv",header=F)
DD2<-read.delim("DKO_D2.csv",header=F)
DD3<-read.delim("DKO_D3.csv",header=F)
DG1<-read.delim("DKO_G1.csv",header=F)
DG2<-read.delim("DKO_G2.csv",header=F)
DG3<-read.delim("DKO_G3.csv",header=F)
CD1<-read.delim("CGAS_D1.csv",header=F)
CD2<-read.delim("CGAS_D2.csv",header=F)
CD3<-read.delim("CGAS_D3.csv",header=F)
CG1<-read.delim("CGAS_G1.csv",header=F)
CG2<-read.delim("CGAS_G2.csv",header=F)
CG3<-read.delim("CGAS_G3.csv",header=F)
RD1<-read.delim("RIG_D1.csv",header=F)
RD2<-read.delim("RIG_D2.csv",header=F)
RD3<-read.delim("RIG_D3.csv",header=F)
RG1<-read.delim("RIG_G1.csv",header=F)
RG2<-read.delim("RIG_G2.csv",header=F)
RG3<-read.delim("RIG_G3.csv",header=F)


#
#Make a dataframe of all read count tables
counts.org <- data.frame(
  row.names=gsub("\\..*","",ND1[,1]),
  ND1=ND1[,2],ND2=ND2[,2],ND3=ND3[,2],NG1=NG1[,2],NG2=NG2[,2],NG3=NG3[,2],
  DD1=DD1[,2],DD2=DD2[,2],DD3=DD3[,2],DG1=DG1[,2],DG2=DG2[,2],DG3=DG3[,2],
  CD1=CD1[,2],CD2=CD2[,2],CD3=CD3[,2],CG1=CG1[,2],CG2=CG2[,2],CG3=CG3[,2],
  RD1=RD1[,2],RD2=RD2[,2],RD3=RD3[,2],RG1=RG1[,2],RG2=DG2[,2],RG3=RG3[,2]
)

# remove last five statistics
counts <- head(counts.org,-5)


#edgeR####
# sample info
meta <-data.frame(
  row.names=colnames(counts),
  treatment=c(rep(c(rep("DMSO",3),rep("GSK",3)),4)),
  replicate=
    c("1","2","3","1","2","3","4","5","6","4","5","6",
      "7","8","9","7","8","9","10","11","12","10","11","12"),
  genotype=c(rep("NTC",6),rep("DKO",6),
             rep("CGAS",6),rep("RIG",6))
  
)

group<- factor(paste(meta$treatment,meta$genotype,sep="."))
cbind(meta,group=group)
treatment<-factor(meta$treatment)
replicate<-factor(meta$replicate)
genotype<-factor(meta$genotype)

# paired t-test model: common effect of treatment considering paired samples
design<-model.matrix(~0+group)
colnames(design)<-levels(group)


# First step of edgeR
b<-DGEList(counts=counts, group=group)
keep<- filterByExpr(b, min.count=100)
b <- b[keep, , keep.lib.sizes=F]
b<- calcNormFactors(b)
# Create a PCA plot equivalent
#tiff("PCA.tiff")
plotMDS(b, col=rep(1:2),each=12,cex=0.4)
#dev.off()



# general linear model
b.glm<-estimateDisp(b,design)

# how well does the model fit?
plotBCV(b.glm)


# Determine DGEs using Likelihood Ratio test
b.glm.fit<-glmFit(b.glm,design)

my.contrasts<- makeContrasts(
  NTC.GvD= GSK.NTC-DMSO.NTC,
  DKO.GvD= GSK.DKO-DMSO.DKO,
  CGAS.GvD= GSK.CGAS-DMSO.CGAS,
  RIG.GvD = GSK.RIG-DMSO.RIG,
  NTC.GvRIG.G = (GSK.NTC-DMSO.NTC)-(GSK.RIG-DMSO.RIG),
  NTC.GvCGAS.G = (GSK.NTC-DMSO.NTC)-(GSK.CGAS-DMSO.CGAS),
  NTC.GvDKO.G = (GSK.NTC-DMSO.NTC)-(GSK.DKO-DMSO.DKO),
  
  levels=design
)


#### Change ENSEMBL to SYMBOL
rownames(b.glm.fit) <- unname(mapIds(org.Mm.eg.db, keys=rownames(b.glm.fit), keytype = "ENSEMBL", column="SYMBOL", columns="SYMBOL", multiVals="asNA"))
# remove non-unique rownames/gene symbols
omitNArows <- na.omit(rownames(b.glm.fit))
b.glm.fit <- b.glm.fit[omitNArows,]

# DGE from a given comparison
NTC.GvD.lrt <- glmLRT(b.glm.fit, contrast=my.contrasts[,"NTC.GvD"])
DKO.GvD.lrt <- glmLRT(b.glm.fit, contrast=my.contrasts[,"DKO.GvD"])
CGAS.GvD.lrt <- glmLRT(b.glm.fit,contrast=my.contrasts[,"CGAS.GvD"])
RIG.GvD.lrt <- glmLRT(b.glm.fit,contrast=my.contrasts[,"RIG.GvD"])
NTC.GvRIG.G.lrt <- glmLRT(b.glm.fit, contrast=my.contrasts[,"NTC.GvRIG.G"])
NTC.GvCGAS.G.lrt <- glmLRT(b.glm.fit, contrast=my.contrasts[,"NTC.GvCGAS.G"])
NTC.GvDKO.G.lrt <- glmLRT(b.glm.fit, contrast=my.contrasts[,"NTC.GvDKO.G"])



# plot MD
#tiff("combined_WT_DMSOvGSK.tiff",units='in',width=5,height=5,res=1200)
#b.wt.DvG.md<-plotMD(NTC.GvD.lrt, xlim=c(-5,15),cex=0.5,main=NULL)
#dev.off()
#tiff("combined_RLR_DMSOvGSK.tiff",units='in',width=5,height=5,res=1200)
#b.mut.DvG.md<-plotMD(b.dko.gskvdmso.lrt, xlim=c(-5,15),cex=0.5,main=NULL)
#dev.off()
#####


summary(decideTests(NTC.GvD.lrt))
summary(decideTests(DKO.GvD.lrt))
summary(decideTests(CGAS.GvD.lrt))
summary(decideTests(RIG.GvD.lrt))
summary(decideTests(NTC.GvRIG.G.lrt))
summary(decideTests(NTC.GvCGAS.G.lrt))
summary(decideTests(NTC.GvDKO.G.lrt))



# sort by logFC all genes not just DGEs
b.ntc.gskvdmso.dge <- topTags(NTC.GvD.lrt, n=803, sort.by="PVal", p.value=0.01, adjust.method="fdr")
b.rig.gskvdmso.dge <- topTags(RIG.GvD.lrt, n=3641, sort.by="PVal",p.value=0.00001, adjust.method="fdr")
b.cgas.gskvdmso.dge <- topTags(CGAS.GvD.lrt, n=2134, sort.by="PVal",p.value=0.01, adjust.method="fdr")
b.dko.gskvdmso.dge <- topTags(DKO.GvD.lrt, n=1147, sort.by="PVal",p.value=0.01, adjust.method="fdr")

b.ntc.gskvdmso.dge.fcsort <- b.ntc.gskvdmso.dge[order(b.ntc.gskvdmso.dge$table$logFC,decreasing=TRUE),]
b.rig.gskvdmso.dge.fcsort <- b.rig.gskvdmso.dge[order(b.rig.gskvdmso.dge$table$logFC,decreasing=TRUE),]
b.cgas.gskvdmso.dge.fcsort <- b.cgas.gskvdmso.dge[order(b.cgas.gskvdmso.dge$table$logFC,decreasing=TRUE),]
b.dko.gskvdmso.dge.fcsort <- b.dko.gskvdmso.dge[order(b.dko.gskvdmso.dge$table$logFC,decreasing=TRUE),]

# subset positive FC 
b.ntc.gskvdmso.dge.up<-b.ntc.gskvdmso.dge.fcsort[b.ntc.gskvdmso.dge.fcsort$table$logFC>0,]
b.rig.gskvdmso.dge.up<-b.rig.gskvdmso.dge.fcsort[b.rig.gskvdmso.dge.fcsort$table$logFC>0,]
b.cgas.gskvdmso.dge.up<-b.cgas.gskvdmso.dge.fcsort[b.cgas.gskvdmso.dge.fcsort$table$logFC>0,]
b.dko.gskvdmso.dge.up<-b.dko.gskvdmso.dge.fcsort[b.dko.gskvdmso.dge.fcsort$table$logFC>0,]

# set colour
colfunc<-colorRampPalette(c("blue","black","red"))


# Volcano plots
w.res <- topTags(b.ntc.gskvdmso.lrt, sort.by="logFC", n="Inf")
w.lab <- rownames(topTags(b.ntc.gskvdmso.lrt,sort.by="PVal",n=50))
svg("WT_volcano.svg", family='sans', width=11, height=7)
EnhancedVolcano(data.frame(w.res),
                lab=rownames(w.res),
                x='logFC',
                y='FDR',
                selectLab=w.lab,
                pCutoff=0.05,
                FCcutoff=F,
                gridlines.minor=F,
                gridlines.major=F,
                legendPosition='none',
                ylab="-logFDR",
                title=NULL,
                subtitle=NULL,
                labSize=5,
                xlim=c(-4,4),
                ylim=c(0,18),
                axisLabSize=18,
                drawConnectors=F,
                col=c('black','black','black','red'),
                colAlpha=1
)
dev.off()

r.res <- topTags(b.rlr.gskvdmso.lrt, sort.by="logFC", n="Inf")
#r.lab <- rownames(topTags(b.rlr.gskvdmso.lrt,sort.by="PVal",n=50))
svg("Mut_volcano.svg", family='sans', width=11, height=7)
EnhancedVolcano(data.frame(r.res),
                lab=rownames(r.res),
                x='logFC',
                y='FDR',
                selectLab="",
                pCutoff=0.05,
                FCcutoff=F,
                gridlines.minor=F,
                gridlines.major=F,
                legendPosition='none',
                ylab="-logFDR",
                title=NULL,
                subtitle=NULL,
                labSize=NULL,
                xlim=c(-2.25,2.25),
                ylim=c(0,18),
                axisLabSize=18,
                drawConnectors=F,
                col=c('black','black','black','red'),
                colAlpha=1
)
dev.off()



## deprecated - down regulated genes separately
#####
# dn regulated genes together and add rowsep 3 colors
b.dn.combined <- c(b.dn.in.wt.not.rlr, b.dn.in.both, b.dn.in.rlr.not.wt)
svglite::svglite(file="dn.combined_heatmap-final.svg",width=7,height=6,
                 system_fonts='Arial',bg='white')
heatmap.2(cpm(b[,c(1:12)])[b.dn.combined,],
          col=colfunc2(8), scale="row", Rowv=NA, Colv=NA, 
          cexCol=1,srtCol=0, adjCol=c(0.5,0), labRow  = NA,
          density.info = "none",trace="none", dendrogram = "none",
          symkey=FALSE,symbreaks=TRUE,revC = FALSE,
          lhei=c(1.2,5), lwid=c( 5, 10),margins=c(1,1),cexRow=1,
          rowsep=c(464,642))
dev.off()
#####




### VennDiagram of common vs. uniquely up/downregulated genes
# commonly up and downregulated genes upon GSK343 treatment
#wt.gvd.lrt<-glmLRT(b.glm.fit.noNA, contrast=my.contrasts[,"WT.GvD"])
dif.wt.gvd<-topTags(b.wt.gskvdmso.lrt,n=Inf,p=0.05)$table
up.dif.wt.gvd <- row.names(dif.wt.gvd[dif.wt.gvd$logFC>0,])
down.dif.wt.gvd <- row.names(dif.wt.gvd[dif.wt.gvd$logFC<0,])

#rlr.gvd.lrt<-glmLRT(b.glm.fit.noNA, contrast=my.contrasts[,"RLR.GvD"])
dif.rlr.gvd<-topTags(b.rlr.gskvdmso.lrt,n=Inf,p=0.05)$table
up.dif.rlr.gvd <- row.names(dif.rlr.gvd[dif.rlr.gvd$logFC>0,])
down.dif.rlr.gvd <- row.names(dif.rlr.gvd[dif.rlr.gvd$logFC<0,])

venn.up<- venn.diagram(
  
  x=list(up.dif.wt.gvd,up.dif.rlr.gvd),
  category.names=c("",""),
  #resolution=600,
  filename= NULL,
  euler.d=T,
  output=T,
  scaled=T,
  cex=2.5,
  lwd=4,
  col="transparent",
  
  fill=c("blue","red"),
  ext.dist=0.05,
  main.fontfamily = 'sans'
)
ggsave(venn.up, file="venn_up.svg",device="svg",height=4,width=4) 
dev.off()

venn.down <- venn.diagram(
  x=list(down.dif.wt.gvd,down.dif.rlr.gvd),
  category.names=c("",""),
  #resolution=600,
  filename= NULL,
  euler.d=T,
  output=T,
  scaled=T,
  cex=2.5,
  lwd=4,
  col="transparent",
  
  fill=c("blue","red"),
  ext.dist=0.05,
  main.fontfamily='sans'
)

ggsave(venn.down, file="venn_down.svg",device="svg",height=4,width=4)
dev.off()






# Creating files for GSEA (Broad) 
#####
# a revelation: gene enrichment analysis uses the whole normalized count table, not just the DGE list

cpm_NTC.GvD<-cpm(b[,c(1,2,3,4,5,6)])
#rownames(table)<- unname(mapIds(org.Mm.eg.db,rownames(table),"SYMBOL","ENSEMBL"))
write.table(cpm_NTC.GvD,file="cpm_NTC.GvD.txt",quote=F,sep="\t")

cpm_DKO.GvD<-cpm(b[,c(7,8,9,10,11,12)])
write.table(cpm_DKO.GvD,file="cpm_DKO.GvD.txt",quote=F,sep="\t")

cpm_CGAS.GvD<-cpm(b[,c(13,14,15,16,17,18)])
write.table(cpm_CGAS.GvD,file="cpm_CGAS.GvD.txt",quote=F,sep="\t")

cpm_RIG.GvD<-cpm(b[,c(19,20,21,22,23,24)])
write.table(cpm_RIG.GvD,file="cpm_RIG.GvD.txt",quote=F,sep="\t")
#####


## Expression of Cgas, Ifih1 and Ddx58 
b.target.not.removed<-DGEList(counts=b.counts, group=group)
b.target.not.removed<- calcNormFactors(b.target.not.removed)
cpm(b.target.not.removed)[c("ENSMUSG00000026896","ENSMUSG00000040296","ENSMUSG00000032344"),]
