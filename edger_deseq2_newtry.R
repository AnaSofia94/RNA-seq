#IMPORT DATA
library(tximport)
source("https://bioconductor.org/biocLite.R")
biocLite("tximport")
source("https://bioconductor.org/biocLite.R")
biocLite("rhdf5")
load("C:/Users/Ana Sofia/Desktop/all/tese/RNA-seq/raw_data_kalisto/gencode.vM14_geneInfo.RData")

#-TXI-IMPORT-####
tx2gene <- geneInfo[,c(1:2)]
inputDir <- "/Users/Ana Sofia/Desktop/all/tese/RNA-seq/raw_data_kalisto" #directory with Kallisto folders
samples <- c("kal_empty_2_1_2","kal_empty_3_1_2","kal_empty_4_1_2","kal_TV2_2_1_2","kal_TV2_3_1_2","kal_TV2_4_1_2") #vector with folder names
files <- file.path(inputDir, samples,"abundance.tsv")
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene,txOut=T,countsFromAbundance = "lengthScaledTPM")
txi_gene <- summarizeToGene(txi, tx2gene,ignoreTxVersion = T) #summaryze transcripts to genes
readCounts <- txi_gene$counts
readlenght <- txi_gene$length

head(tx)
source("https://bioconductor.org/biocLite.R")
source("https://bioconductor.org/biocLite.R")
biocLite("pasilla")

biocLite("apeglm")
biocLite("ashr")
biocLite("DESeq")
library("DESeq2")
library("pasilla")

coldata<- data.frame(condition=factor(rep(c("kal_empty", "kal_TV2"), each=3)))
rownames(sampleTable)<- colnames(txi_gene$counts)

#Create a DESeqDataset(dds) from the count table and corresponding metadata

dds2<- DESeqDataSetFromTximport(txi_gene, colData = coldata, ~condition)
colData(dds2)$condition<- factor(colData(dds2)$condition, levels = c("oe", "e"))

#The standard differential expression analysis steps are wrapped into a single function, DESeq.

dds2<- DESeq(dds2)

res<- results(dds2, name)

res<- res[order(res$padj),]
head(res)

summary(res)

vsd<- varianceStabilizingTransformation(res)
plotPCA(vsd, intgroup=c("condition"))

plotMA(dds2)
plotMA(res)

#To see how many many adjusted p-values were less than 0.1?

sum(res$padj < 0.1, na.rm=TRUE) #5847

#Results function automatically performs independent filtering based on the mean of normalized 
#counts for each gene, optimizing the number of genes which will have an adjusted p value below 
#a given FDR cutoff, alpha. By default the argument alpha is set to 0.1. To adjust pvalue cutoff, 
#alpha must be set to that value:
res05 = results(dds, alpha=0.05)

# Principal component analysis (PCA) plot of the samples

#Related to the distance matrix is the PCA plot, which shows the samples in the 2D plane spanned by their
#first two principal components (subset and replicates). This type of plot is useful for visualizing 
#the overall effect of experimental covariates and batch effects.

par(mar=c(5,10,1,2))
plotPCA(vsd, intgroup = "condition")



plotDispEsts(cds2)

res<- results(cds2)
resultsNames(cds2)


#-deseq-####

resLFC2<- lfcShrink(cds2, coef=2)
resLFC2



library("BiocParallel")
register(MulticoreParam(4))

resOrdered<- res[order(res$pvalue), ]
summary(res)

#count how many adjusted p-values are less than .1
sum(res$padj < .1, na.rm=TRUE)


res05<- results(cds2, alpha = .5)
summary(res05)

sum(res05$padj < .5, na.rm = TRUE)

# the use of IHW for p value adjustment of DESeq2 results

source("https://bioconductor.org/biocLite.R")
biocLite("IHW")
library("IHW")

resIHW<- results(cds2, filterFun = ihw)
summary(resIHW)

sum(resIHW$padj<.1, na.rm = TRUE)

metadata(resIHW)$ihwResult

plotMA(res, ylim = c(-2,2))
plotMA(resLFC, ylim = c(-2,2))

resApe<- lfcShrink(cds2, coef = 2, type = "apeglm")
resAsh<- lfcShrink(cds2, coef = 2, type = "ashr")

par(mfrow=c(1,3), mar=c(4,4,2,1))

xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="normal")
plotMA(resApe, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

plotCounts(cds2, gene=which.min(res$padj), intgroup = "condition")

d <- plotCounts(cds2, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

mcols(res)$description

write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")

resSig <- subset(resOrdered, padj < 0.1)
resSig

colData(dds2)

ddsMF <- dds2
levels(ddsMF$condition)

levels(ddsMF$condition) <- sub("-.*", "", levels(ddsMF$condition))
levels(ddsMF$condition)

design(ddsMF) <- formula(~ condition)
ddsMF <- DESeq(ddsMF)

resMF <- results(ddsMF)
head(resMF)

vsd <- vst(dds2, blind=FALSE)
rld <- rlog(dds2, blind=FALSE)
head(assay(vsd), 3)

# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


plotPCA(vsd, intgroup="condition")

pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


dds2 <- estimateSizeFactors(dds2)
dds2 <- estimateDispersions(dds2)
dds2 <- nbinomWaldTest(dds2)

results(dds2)

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds2)[["cooks"]]), range=0, las=2)

plotDispEsts(dds2)


ddsCustom <- dds2
useForMedian <- mcols(ddsCustom)$dispGeneEst > 1e-7
medianDisp <- median(mcols(ddsCustom)$dispGeneEst[useForMedian],
                     na.rm=TRUE)
dispersionFunction(ddsCustom) <- function(mu) medianDisp
ddsCustom <- estimateDispersionsMAP(ddsCustom)

metadata(res)$alpha

metadata(res)$filterThreshold

plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)

resNoFilt <- results(dds2, independentFiltering=FALSE)
addmargins(table(filtering=(res$padj < .1),
                 noFiltering=(resNoFilt$padj < .1)))


library(DESeq)

par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds2, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds2, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds2, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds2, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)

plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()

mcols(dds2,use.names=TRUE)[1:4,1:4]
substr(names(mcols(dds2)),1,10) 
mcols(mcols(dds2), use.names=TRUE)[1:4,]

head(assays(dds2)[["mu"]])
head(assays(dds2)[["cooks"]])
head(dispersions(dds2))
head(mcols(dds2)$dispersion)
sizeFactors(dds2)
head(coef(dds2))
attr(dds2, "betaPriorVar")
priorInfo(resLFC)
priorInfo(resApe)
priorInfo(resAsh)

dispersionFunction(dds2)
attr(dispersionFunction(dds2), "dispPriorVar")

normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds2) <- normFactors

W <- res$stat
maxCooks <- apply(assays(dds2)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds2)
p <- 3
abline(h=qf(.99, p, m - p))

plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

dgeFilt$samples$group <- dgeFilt$samples$condition
plotSmear(dgeFilt, de.tags = rownames(selectedFilt))
volcanoData <- cbind(resLRTfilt$table$logFC, -log10(resLRTfilt$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
DEGs <- resLRTfilt$table$FDR < 0.05 & abs(resLRTfilt$table$logFC) > 1
point.col <- ifelse(DEGs, "red", "black")
plot(volcanoData, pch = 16, col = point.col, cex = 0.5)

#--------------------------




#-edgr-####
library(edgeR)

cts <- txi_gene$counts
normMat <- txi_gene$length
normMat <- normMat/exp(rowMeans(log(normMat)))

o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)
readsf<- as.data.frame(y$offset)
names(readsf)<-c( "empty_2", "empty_3", "empty_4", "TV2_2", "TV2_3", "TV2_4")

#Filter weakly expressed and noninformative(e.g., non-aligned) features:

noint = rownames(readsf) %in% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
cpms = cpm(readsf)
keep=rowSums(cpms>1)>=3 & !noint #remove featres without at least 1 read pr million in n (3) of the samples, where n is the size of smallest group of replicates.

head(CountDF)
counts<- CountDF[,c(1:ncol(CountDF))]
head(counts)



counts=counts[keep,]

head(counts) 
tail(counts)
dim(counts)

DGE = DGEList(counts = counts, group = group)

#Estimate normalization factors using:

DGE = calcNormFactors(DGE)
DGE$samples

plotMDS(DGE, main="MDS Plot for Count Data", labels=colnames(DGE$counts))

#Estimate tagwise dispersion:
DGE = estimateCommonDisp(DGE)
DGE = estimateTagwiseDisp(DGE)

#Create a visual representation of the mean-variance relationship using the plotMeanVar and plotBCV functions:

plotMeanVar(DGE, show.tagwise.vars = TRUE, NBline = TRUE) #explore the mean-variance relationship, where each
#dot represents the estimated mean and variance of each gene, with binned variances as well as the trended 

#common dispersion overlaid.
plotBCV(DGE) #illustrated the relationship of biological coefficient of variaion vs mean log CPM.

#Test for differential expression as follows:

de = exactTest(DGE, pair = c("Empty", "TV2")) 
summary(decideTests(de, p.value=0.1))

#Use the topTags function to present a tabular summary of the differential expression
#statistics, ordered by pvalue:
tt = topTags(de, n = nrow(DGE))
head(tt$table)

#Inspect the depth-adjusted reads per milion for the 5 top differentially expressed genes:
nc = cpm(d, normalized.lib.sizes = TRUE)
rn = rownames(tt$table)

head(nc[rn,order()])

#Create a graphical summary such as an M (log-fold change) versus A (log-average expression) 
#plot here showing the genes selected as differentially expressed (with a 10% false discovery rate:

deg = rn[tt$table$FDR < .1] #para restringirmos os resultados em funçao do p-value:  substituir FDR por PValue]


plotSmear(DGE, de.tags = deg) #plots the logFC against the log counts per million.

write.csv(tt$table,file="toptags_edgeR.csv")

with(tt$table, plot(logFC, -log10(FDR), pch=20, main="tv2/Empty", xlim=c(-10,10)))
with(tt$table, plot(tt$table$logFC, -log10(tt$table$FDR), pch=20, main="tv2/Empty", xlim=c(-10,8)))


