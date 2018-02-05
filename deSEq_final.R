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

sampleTable<- data.frame(condition=factor(rep(c("kal_empty", "kal_TV2"), each=3)))
rownames(sampleTable)<- colnames(txi$counts)

#Create a DESeqDataset(dds) from the count table and corresponding metadata

library("DESeq2")
biocLite("DESeq2")

colData<- data.frame(condition=factor(rep(c("kal_empty", "kal_TV2"), each=3)))
rownames(colData)<- colnames(txi_gene$counts)

dds2<- DESeqDataSetFromTximport(txi_gene, colData = colData, ~condition)


#use only if the dds reveals NAs in the condition 
#dds <- dds[ , dds$condition %in% c("empty2","empty3", "empty4", "TV2_2", "TV2_3", "TV2_4")]

#filter to remove rows with 0 oe 1 read 

keep<- rowSums(counts(dds2))>=10
dds2<- dds2[keep, ]

dds<- dds[rowSums(counts(dds))>1, ]

#remove levels which do not have samples in the current DataSeq
dds$condition<- droplevels(dds$condition)


##The standard differential expression analysis steps are wrapped into a single function, DESeq
dds<- DESeq(dds2)
dds
dds$condition

mcols(dds, use.names = TRUE)

##Results tables are generated using the function results, which extracts a results table with 
#log2 fold changes, p values and adjusted p values. The comparison can be ordered using the 
#contrast arguments of results. In this case, logFC=(OE/EMPT)

res<- results(dds) #just to have and idea of the data
res
summary(res)

#specifie the coefficient or contrast to build results table

res<- results(dds, name = "condition_kal_TV2_vs_kal_empty")

res
resultsNames(dds)

#Log fold change shrinkage for ranking and visualization 
#moderation of log2 fold change split into a separate function 

resLFC<- lfcShrink(dds, coef="condition_kal_TV2_vs_kal_empty")
resLFC

#order results table by the smallest adjusted p value:

res0ordered<- res[order(res$pvalue), ]

summary(res0ordered) 

#see how many adjusted pvalues are less than 0.1
sum(res$padj< 0.5, na.rm = TRUE) #10558

#function results automatically performs independent filtering based on the mean of normalized counts for
#each gene, optimizing the number of genes which will have an adjusted p value bellow a given FDR cutoff, alpha
#by default argument alpha is set to 0.1. 

res0.5<- results(dds, alpha = 0.05)
summary(res0.5)

sum(res0.5$padj< 0.05, na.rm = TRUE) #4952

#use of IHW for p value adjustment of DESeq2 results

library("IHW")

RESIHW<- results(dds, filterFun = ihw)
summary(RESIHW)

sum(RESIHW$padj< 0.1, na.rm = TRUE) #5940

metadata(RESIHW)$ihwResult

#MA plot: log2 fold changes attributable to a given variable over the mean of normalized counts. 
#Points will be colored red if the adjusted p value is less than 0.1.

View(res)
plotMA(res, ylim=c(-2,2)) #nao da correto


#more usefull visualizing the MA-plot for the shrunken log2 fold change, which removes the noise associated
#with log2 folds changes from low count genes without requiring arbitrary filtering thresholds

plotMA(resLFC, ylim=c(-2,2))

#after calling plotMA, one can use the function identify to interactively detect the row number of individual
#genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices

idx<- identify(res$baseMean, res$log2FoldChange)

#alternative shrinkage estimators 

resultsNames(dds)

resApe<- lfcShrink(dds, coef = 2, type = "apeglm")
resAsh<- lfcShrink(dds, coef = 2, type = "ashr")

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim<- c(1,1e5); ylim<- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="normal")
plotMA(resApe, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

#Plot Counts

#count the reads for a single gene across the groups. Normalize the counts by sequencing the depth 
#and adds a pseudocount of 1/2 to allow for log scale plotting. The counts are grouped by the variables
#in intgroup, where more than one can be specified. 
#Here we specifie the gene which had the smallest p value from the results table created above. 

par(mar=c(5.1,3.1,2.1,2.1))

plotCounts(dds, gene=which.min(res$padj), intgroup = "condition")

#costumizing plotting

d<- plotCounts(dds, gene=which.min(res$padj), intgroup = "condition", 
               returnData = TRUE)

library("ggplot2")


par(mar=c(4,2,1,1))

ggplot(d, aes(x = condition, y=count)) +
  geom_point(position = position_jitter(w=0.1, h=0)) +
  scale_y_log10(breaks=c(25,100,400))

#more information on results columns: 

mcols(res)$description

#log2 fol change for -1 for condition kal TV2 vs Kal empty means that the treatment induces a multiplicative
#change in observed gene expression level of 2^-1=0.5 compared to the untretaed condition. 

write.csv(as.data.frame(res0ordered), file="condition_kal_TV2_vs_kal_empty")

resSig<- subset(res0ordered, padj< 0.5)
resSig

#Multi factor design 

colData(dds)

#create a copy of deseqdataset, to run a multi-factor design

ddsMF<- dds

#change the levels of type so it only contains letters

levels(ddsMF$condition)

#account the different types of sequencing and a get clearer picture of the differences attribute to the treatment
#as condition is a variable of interest, we put it at the end of the formula. 

design(ddsMF)<- formula(~ type + condition)

ddsMF<- DESeq(ddsMF)

#results of the function

resMF<- results(ddsMF)
head(resMF)

#data transformation and visualization
#test for de, we operato on raw counts and use discrete distribution, but to work with clustering or visualiing
#it might be useful to transform versions of the count data

vsd<- vst(dds, blind = FALSE)
rld<- rlog(dds, blind = FALSE)
head(assay(vsd),3)


#effeects of transformation on the variance

#plot the standard deviation of the transformed data, across samples, against means 

#this gives log2(n+1)

ntd<- normTransform(dds)
head(assay(ntd))

library("vsn")

meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
meanSdPlot(assay(ntd), main="ntd" )
meanSdPlot(assay(vsd), main="vsd" )
meanSdPlot(assay(rld), main="rld" )

#data qualitiy assessment by sample clustering and visualization

#heatmap of the count matrix

library("pheatmap")
select<- order(rowMeans(counts(dds, normalized=TRUE)),
               decreasing = TRUE)[1:20]

df<- as.data.frame(colData(dds)[, c("condition")])

pheatmap(assay(ntd)[select,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols=FALSE, annotation_col = df)

pheatmap(assay(vsd)[select,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols=FALSE, annotation_col = df)

pheatmap(assay(rld)[select,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols=FALSE, annotation_col = df)


#heatmap of the sample-to-sample distances

#dist function to the transpose of the transformed count matrix to get sample-to-sample distances

sampleDists<- dist(t(assay(vsd)))

#heatmap gives us an overview over similarities and dissimilarities between samples. 

library("RcolorBrewer")

sampleDistMatrix<- as.matrix(sampleDists)

rownames(sampleDistMatrix)<- paste(vsd$condition)
colnames(sampleDistMatrix)<- NULL
colors<- colorRampPalette(rev(brewer.pal(9, "Blues")) ) (225)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors)

#Principal component plot of the samples

#shows the samples in the 2D plane spanned by their first two principal components
#useful for visualizing the overall effect of experimental covariates and batch effects

plotPCA(vsd, intgroup="condition")

#other way

pcaData<- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar<- round(100* attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color="condition")) + 
  geom_point(size=3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()


#variations to the standard workflow

dds<- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- nbinomWaldTest(dds)


#Interactios:can be added to the desing formula in order to test, for example log2 fold change attributable
#to a given condition is different based on another factor

dds$group<- factor(paste0(dds$condition))
design(dds)<- ~ group
dds<- DESeq(dds)
resultsNames(dds) #group_kal_TV2_vs_kal_empty

results(dds, contrast =c("group", "TV2", "empty"))

# Cook's distance is a measure of how much a single sample is influencing the fitted coefficients 
#for a gene, and a large value of Cook's distance is intended to indicate an outlier count.

par(mar=c(8,5,2,2))
boxplot(log2(assay(dds)[["cooks"]]), range=-1, las=2)

#dispersion plot and fitting alternatives
#final estimates shrunk from the gene-wise etimates towards the fitted estimates. 

plotDispEsts(dds)

#Tests log2 fold change above or bellow a threshold

par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim<- c(-2.5, 2.5)
resGA<- results(dds, lfcThreshold = .5, altHypothesis = "greaterAbs")
resLA<- results(dds, lfcThreshold = .5, altHypothesis = "lessAbs")
resG<- results(dds, lfcThreshold = .5, altHypothesis = "greater")
resL<- results(dds, lfcThreshold = .5, altHypothesis = "less")
drawlines<- function() abline(h=c(-.5,.5), col="dodgerblue", lwd=2)
plotMA(resGA, ylim=ylim); drawlines()
plotMA(resLA, ylim=ylim); drawlines()
plotMA(resG, ylim=ylim); drawlines()
plotMA(resL, ylim=ylim); drawlines()


#access to all calculated values 

mcols(dds, use.names = TRUE) [1:4, 1:4]

substr(names(mcols(dds)),1,10)

mcols(mcols(dds), use.names = TRUE) [1:4, ]

head(assays(dds)[["mu"]])

head(assays(dds)[["cooks"]])

head(dispersions(dds))

head(mcols(dds)$dispersion)

sizeFactors(dds)

head(coef(dds))

attr(dds, "betaPriorVar")

priorInfo(resLFC)


#-----------------------------------another PLOT-####
library(DESeq2)


ddsm<- DESeq(dds)



res2<- results(ddsm, name = "group_kal_TV2_vs_kal_empty")

summary(res) #para ver summary com outro pvalue: summary(res,alpha=x)

#To see how many many adjusted p-values were less than 0.1?
sum(res2$padj < 0.1, na.rm=TRUE)
#5844

#Results function automatically performs independent filtering based on the mean of normalized 
#counts for each gene, optimizing the number of genes which will have an adjusted p value below 
#a given FDR cutoff, alpha. By default the argument alpha is set to 0.1. To adjust pvalue cutoff, 
#alpha must be set to that value:

res05 = results(dds, name = "group_kal_TV2_vs_kal_empty", alpha=0.05)


#sample clustering and visualization

#Another use of the transformed data is sample clustering. Here, dist function is applied to  
#transpose the transformed count matrix to get sample-to-sample distances, which requires vsd(variance-stabilizing transformation)
#A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between
#samples:

vsd = varianceStabilizingTransformation(ddsm)

sampleDists = dist(t(assay(vsd)))
library("RColorBrewer")

sampleDistMatrix = as.matrix(sampleDists)

rownames(sampleDistMatrix) =colnames(ddsm)
colnames(sampleDistMatrix) = colnames(ddsm)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
library(pheatmap)
pheatmap (sampleDistMatrix, 
          clustering_distance_rows=sampleDists,
          clustering_distance_cols=sampleDists, 
          col=colors) 

# Principal component analysis (PCA) plot of the samples

#MA-plot shows the log2 fold changes attributable to a given variable over the mean of normalized counts 
#for all samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. 
#Points which fall out of the window are plotted as open triangles pointing either up or down.

plotMA(res05, ylim=c(-10,10), main="Empty vs. OVerExpressed")


plot(res$baseMean, res$log2FoldChange,
     main="MA plot",
     xlab= "Mean of normalized counts",
     ylab = "Log2FC (OVerExpressed/Empty )",
     col=densCols(res$log2FoldChange, res$baseMean),
     grid(lty="solid",col="grey"),
     pch=20)
abline(col="black", h=0)
abline(col="red", h=c(1, -1))

#MA-plot for the shrunken log2 fold changes: Remove the noise associated with log2 fold changes from 
#low count genes without requiring arbitrary filtering thresholds.

resLFC = lfcShrink(dds, coef=2)
plotMA(resLFC, ylim=c(-6,6), main="Empty vs. OVerExpressed")

#Plotting the dispersion estimates is a useful diagnostic. The dispersion plot below is 
#typical, with the final estimates shrunk from the gene-wise estimates towards the fitted 
#estimates. Some gene-wise estimates are flagged as outliers and not shrunk towards the 
#fitted value. The amount of shrinkage can be more or less than seen here, depending on the 
#sample size, the number of coefficients, the row mean and the variability of the gene-wise estimates.

plotDispEsts(dds)



#Examine the counts of reads for a single gene across the groups: plotCounts normalizes counts by sequencing 
#depth and adds a pseudocount of 1/2 to allow for log scale plotting. The counts are grouped by the 
#variables in intgroup, where more than one variable can be specified. which.min specify the gene which 
#had the smallest p adj value from the results table (res) created above. You can select the gene to plot 
#by rowname or by numeric index.

plotCounts(dds, gene=which.min(res$padj), intgroup="", main="")

