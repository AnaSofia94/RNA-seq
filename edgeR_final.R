###################################################################################
#################################25/01/2018#######################################

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

#put the header in the table
CountDF<- as.data.frame(readCounts)

names(CountDF)<-c("empty_2", "empty_3", "empty_4", "TV2_2", "TV2_3", "TV2_4")

library("edgeR")

lengthDF<- as.data.frame(readlenght)
names(lengthDF)<-c("empty_2", "empty_3", "empty_4", "TV2_2", "TV2_3", "TV2_4")

rawCountTable<- CountDF
sampleInfo <- read.table("SampleInfo.csv", header=TRUE, sep=",", row.names = 1)

nrow(rawCountTable)

dgeFull <- DGEList(rawCountTable, remove.zeros = TRUE) #removed 19598 rows

dgeFull
dgeFull$counts

dgeFull$samples

# information about sample genotypes and conditions

dgeFull$samples$condition<- relevel(sampleInfo$Condition, ref = "wt")

dgeFull$samples$genotype<- sampleInfo$Genotype
dgeFull

#Data exploration and quality assessment

#exploratory analysis performed on log2 transformed counts to avoid problems due to skewness of the count distribution

pseudocounts<- log2(dgeFull$counts + 1)
head(pseudocounts)

#histogram for pseudo-counts

hist(pseudocounts[, "empty_2"], main="", xlab="counts")

#Boxplot for pseudo-counts
par(mar=c(8,4,1,2))
boxplot(pseudocounts, col="gray", las=3, cex.names=1)

#MA plots betweem the first two WT samples (using limma package); the other plot between two othes samples(first vs third)

limma::plotMA(pseudocounts[, 1:2], xlab="M", ylab="A", main="")
abline(h=0, col="red")

limma::plotMA(pseudocounts[, 1:3], xlab="M", ylab="A", main="")
abline(h=0, col="red")

#MA plots betweem the first two OE samples
limma::plotMA(pseudocounts[, 4:5], xlab="M", ylab="A", main="")
abline(h=0, col="red")

limma::plotMA(pseudocounts[, 4:6], xlab="M", ylab="A", main="")
abline(h=0, col="red")

#MDS for pseudo-counts
#MDS samples are similiar to PCA when Euclidean distances are used to assess the distances between samples.In limma package 
#only 500 top genes (the most variable genes across samples).

library(RColorBrewer)


library(RColorBrewer)

colConditions<- brewer.pal(3, "Set2")
colConditions<- colConditions[match(sampleInfo$Condition,
                                    levels(sampleInfo$Condition))]

pchGenotypes <- c(8,15,16)[match(sampleInfo$Condition,
                                 levels(sampleInfo$Condition))]

plotMDS(pseudocounts, pch=pchGenotypes, col=colConditions)
legend("right", lwd=2, col=brewer.pal(3, "Set2")[1:3],
       legend = levels(sampleInfo$Condition))
legend("bottomright", pch = c(8,15,16),
       legend=levels(sampleInfo$Genotype))

#heatmap for pseudo-counts (using mixOmics package)

###################################################################################
#################################28/01/2018#######################################


library(mixOmics)
library(VennDiagram)
library(HTSFilter)

sampleDists<- as.matrix(dist(t(pseudocounts)))
sampleDists

#write table samplesDist
write.csv(sampleDists,file="sampleDists.csv")

cimColor<- colorRampPalette(rev(brewer.pal(9, "Reds")))(16)
cim(sampleDists, color=cimColor, symkey=FALSE, row.cex= 0.7, col.cex= 0.7)


#Normalization
#Compute normalization factors

dgeFull<- calcNormFactors(dgeFull, method = "TMM")
dgeFull

write.csv(dgeFull$counts,file="dgeFull_counts.csv")

write.csv(dgeFull$samples,file="dgeFull_samples.csv")

tail(dgeFull$counts)


#Normalized counts exploratory analysis

normCounts<- cpm(dgeFull)
pseudoNormCounts<- cpm(dgeFull, log = TRUE, prior.count = 1)
par(mar=c(8,4,1,2))
boxplot(pseudoNormCounts, col="red", las=3, cex.names=2)

plotMDS(pseudoNormCounts, pch = pchGenotypes, col= colConditions)
legend("right", lwd=2, col = brewer.pal(3, "Set2")[1:2],
       legend = levels(sampleInfo$Condition))
legend("bottomright", pch=c(8,15,16), 
       legend = levels(sampleInfo$Genotype))


#Differential Analysis-####

#compare the results of different types of approaches to obtain genes which are differentially expressed between the 
#wild type tomatoes and the mutants:
#standard NB exact test between two conditions;
#GLM with the plant and genotype effects 

#differences between the WT and M is tested using NB test. using DGEList with group and using the same normalization 
#factors than in dgeFull;

#estimating the dispersion for this object with the functions estimateCommonDisp and estimateTagwiseDisp

#performing the test with the function exactTest

dgeFull.group<- DGEList(rawCountTable, remove.zeros = TRUE, group = dgeFull$samples$condition) #Removing 19598 rows with all zero counts

dgeFull.group$samples$norm.factors<- dgeFull$samples$norm.factors
dgeFull.group

#Estimate dispersion

dgeFull.group<- estimateCommonDisp(dgeFull.group)
dgeFull.group<- estimateCommonDisp(dgeFull.group)
dgeFull.group

#dispersion in a more robust way, that can be used if the previous approach seems to fail.This function 
#is estimateGLMRobustDisp. The quality of the variety estimation can be assessed with the BCV versus average 
#log CPM plot(that plots 0 versus the average normalized count for all genes)

plotBCV(dgeFull.group)

#Perform the test

dgeExactTest<- exactTest(dgeFull.group)
dgeExactTest

#p-values are corrected with the function topTags

resExactTest<- topTags(dgeExactTest, n=nrow(dgeExactTest$table))
head(resExactTest$table)

#p-value and (BH) adjusted p-value distribution can be assessed with:

par(mfrow=c(1,2))
hist(resExactTest$table$PValue, xlab = "p-value", main = "raw p-values")
hist(resExactTest$table$FDR, xlab="p-value", main = "adjusted p-values")

#genes with a FDR smaller than 5% and a log Fold Change larger than 1 or smaller than -1 are extracted

selectedET<- resExactTest$table$FDR < 0.05 & abs(resExactTest$table$logFC) > 1

selectedET<- resExactTest$table[selectedET, ]
nrow(selectedET) #1807

head(selectedET)

#these method shows that 1807 genes are found differential with these method. The column logFC can be used to found
#up/down regulated genes in the M

selectedET$upDown<- factor(ifelse(selectedET$logFC>0, "up", "down"))
head(selectedET)

write.table(selectedET, file="tomatoDEG.csv", sep = ",")

#Second approach: GLM with condition and genotype effects

#estimate dispersion 

design.matrix<- model.matrix(~ dgeFull$samples$condition +
                               dgeFull$samples$genotype)

design.matrix

#Common, trended and then tagwise dispersions can be estimated with:

dgeFull<- estimateDisp(dgeFull, design.matrix)
dgeFull

#the qualitiy of the variability estimation can be assessed with the BCV versus average log CPM plot ( that plots
#0 versus the average normalized count for all genes):

plotBCV(dgeFull)

#fil GLM and perform the test

fit<- glmFit(dgeFull, design.matrix)
fit

#test performed with a log-ratio test. To test differential genes between WT and M is equivalent to testing the
#nullity of the second coefficient

dgeLRTtest<- glmLRT(fit, coef = 2)
dgeLRTtest

#testing differential genes between clone number 2 and 3 is equivalent to testing the equality of coefficients 3 and 4:

contrasts<- rep(0, ncol(design.matrix))
contrasts[3]<- 1
contrasts[4]<- 1
dgeLRTtest2<- glmLRT(fit, contrast = contrasts)
dgeLRTtest2

#extracted as previously using the the function topTags

resLRT<- topTags(dgeLRTtest, n=nrow(dgeFull$counts))
head(resLRT$table)

selectedLRT<- resLRT$table$FDR < 0.05 & abs(resLRT$table$logFC) > 1
selecterLRT<- resLRT$table[selectedLRT, ]
nrow(selecterLRT) #1838

head(selecterLRT)


#comparison
#venn diagram comparing the two approaches

vd<- venn.diagram(x=list("Exact test"= rownames(selectedET),
                         "GLM"= rownames(selectedLRT)),
                  fill= brewer.pal(3, "Set2")[1:2], filename = NULL)
grid.draw(vd)


#filtering;Differential analysis after independent filtering

source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

library(HTSFilter)

dgeFilt<- HTSFilter(dgeFull)$filteredData
dgeFilt

fit<- glmFit(dgeFilt, design.matrix)
dgeLRTfilt<- glmLRT(fit, coef = 2)
resLRTfilt<- topTags(dgeLRTfilt, n=nrow(dgeFilt$counts))
selectedFilt<- resLRTfilt$table$FDR< 0.05 & abs(resLRTfilt$table$logFC)>1
selectedFilt<- resLRTfilt$table[selectedFilt, ]
nrow(selectedFilt)

head(selectedFilt)

#Comparison

vd<- venn.diagram(x=list("No filtering"= rownames(selectedLRT),
                         "Filtering"= rownames(selectedFilt)),
                  fill=brewer.pal(3, "Set2")[1:2], filename = NULL)

#MA plot with differentially expressed genes 
#to create a MA plot between oe and m, the entry $sample$group of the DGEList object must be filled with 
#the indication of what two groups are

dgeFilt$samples$group<- dgeFilt$samples$condition
plotSmear(dgeFilt, de.tags = row.names(selectedFilt))

#Volcano Plot

volcanoData<- cbind(resLRTfilt$table$logFC, -log10(resLRTfilt$table$FDR))
colnames(volcanoData)<- c("logFC", "negLogPval")
DEGs<- resLRTfilt$table$FDR < 0.05 & abs(resLRTfilt$table$logFC)>1
point.col<- ifelse(DEGs, "red", "black")
plot(volcanoData, pch=16, col=point.col, cex=0.5)

#heatmap

selY<- cpm(dgeFilt, log = TRUE, prior.count = 1)
selY<- selY[match(rownames(selectedFilt), rownames(dgeFilt$counts)),]
finalHM<- cim(t(selY), color = cimColor, symkey = FALSE, row.cex = 0.7,
              col.cex = 0.1)

plot(finalHM$ddc, leaflab = "none")
abline(h=8, lwd=2, col="pink")

geneClust<- cutree(as.hclust(finalHM$ddc), h=10)
head(geneClust)

length(unique(geneClust))

names(which(geneClust==1))


#-another aproach-####
###################################################################################
#################################5/02/2018#######################################



#data frame from the counts of tximport

CountDF<- as.data.frame(readCounts)

names(CountDF)<-c("empty_2", "empty_3", "empty_4", "TV2_2", "TV2_3", "TV2_4")

#starting edger
cts <- txi_gene$counts
normMat <- txi_gene$length
normMat <- normMat/exp(rowMeans(log(normMat)))

o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)
readsf<- as.data.frame(y$offset)
names(readsf)<-c( "empty_2", "empty_3", "empty_4", "TV2_2", "TV2_3", "TV2_4")

#Filter weakly expressed and noninformative(e.g., non-aligned) features:

head(CountDF)
counts<- CountDF[,c(1:ncol(CountDF))]
head(counts)

dim( counts )
colSums( counts ) # Library Sizes

colSums( counts ) / 1e06 # Library Sizes in millions of reads
table( rowSums( counts ) )[ 1:30 ] # Number of genes with low counts

#Creates a dgeList

group<- c(rep("Empty",3), rep("TV2",3))
cds<- DGEList(counts, group = group)
names(cds)

head(cds$counts) #original count matrix
cds$samples #summary of my samples
sum(cds$all.zeros) # how many genes have 0 counts across all genes

#filter reads with low count 
cds<- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) >1) >=3,]
dim(cds)

cds<- calcNormFactors(cds)
cds$samples


#effective library sizes
cds$samples$lib.size * cds$samples$norm.factors

#Multi-Dimension Plot

#view the plot

plotMDS(cds, main="MDS Plot for Count Data", labels=colnames(cds$counts))

#Outplut plot as a pdf
pdf("MDS_plot_.pdf", width = 7, height = 7)
plotMDS(cds, main ="MDS Plot for count Data", labels= colnames(cds$counts) )
dev.off()

#estimating dispersion 
cds <- estimateCommonDisp(cds)
names(cds)

#the estimate
cds$common.dispersion

#parameterization for the variance of the negative binomial
sqrt(200) #poisson sd
sqrt(200 +200^2 * cds$common.dispersion) #negative binomial sd
sqrt(200 + 200^2 * cds$common.dispersion)/sqrt(200) #MB sd is over to times larger

##default settings
cds <- estimateTagwiseDisp(cds, prior.df=10)
names(cds)
summary(cds$tagwise.dispersion)

#More shrinkage/sqeezing toward the commmon
cds <- estimateTagwiseDisp(cds, prior.df = 25)
summary(cds$tagwise.dispersion)

#the recommended setting for this data is the defaul of 10. 
cds<- estimateTagwiseDisp(cds, prior.df = 10)

#Mean-Variance Plot
#Create a visual representation of the mean-variance relationship using the plotMeanVar and plotBCV functions:
meanVarPlot<- plotMeanVar(cds, show.raw.vars = TRUE,
                          show.tagwise.vars = TRUE,
                          show.binned.common.disp.vars = FALSE,
                          dispersion.method= "qml", NBline=TRUE,
                          nbins=100,
                          pch=16, 
                          xlab="Mean Expression (Log10 Scale)",
                          ylab="Variance (Log10 Scale)",
                          main= "Mean-Variance Plot")

plotBCV(cds, xlab = "Average log CPM", ylab = "Variance",
        pch=16,cex=0.2, col.common = "red", col.trend = "blue", col.tagwise = "black") 

plotBCV(cds)


cds$samples
cds$samples$group

#pairwise comparison between groups
et<- exactTest(cds, pair = c("Empty", "TV2"))
topTags(et)

y<-cds
y$samples$group<- relevel (y$samples$group, ref= "TV2")
levels(y$samples$group)

#Test for differential expression as follows(OE/E):

de<- exactTest(cds, pair = c("TV2", "Empty"))
summary(decideTests(de, p.value = 0.1))

#Use the topTags function to present a tabular summary of the differential 
#expression
#statistics, ordered by pvalue:

tt<- topTags(de, n=nrow(cds))
head(tt$table)

#Inspect the depth-adjusted reads per milion for the 5 differentially 
#expressed genes:
nc=cpm(cds, normalized.lib.sizes=TRUE)
rn<- rownames(tt$table)
head(nc[rn, order(sampleInfo$Genotype)])

#create a graphical summary such as an M (log-fold change) vs A (log-average expression)
#plot here showing the genes selected as differentially expressed (with a 10% false discovery rate):

deg=rn[tt$table$FDR< .1]
plotSmear(cds, de.tags = deg) #plots the logFC against the log counts per million.


write.csv(tt$table,file="toptags_edgeR.csv")

#To convert the Ensemble ID to gene name:

ensemble.all<- read.table("gencode.v27.annotation.gtf", sep = "\t", header = F)

ensembl.id.names= data.frame(sapply(strsplit(sapply(strsplit(as.character(ensemble.all[,9]), "gene_name"), "[", 2), ";"), "[", 1))



