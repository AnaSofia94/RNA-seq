####################################################################################################
#Cell-line-#### 
# BaF3 ==> paired end 
# Samples:
#empty_2_1.fq
#empty_2_2.fq
#empty_3_1.fq
#empty_3_2.fq
#empty_4_1.fq
#empty_4_2.fq
#TV2_2_1.fq
#TV2_2_2.fq
#TV2_3_1.fq
#TV2_3_2.fq
#TV2_4_1.fq
#TV2_4_2.fq

#UCSC Genome Browser assembly ID: mm10
#Sequencing/Assembly provider ID: Genome Reference Consortium Mouse Build 38 (GCA_000001635.2)
#Assembly date: Dec. 2011
#Accession ID: GCA_000001305.2
#NCBI Genome ID: 52 (Mus musculus)
#NCBI Assembly ID: 327618 (GRCm38/GCA_000001635.2)
#NCBI BioProject ID: 20689

#Download the total files from the sequencing company-####
wget --ftp-user=20171106F17FTSEUHT1267 ftp://20171106F17FTSEUHT1267:MOUgmkE@cdts-wh.genomics.cn/F17FTSEUHT1267_MOUgmkE/*
wget --ftp-user=20171106F17FTSEUHT1267 ftp://20171106F17FTSEUHT1267:MOUgmkE@cdts-wh.genomics.cn/F17FTSEUHT1267_MOUgmkE/IGV/bam/*  
  
  #confirm if the download is correct-####
#empty_2_1.fq:md5sum -c "83ddf74428afa89fd72ef22bbff60436 latest/Clean_Data/empty_2_1.fq.gz"
#empty_2_2.fq:md5sum -c "0cd9441867a13fbc152aee751d54462b latest/Clean_Data/empty_2_2.fq.gz"
#empty_3_1.fq:md5sum -c "94f9ce97b83ab61a49c3a67d3bc8b8d9 latest/Clean_Data/empty_3_1.fq.gz"
#empty_3_2.fq:md5sum -c "be7926dea3091b98e5fd1bda35bae1bb latest/Clean_Data/empty_3_2.fq.gz"
#empty_4_1.fq:md5sum -c "9ff079987cedb9551151d0df3eb5196d latest/Clean_Data/empty_4_1.fq.gz"
#empty_4_2.fq:md5sum -c "3061f1e27e78c215704f51b4fb43d5a2 latest/Clean_Data/empty_4_2.fq.gz"
###TV2_2_1.fq:md5sum -c "ad12df71ca5fb934b724a9b074fad863 latest/Clean_Data/TV2_2_1.fq.gz"
###TV2_2_2.fq:md5sum -c "1ebe68d34f26b680cab564bddaeb13e0 latest/Clean_Data/TV2_2_2.fq.gz"
###TV2_3_1.fq:md5sum -c "6942e690a59998a81362ab49ce5e62d9 latest/Clean_Data/TV2_3_1.fq.gz"
###TV2_3_2.fq:md5sum -c "23903faa955dc80fe47e8836e400e1d4 latest/Clean_Data/TV2_3_2.fq.gz"
###TV2_4_1.fq:md5sum -c "17dc9300efb4d2fd2965493a3c1dc803 latest/Clean_Data/TV2_4_1.fq.gz"
###TV2_4_2.fq:md5sum -c "6f2fd3aded15ec908bea05558c00ff49 latest/Clean_Data/TV2_4_2.fq.gz"

####################################################################################################
###############################20/11/2017###############################################

#Fastq to samples with parallel-####

#--------------------------------#
#------------- SBATCH File Step 1 -  FASTQ and check quality ####
#--------------------------------#

#!/bin/bash
#SBATCH --job-name=Par_fastq
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=4
#SBATCH --nodeslist=compute-13 
#SBATCH --workdir=/mnt/beegfs/scratch/JBARATA/anasofiamoreira/MasterAnaS/latest

### ----------------------------------------------------------- INPUTS --------------------------------------------- ###

## working directory
export wk=/mnt/beegfs/scratch/JBARATA/anasofiamoreira/MasterAnaS/latest

##export directories in clean_data
export Clean_Data=$wk/Clean_data/
  export fastqc =$wk/Clean_Data/fastqc/
  
  
  srun="srun --nodes=1 --ntasks=1 --cpus=$SLURM_CPUS_PER_TASK"
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog $wk/$SLURM_JOB_ID.log"

## Check the quality of the reads with fastqc
ls $Clean_Data/*.fq | $parallel '$srun shifter --image=argrosso/htspreprocessing:0.1.1 fastqc -o $fastqc_dir -t $SLURM_CPUS_PER_TASK {}'

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


#-Again-Pipe-####
#!/bin/bash
#SBATCH --job-name=Par_fastq
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=4
#SBATCH --nodelist=compute-13
#SBATCH --workdir=/mnt/beegfs/scratch/JBARATA/anasofiamoreira/MasterAnaS/latest

## working directory
export wk=/mnt/beegfs/scratch/JBARATA/anasofiamoreira/MasterAnaS/latest
export LANG="" #deixa de dar warning perl

##export directories in clean_data
export Clean_Data=$wk/Clean_Data
export fastqc=$wk/Clean_Data/fastqc


export srun="srun --nodes=1 --ntasks=1 --cpus=$SLURM_CPUS_PER_TASK"
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog $wk/$SLURM_JOB_ID.log"

## Check the quality of the reads with fastqc
time ls $Clean_Data/*.fq | $parallel '$srun shifter --image=argrosso/htspreprocessing:0.1.1' 

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitc$"

#-kalisto_Alignment-####

#!/bin/bash
#SBATCH --job-name=Par_kalisto
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=4
#SBATCH --nodelist=compute-13
#SBATCH --workdir=/mnt/beegfs/scratch/JBARATA/anasofiamoreira/MasterAnaS/latest

## working directory
export wk=/mnt/beegfs/scratch/JBARATA/anasofiamoreira/MasterAnaS/latest
export LANG="" #deixa de dar warning perl

export TRANSCRIPTOME=/mnt/nfs/lobo/IMM-NFS/genomes/mm10/Sequence/KallistoIndex/transcriptome.idx

##export directories in clean_data
export Clean_Data=$wk/Clean_Data
export fastqc=$wk/Clean_Data/fastqc
export kalisto=$wc/Clean_Data/kalisto


export srun="srun --nodes=1 --ntasks=1 --cpus=$SLURM_CPUS_PER_TASK"
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog $wk/$SLURM_JOB_ID.log"

## Check the quality of the reads with fastqc
time ls $Clean_Data/*.fq | $parallel '$srun shifter --image=docker:argrosso/kallisto:latest  kallisto quant -i $TRANSCRIPTOME -o $kalisto $Clean_Data/*.fq -t $SLURM_CPUS_PER_TASK {}' 

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitc$"


export wk=/mnt/beegfs/scratch/JBARATA/anasofiamoreira/MasterAnaS/latest
export LANG="" #deixa de dar warning perl

#-Code right Par (in the end)-####
##export directories in clean_data
export Clean_Data=$wk/Clean_Data
export fastqc=$wk/Clean_Data/fastqc


export srun="srun --nodes=1 --ntasks=1 --cpus=$SLURM_CPUS_PER_TASK"
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog $wk/$SLURM_JOB_ID.log"

## Check the quality of the reads with fastqc
time ls $Clean_Data/*.fq | $parallel '$srun shifter --image=argrosso/htspreprocessing:0.1.1 fastqc -o $fastqc/ -t $SLURM_CPUS_PER_TASK {}'

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

####################################################################################################
###############################22/11/2017###############################################

#-Par_kalisto-####
#!/bin/bash
#SBATCH --job-name=Par_kalisto
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=4
#SBATCH --nodelist=compute-13
#SBATCH --workdir=/mnt/beegfs/scratch/JBARATA/anasofiamoreira/MasterAnaS/latest

## working directory
export wk=/mnt/beegfs/scratch/JBARATA/anasofiamoreira/MasterAnaS/latest
export LANG="" #deixa de dar warning perl

export TRANSCRIPTOME=/mnt/nfs/lobo/IMM-NFS/genomes/mm10/Sequence/KallistoIndex/transcriptome.idx

##export directories in clean_data
export Clean_Data=$wk/Clean_Data
export fastqc=$wk/Clean_Data/fastqc
export kalisto=$wc/Clean_Data/kalisto


export srun="srun --nodes=1 --ntasks=1 --cpus=$SLURM_CPUS_PER_TASK"
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog $wk/$SLURM_JOB_ID.log"

## Check the quality of the reads with fastqc
time ls $Clean_Data/*.fq | $parallel '$srun shifter --image=docker:argrosso/kallisto:latest  kallisto quant -i $TRANSCRIPTOME -o $kalisto $Clean_Data/*.fq -t $SLURM_CPUS_PER_TASK {}' 

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitc$"

#-kalisto-####

#!/bin/bash
#SBATCH --job-name=kal_rna
#SBATCH --time=4:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/kallisto:latest
#SBATCH --workdir=/mnt/beegfs/scratch/JBARATA/anasofiamoreira/MasterAnaS/latest/Clean_Data

export TRANSCRIPTOME=/mnt/nfs/lobo/IMM-NFS/genomes/mm10/Sequence/KallistoIndex/transcriptome.idx

export kalisto= $wk/kalisto

srun shifter kallisto quant -i $TRANSCRIPTOME -o kal_empty_2_1_2/$kalisto -t 10 -b 100 empty_2_1.fq empty_2_2.fq
srun shifter kallisto quant -i $TRANSCRIPTOME -o kal_empty_3_1_2/$kalisto -t 10 -b 100 empty_3_1.fq empty_3_2.fq
srun shifter kallisto quant -i $TRANSCRIPTOME -o kal_empty_4_1_2/$kalisto -t 10 -b 100 empty_4_1.fq empty_4_2.fq
srun shifter kallisto quant -i $TRANSCRIPTOME -o kal_TV2_2_1_2/$kalisto -t 10 -b 100 TV2_2_1.fq TV2_2_2.fq
srun shifter kallisto quant -i $TRANSCRIPTOME -o kal_TV2_3_1_2/$kalisto -t 10 -b 100 TV2_3_1.fq TV2_3_2.fq
srun shifter kallisto quant -i $TRANSCRIPTOME -o kal_TV2_4_1_2/$kalisto -t 10 -b 100 TV2_4_1.fq TV2_4_2.fq


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derive$"

####################################################################################################
###############################23/11/2017###############################################

#upload the data from the mobExterm

empty_2_1_2<- read.table("kal_empty_2_1_2.tsv",header = TRUE);
empty_3_1_2<- read.table("kal_empty_3_1_2.tsv",header = TRUE);
empty_4_1_2<- read.table("kal_empty_4_1_2.tsv",header = TRUE);

TV2_2_1_2<- read.table("kal_TV2_2_1_2.tsv",header = TRUE);
TV2_3_1_2<- read.table("kal_TV2_3_1_2.tsv",header = TRUE);
TV2_4_1_2<- read.table("kal_TV2_4_1_2.tsv",header = TRUE);

#Prepare the data table

#-----------empty data 

empty_2_1_2 <- as.matrix(read.delim("kal_empty_2_1_2.tsv", sep= "\t", header=T))

ID_empty_2_1_2 <- lapply(strsplit(empty_2_1_2[,1],"\\|"),function(x) x[1:2])

empty_2_1_2[,1] <- as.vector(unlist(lapply(ID_empty_2_1_2,function(x) x[1])))
empty_2_1_2 <- cbind(empty_2_1_2,"Gene_ID"=c(as.vector(unlist(lapply(ID_empty_2_1_2,function(x) x[2])))))

GeneID <- unique(empty_2_1_2[,"Gene_ID"])

Gene_Em2 <- matrix(data=NA,ncol=2)

for(a in 1:length(GeneID)){
  curGene <- GeneID[a]
  curTable <- matrix(empty_2_1_2[which(empty_2_1_2[,"Gene_ID"]==curGene),],ncol = ncol(empty_2_1_2)); colnames(curTable) <- colnames(empty_2_1_2)
  TPM <- max(as.numeric(curTable[,"tpm"]))
  Gene_Em2 <- rbind(Gene_Em2,c(curGene,TPM))
}

#--------------------------------------------

empty_3_1_2 <- as.matrix(read.delim("kal_empty_3_1_2.tsv",sep="\t",header=T))

ID <- lapply(strsplit(empty_3_1_2[,1],"\\|"),function(x) x[1:2])

empty_3_1_2[,1] <- as.vector(unlist(lapply(ID,function(x) x[1])))
empty_3_1_2 <- cbind(empty_3_1_2,"Gene_ID"=c(as.vector(unlist(lapply(ID,function(x) x[2])))))

GeneID <- unique(empty_3_1_2[,"Gene_ID"])

Gene_Em3 <- matrix(data=NA,ncol=2)

for(a in 1:length(GeneID)){
  curGene <- GeneID[a]
  curTable <- matrix(empty_3_1_2[which(empty_3_1_2[,"Gene_ID"]==curGene),],ncol = ncol(empty_3_1_2)); colnames(curTable) <- colnames(empty_3_1_2)
  TPM <- max(as.numeric(curTable[,"tpm"]))
  Gene_Em3 <- rbind(Gene_Em3,c(curGene,TPM))
}

####################################################################################################
###############################3/12/2017###############################################

#----------empty

empty_4_1_2 <- as.matrix(read.delim("kal_empty_4_1_2.tsv",sep="\t",header=T))

ID <- lapply(strsplit(empty_4_1_2[,1],"\\|"),function(x) x[1:2])

empty_4_1_2[,1] <- as.vector(unlist(lapply(ID,function(x) x[1])))
empty_4_1_2 <- cbind(empty_4_1_2,"Gene_ID"=c(as.vector(unlist(lapply(ID,function(x) x[2])))))

GeneID <- unique(empty_4_1_2[,"Gene_ID"])

GeneEmp4 <- matrix(data=NA,ncol=2)

for(a in 1:length(GeneID)){
  curGene <- GeneID[a]
  curTable <- matrix(empty_4_1_2[which(empty_4_1_2[,"Gene_ID"]==curGene),],ncol = ncol(empty_4_1_2)); colnames(curTable) <- colnames(empty_4_1_2)
  TPM <- max(as.numeric(curTable[,"tpm"]))
  GeneEmp4 <- rbind(GeneEmp4,c(curGene,TPM))
}

#-----------TV2_2_1_2

TV2_2_1_2 <- as.matrix(read.delim("kal_TV2_2_1_2.tsv",sep="\t",header=T))

ID <- lapply(strsplit(TV2_2_1_2[,1],"\\|"),function(x) x[1:2])

TV2_2_1_2[,1] <- as.vector(unlist(lapply(ID,function(x) x[1])))

TV2_2_1_2 <- cbind(TV2_2_1_2,"Gene_ID"=c(as.vector(unlist(lapply(ID,function(x) x[2])))))

GeneID <- unique(TV2_2_1_2[,"Gene_ID"])

GeneTV2_2 <- matrix(data=NA,ncol=2)

for(a in 1:length(GeneID)){
  curGene <- GeneID[a]
  curTable <- matrix(TV2_2_1_2[which(TV2_2_1_2[,"Gene_ID"]==curGene),],ncol = ncol(TV2_2_1_2)); colnames(curTable) <- colnames(TV2_2_1_2)
  TPM <- max(as.numeric(curTable[,"tpm"]))
  GeneTV2_2 <- rbind(GeneTV2_2,c(curGene,TPM))
}

#-----------TV2_3_1_2

TV2_3_1_2 <- as.matrix(read.delim("kal_TV2_3_1_2.tsv",sep="\t",header=T))


ID <- lapply(strsplit(TV2_3_1_2[,1],"\\|"),function(x) x[1:2])

TV2_3_1_2[,1] <- as.vector(unlist(lapply(ID,function(x) x[1])))
TV2_3_1_2 <- cbind(TV2_3_1_2,"Gene_ID"=c(as.vector(unlist(lapply(ID,function(x) x[2])))))

GeneID <- unique(TV2_3_1_2[,"Gene_ID"])

GeneTV2_3 <- matrix(data=NA,ncol=2)

for(a in 1:length(GeneID)){
  curGene <- GeneID[a]
  curTable <- matrix(TV2_2_1_2[which(TV2_2_1_2[,"Gene_ID"]==curGene),],ncol = ncol(TV2_2_1_2)); colnames(curTable) <- colnames(TV2_2_1_2)
  TPM <- max(as.numeric(curTable[,"tpm"]))
  GeneTV2_3 <- rbind(GeneTV2_3,c(curGene,TPM))
}

#-----------TV2_4_1_2

TV2_4_1_2 <- as.matrix(read.delim("kal_TV2_4_1_2.tsv",sep="\t",header=T))

ID <- lapply(strsplit(TV2_4_1_2[,1],"\\|"),function(x) x[1:2])

TV2_4_1_2[,1] <- as.vector(unlist(lapply(ID,function(x) x[1])))
TV2_4_1_2 <- cbind(TV2_4_1_2,"Gene_ID"=c(as.vector(unlist(lapply(ID,function(x) x[2])))))

GeneID <- unique(TV2_4_1_2[,"Gene_ID"])

GeneTV2_4 <- matrix(data=NA,ncol=2)

for(a in 1:length(GeneID)){
  curGene <- GeneID[a]
  curTable <- matrix(TV2_4_1_2[which(TV2_4_1_2[,"Gene_ID"]==curGene),],ncol = ncol(TV2_4_1_2)); colnames(curTable) <- colnames(TV2_4_1_2)
  TPM <- max(as.numeric(curTable[,"tpm"]))
  GeneTV2_4 <- rbind( GeneTV2_4 ,c(curGene,TPM))
}

####################################################################################################
###############################14/12/2017###############################################
#-Analyse gene expression(edgeR)-####

#load the package

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)


#read in the count file created from htseq-count

Data_Emp2 <- data.frame(empty_2_1_2)
Data_Emp3 <- data.frame(empty_3_1_2)
Data_Emp4 <- data.frame(empty_4_1_2)

Data_TV2_2 <- data.frame(TV2_2_1_2)
Data_TV2_3 <- data.frame(TV2_3_1_2)
Data_TV2_4 <- data.frame(TV2_4_1_2)

# the counts with the row names as the gene ids and the column names as the
#sample ids

cout_Emp2 <- Data_Emp2[, -c(1, ncol(Data_Emp2))]
rownames (cout_Emp2 ) <- Data_Emp2[ , 1 ] # gene names in the raws

cout_Emp3 <- Data_Emp3[, -c(1, ncol(Data_Emp3))]
rownames (cout_Emp3 ) <- Data_Emp3[ , 1 ] # gene names in the raws

cout_Emp4 <- Data_Emp4[, -c(1, ncol(Data_Emp4))]
rownames (cout_Emp4 ) <- Data_Emp4[ , 1 ] # gene names in the raws

cout_TV2_2 <- Data_TV2_2[, -c(1, ncol(Data_TV2_2))]
rownames (cout_TV2_2 ) <- Data_TV2_2[ , 1 ]

cout_TV2_3 <- Data_TV2_3[, -c(1, ncol(Data_TV2_3))]
rownames (cout_TV2_3 ) <- Data_TV2_3[ , 1 ]

cout_TV2_4 <- Data_TV2_4[, -c(1, ncol(Data_TV2_4))]
rownames (cout_TV2_4 ) <- Data_TV2_4[ , 1 ]


#summaries

dim(cout_Emp2) 
dim(cout_Emp3) 
dim(cout_Emp4) 
dim(cout_TV2_2) 
dim(cout_TV2_3) 
dim(cout_TV2_4) 

#some of the columns aren't numerical so in order to transform them into numeric:

install.packages("ISLR")
library(ISLR)

#change the factor variables into numeric variables 

cor(cout_Emp2) #see the type of variables

indx <- sapply(cout_Emp2, is.factor)
cout_Emp2[indx] <- lapply(cout_Emp2[indx], function(x) as.numeric(as.character(x))) #transform the variables in numeric


colSums(cout_Emp2)/1e06 #library sizes in millions of reads

#number of genes with low counts

table(rowSums(cout_Emp2))[1:30]

#building the edgeR object 

####################################################################################################
###############################9/1/2018###############################################

install.packages("ISLR")
library(ISLR)

####EMPTY TABLES-from factor to numeric-###

#change the factor variables into numeric variables 

cor(cout_Emp2) #see the type of variables
cor(cout_Emp3)
cor(cout_Emp4)
cor(cout_TV2_2)
cor(cout_TV2_3)
cor(cout_TV2_4)


indx3 <- sapply(cout_Emp3, is.factor)
cout_Emp3[indx3] <- lapply(cout_Emp3[indx3], function(x) as.numeric(as.character(x))) #transform the variables in numeric

indx4 <- sapply(cout_Emp4, is.factor)
cout_Emp4[indx4] <- lapply(cout_Emp4[indx4], function(x) as.numeric(as.character(x)))

indx_TV2 <- sapply(cout_TV2_2, is.factor)
cout_TV2_2[indx_TV2] <- lapply(cout_TV2_2[indx_TV2], function(x) as.numeric(as.character(x))) #transform the variables in numeric

indx_TV3 <- sapply(cout_TV2_3, is.factor)
cout_TV2_3[indx_TV3] <- lapply(cout_TV2_3[indx_TV3], function(x) as.numeric(as.character(x)))

indx_TV4 <- sapply(cout_TV2_4, is.factor)
cout_TV2_4[indx_TV4] <- lapply(cout_TV2_4[indx_TV3], function(x) as.numeric(as.character(x)))


#PUT THE DATA IN MILLIONS OF READS

cout_Emp2[,2:4]

(colSums(cout_Emp2[,1:4]))/1e06
(colSums(cout_Emp3[,1:4]))/1e06
(colSums(cout_Emp4[,1:4]))/1e06
(colSums(cout_TV2_2[,1:4]))/1e06
(colSums(cout_TV2_3[,1:4]))/1e06
(colSums(cout_TV2_4[,1:4]))/1e06 #library sizes in millions of reads


#number of genes with low counts
table(rowSums(cout_Emp2))[1:30]




###################################################################################
#################################22/01/2018#######################################

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

colnames(readCounts)<- samples

#put the header in the table
CountDF<- as.data.frame(readCounts)

names(CountDF)<-c("empty_2", "empty_3", "empty_4", "TV2_2", "TV2_3", "TV2_4")



lengthDF<- as.data.frame(readlenght)
names(lengthDF)<-c("empty_2", "empty_3", "empty_4", "TV2_2", "TV2_3", "TV2_4")


#-code from:bioconductor-####
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
library("edgeR")

normFactors <- readlenght/exp(rowMeans(log(readlenght)))
o <- log(calcNormFactors(readCounts/norMat)) + log(colSums(readCounts/normMat))
y <- DGEList(readCounts)
y$offset <- t(t(log(normMat)) + o)
readsf<- as.data.frame(y$offset)
names(readsf)<-c( "empty_2", "empty_3", "empty_4", "TV2_2", "TV2_3", "TV2_4")


#-dataframe:CountDF-####

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

#-------------------------------------------------24/01/2018-----------

d<- calcNormFactors(cds)
d$samples
par(mfrow= c(1,2))

#before normalization 
maPlot(cds$counts[,1], cds$counts[,2], normalize = TRUE, pch=19, 
       cex= 0.4, ylim=c(-8,8))
grid(col="blue")
abline(h=log2(cds$samples$norm.factors[2]/d$samples$norm.factors[1]),
       col="red", lwd=4)
eff.libsize<- cds$samples$lib.size * cds$samples$norm.factors

#after normalization
maPlot(cds$counts[,1]/eff.libsize[1], d$counts[,2]/eff.libsize[2],
       normalize = FALSE, pch=19, cex= 0.4, ylim=c(-8,8))
grid(col="blue")


#-----------------------------------------

rawCountTable<- CountDF
sampleInfo <- read.table("SampleInfo.csv", header=TRUE, sep=",", row.names = 1)

nrow(rawCountTable)

dgeFull <- DGEList(rawCountTable, remove.zeros = TRUE)

dgeFull
dgeFull$counts

dgeFull$samples

dgeFull$samples$condition<- relevel(sampleInfo$Condition, ref = "wt")

dgeFull$samples$genotype<- sampleInfo$Genotype
dgeFull

#Data exploration and quality assessment

#exploratory analysis performed on log2-transformed counts to avoid 
#problems duo to skewness of the count distribution

pseudoCounts <- log2(dgeFull$counts+1)
head(pseudoCounts)

#histogram for pseudo-counts

hist(pseudoCounts[, "empty_2"], main="", xlab="counts")

#Bloxplot for pseudo-counts

par(mar=c(8,4,1,2))
boxplot(pseudoCounts, col="gray", las=3, cex.names=1)

#MA-plots between first two WT samples (using  limma package)

limma::plotMA(pseudoCounts[,1:2], xlab="M", ylab="A", main="")
abline(h=0, col="red")

#MDS for pseudo-coounts(limma package)

library(RColorBrewer)

colConditions<- brewer.pal(3, "Set2")
colConditions<- colConditions[match(sampleInfo$Condition,
                                    levels(sampleInfo$Condition))]

pchGenotypes <- c(8,15,16)[match(sampleInfo$Condition,
                                 levels(sampleInfo$Condition))]

plotMDS(pseudoCounts, pch=pchGenotypes, col=colConditions)
legend("center", lwd=2, col=brewer.pal(3, "Set2")[1:2],
       legend = levels(sampleInfo$Condition))
legend("bottomright", pch = c(8,15,16),
       legend=levels(sampleInfo$Genotype))

#Normalization; Compute normalization factors
dgeFull<- calcNormFactors(dgeFull, method = "TMM")
dgeFull

head(dgeFull$counts)
head(dgeFull$samples)

normCounts <- cpm(dgeFull)
pesudoNormCounts <- cpm(dgeFull, log=TRUE, prior.count=1)
par(mar=c(8,4,1,2))
boxplot(pesudoNormCounts, col = "gray", las=3, cex.names=1)

plotMDS(pesudoNormCounts, pch= pchGenotypes, col=colConditions)
legend("top", lwd=2, col=brewer.pal(3, "Set2")[1:2],
       legend= levels(sampleInfo$Condition))
legend("bottomleft", pch= c(8,15,16),
       legend= levels(sampleInfo$Condition))


#new argument object

dgeFull.group<- DGEList(rawCountTable, remove.zeros = TRUE,
                        group = dgeFull$samples$condition) #REMOVED 19598 ROWS WITH ALL ZERO COUNTS

dgeFull.group$samples$norm.factors<- dgeFull$samples$norm.factors
dgeFull.group

#estimate dispersion

dgeFull.group<- estimateCommonDisp(dgeFull.group)
dgeFull.group<- estimateTagwiseDisp(dgeFull.group)
dgeFull.group

#function estimateGLMRoboutDisp used if the previous approach seems to fail
plotBCV(dgeFull.group)

#perform the test: differences between the conditions "WT" AND "OE"

dgeExactTest<- exactTest(dgeFull.group)
dgeExactTest

#topTags corrects the p-values
resExactTest<- topTags(dgeExactTest, n= nrow(dgeExactTest$table))
head(resExactTest$table)

#p-value and (BH) adjusted p-value distribution 

par(mfrow= c(1,2))
hist(resExactTest$table$PValue, xlab = "p-value", main="raw p-values")
hist(resExactTest$table$FDR, xlab = "p-value", main="adjusted p-values")

#genes with a FDR smaller than 5% and a log Fold Change larger than 1 or smallet
#than -1 are extracted

selectedET<- resExactTest$table$FDR < 0.05 & abs(resExactTest$table$logFC)>1
selectedET<- resExactTest$table[selectedET,]
nrow(selectedET) #1536

head(selectedET)

selectedET$updown<- factor(ifelse(selectedET$logFC >0 , "up", "down"))
head(selectedET)


#second approach: GLM with condition and genotype effects

#estimate dispersion

design.matrix<- model.matrix(~ dgeFull$samples$condition +
                               dgeFull$samples$genotype)
design.matrix

dgeFull<- estimateCommonDisp(dgeFull, design.matrix)
dgeFull

#assessment of the quality of the variability estimation
plotBCV(dgeFull)

#fit GLM and perform the test

fit<- glmFit(dgeFull, design.matrix)
fit

#LOG-RATIO TEST
dgeLRtest<- glmLRT(fit, coef = 2)
dgeLRtest

#DGEs extracting using topTags
resLRT <- topTags(dgeLRtest, n= nrow(dgeFull$counts))
head(resLRT$table)

selectedLRT<- resLRT$table$FDR < 0.05 & abs(resLRT$table$logFC)
selectedLRT <- resLRT$table[selectedLRT, ]
nrow(selectedLRT)

head(selectedLRT)

#comparing the two approaches bellow
vd <- venn.diagram(x=list("Exact test"= rownames(selectedET),
                          "GLM"= rownames(selectedLRT)),
                   fill= brewer.pal(3, "Set2")[1:2], filename=NULL)

#exploratory analysis of DGEs

install.packages(c("RColorBrewer","mixOmics"))
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")

biocLite("HTSFilter")
library(HTSFilter)
library(Biobase)

browseVignettes("HTSFilter")

dgeFilt<- HTSFilter(dgeFull)$filteredData
dgeFilt

fit<- glmFit(dgeFilt, design.matrix)
dgeLRTfilt<- glmLRT(fit, coef = 2)
resLRTfilt<- topTags(dgeLRTfilt, n=nrow(dgeFilt$counts))
selectedFilt<- resLRTfilt$table$FDR < 0.05 & abs(resLRTfilt$table$logFC)>1
selectedFilt<- resLRTfilt$table[selectedFilt,]
nrow(selectedFilt)
head(selectedFilt)

#create a MA plot with differentially expressed genes

dgeFilt$samples$group<- dgeFilt$samples$condition
plotSmear(dgeFilt, de.tags = rownames(selectedFilt))

#Volcano Plot

volcanoData<- cbind(resLRTfilt$table$logFC, -log10(resLRTfilt$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
DEGs<- resLRTfilt$table$FDR < 0.05 & abs(resLRTfilt$table$logFC) > 1
point.col<-  ifelse(DEGs, "red", "black")
plot(volcanoData, pch=16, col=point.col, cex=0.5)


#-DESeq2 -####


source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("apeglm")
biocLite("ashr")
library("DESeq2")
library("readr")

sampleInfo <- read.table("SampleInfo.csv", header=TRUE, sep=",", row.names = 1)

head(CountDF)

sampleTable<- data.frame(condition= factor(rep(c("empty", "TV2"), each=3)))
colData<- data.frame(condition=factor(rep(c("kal_empty", "kal_TV2"), each=3)))
rownames(sampleTable)<- colnames(txi_gene$counts)
dds<- DESeqDataSetFromTximport(txi_gene, colData = colData, design= ~condition)

dds<-DESeq(dds)
resultsNames(dds) #lists the coefficients
res<- results(dds, name = "condition_kal_TV2_vs_kal_empty")
resLFC<- lfcShrink(dds, coef = "condition_kal_TV2_vs_kal_empty")

#differential expression analysis

plotMA(res, ylim=c(-2,2))

#shrunken log2 fold changes, removes the noise associated with log2 fold changes

plotMA(resLFC, ylim=c(-2,2))



#idx<- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]

resApe<- lfcShrink(dds, coef=2, type = "apeglm")
resAsh<- lfcShrink(dds, coef=1 , type = "ashr")

#counts reads for a single gene across the groups
plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")

d<- plotCounts(dds, gene = which.min(res$padj), intgroup = "condition",
               returnData = TRUE)

library("ggplot2")
ggplot(d, aes(x=condition, y=count))+ 
  geom_point(position = position_jitter(w=0.1, h=0))+
  scale_y_log10(breaks=c(25,100,400))

source("https://bioconductor.org/biocLite.R")
biocLite("Deseq2") 


mcols(res)$description


colData(dds)

#create a copy of the DESeqDataSet, so that we can rerun the analysis
#using a multi-factor design

ddsMF<- dds

levels(ddsMF$type)

levels(ddsMF$condition)<- sub("-.", "", levels(ddsMF$condition))
levels(ddsMF$condition)


ddsMF<- DESeq(ddsMF)
resMF<- results(ddsMF)

#extracting transformed values

vsd<- vst(dds, blind = FALSE)
rld<- rlog(dds, blind = FALSE)
head(assay(vsd), 3)

#effects of transformations on the variance
ntd<- normTransform(dds)

source("https://bioconductor.org/biocLite.R")
biocLite("vsn")
library("vsn")

meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

#PCA plot of samples

plotPCA(vsd, intgroup= "condition") 

install.packages("ggplot2")

pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar<- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed()

#dispersion plot and fotting alternatives

source("http://bioconductor.org/biocLite.R")
biocVersion()

plotDispEsts(dds)

#independent filtering of results

metadata(res)$alpha

metadata(res)$filterThreshold

plot(metadata(res)$filterNumRej,
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)


resNoFilt<- results(dds, independentFiltering = FALSE)
addmargins(table(filtering=(res$padj < .1),
                 noFiltering=(resNoFilt$padj < .1)))

par(mfrow=c(2,2), mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA<- results(dds, lfcThreshold = .5, altHypothesis = "greaterAbs")
resLA<- results(dds, lfcThreshold = .5, altHypothesis = "lessAbs")
resG<- results(dds, lfcThreshold = .5, altHypothesis = "greater")
resL<- results(dds, lfcThreshold = .5, altHypothesis = "less")
drawLines <- function() abline(h=c(-0.5,0.5), col="dodgerblue", lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()

#Access to calculated values
mcols(dds, use.names = TRUE)[1:4, 1:4] #VOLTAR A REPETIR

substr(names(mcols(dds)),1,10)

head(assay(dds)["cooks"])

head(dispersions(dds))

sizeFactors(dds)

head(coef(dds))

attr(dds, "betaPriorVar")

priorInfo(resLFC)
priorInfo(resApe)
priorInfo(resAsh)
dispersionFunction(dds)
attr(dispersionFunction(dds), "dispPriorVar")
metadata(dds)["version"]

#Count outlier detection

w<- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx<- !is.na(w)
plot(rank(w[idx]), maxCooks[idx], xlab="rank of wald statistic",
     ylab= "maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m<- ncol(dds)
p<- 3
abline(h=qf(.99, p, m-p))

#Independent filtering and multiple testing

plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))


#DESEQ VS EDGER-####

de<- estimateSizeFactors(dds)

de<- estimateDispersions(d)

#DESEQ:file:///C:/Users/Ana%20Sofia/Downloads/DESeq2%20(1).pdf-####

source("https://bioconductor.org/biocLite.R")
biocLite("pasilla")
biocLite("Biobase")
library("pasilla")
library("Biobase")

head(CountDF)
head(sampleInfo)

library("tximport")
library("readr")
library("tximporData")

