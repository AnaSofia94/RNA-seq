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
time ls $Clean_Data/*.fq | $parallel '$srun shifter --image=argrosso/htspreprocessing:0.1.1 

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitc$

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
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitc$

