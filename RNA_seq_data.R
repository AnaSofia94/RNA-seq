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


  