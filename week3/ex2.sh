#!/usr/bin/env bash

### Question 2.1 ###
wget https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz
bwa index sacCer3.fa

#From working with yeast for two years, I know that budding yeast have 16 chromosomes.

### Questions 2.2-2.4 Code ###

#for each sample (numbers given below)
for sample in 09 11 23 24 27 31 35 39 62 63
    do 
        #Aligns reads in FASTQ to reference    
        bwa mem -t 4 -r "@RG\tID:$sample\tSM:$sample" sacCer3.fa A01_$sample.fastq  > A01_${sample}.sam
        #Sort .sam and convert to .bam
        samtools sort -@ 4 -O bam -o A01_${sample}.bam A01_${sample}.sam
        #Indexes .bam for easier processing
        samtools index A01_${sample}.bam
    done


### Questions 2.2-2.6 Answers ### 

#Question 2.2
#@ for reads
grep -v '^@' A01_09.sam | wc -l
#There are 669548 reads

#Question 2.3
#Finds "chrIII" to identify reads to chromosome III
grep -w "chrIII" A01_09.sam | wc -l
#There are 18195 alignments for chromosome III

#Question 2.4
#Loaded A01_09.bam into the IGV app, compared to sacCer3 reference genome 
#4x coverage seems reasonable; there are areas that are covered more (perhaps 7 or 8x) and others where there is no coverage. 

#Question 2.5
#There are three obvious SNPS: one near 113,132, another close to 113,207, and a third at roughly 113,326. There don't appear to
#be any others, as only these three positions have polymorphisms seen across all reads at that position.

#Question 2.6 
#The SNP is at approximately position 825,834. It falls just downstream of the SCC2 gene but is not located within a coding region.
