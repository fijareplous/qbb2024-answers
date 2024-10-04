#!/usr/bin/env bash

### Question 1.1 ###
FASTQC A01_09.fastq

#The sequence length is 76 


### Question 1.2 ###
wc -l A01_09.fastq
#2678192 from previous line, divide by 4:
expr 2678192 / 4
#There are 669548 reads in the file


### Question 1.3
#Used code from question 2.1 to download .fa file for S. cerevisiae reference genome

# "grep -v '^>'" filters out >; "tr -d '\n'" removes new line characters; "wc -c" counts the number of characters in the string
grep -v '^>' sacCer3.fa | tr -d '\n' | wc -c
#Gives an answer of 12157105

#To get coverage, divide answer by number of reads:
expr 12157105 / 2678192
#The average depth coverage is about 4x


### Question 1.4 ###
#-h makes it human readable:
du -h A01_*.fastq | sort -h
#A01_62.fastq is the largest (149M), and A01_27.fastq is the smallest (110M)

## Question 1.5 ###
#Had already run FASTQC in 1.1
#The median appears to be about 35.
#Given that the probability of error can be found by P(error)=10^−(Q/10) where Q is the median quality score, the probability 
#an error at a given base is 3.16x10^-4 or .0316%.
#The ends of the reads have slightly lower quality scores (which makes sense given how sequencing works) compared to the middle,
#but overall there is not much difference. 






