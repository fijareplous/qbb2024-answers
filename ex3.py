#!/usr/bin/env python3

import sys

#Creating allele frequency list to be imported into file later:
allele_frequencies = []

#Create output files for later and open VCF filel:
vcf = open('/Users/cmdb/qbb2024-answers/week3/biallelic.vcf')
af_file = open('/Users/cmdb/qbb2024-answers/week3/AF.txt', 'w')
dp_file = open('/Users/cmdb/qbb2024-answers/week3/DP.txt', 'w')

#Realized in the next part of the assignment that we definitely need headers...
af_file.write('allele_frequency\n')
dp_file.write('read_depth\n')


for line in vcf:
    # Check if the line is a header to skip metadata lines
    if line.startswith('##'):
        continue 
    #Column headers start with #, want to skip
    if line.startswith('#'):
        continue 
    fields = line.rstrip('\n').split('\t')

    #Question 3.2 (AF code)
    info = fields[7]
    af_value = None
    for entry in info.split(';'):
        #Looking at AF specifically 
        if entry.startswith('AF='):
            af_value = entry.split('=')[1]
            #Want to stop when AF encountered
            break 
    if af_value is not None:
        #Add to AF.txt
        af_file.write(af_value + '\n')

    #Question 3.3 (DP code)
    dp_value = None
    #Want everything from 10th onwards
    dp = fields[9:]
    for entry in dp:
        dp_data = entry.split(':')
        #Get the DP value
        dp_value = dp_data[2]
        #Add to DP.txt
        if dp_value is not None: 
                dp_file.write(dp_value + '\n')

vcf.close()
af_file.close()
dp_file.close()


