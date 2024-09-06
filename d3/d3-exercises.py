#!/usr/bin/env python

#File is in the d3 folder under qbb2024-answers 
#(GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct)

import sys

import numpy

#open file 
fs = open(sys.argv[1], mode = "r")

#skip first two lines
fs.readline()
fs.readline()
line = fs.readline()

#Split column heads by tabs and skip first two entries
fields = line.strip("\n").split("\t")
tissues = fields[2:]
#print(tissues)

#Create way to hold genes names, IDs, and expression
gene_names = []
gene_IDs = []
expression = []

#for each line: 1) split line, 2) save field 0 into gene names, 3) save 2+ into expression values
for line in fs:
    fields = line.strip("\n").split("\t")
    gene_IDs.append(fields[0])
    gene_names.append(fields[1])
    expression.append(fields[2:])

fs.close()

#Part 2 (Group Work):
tissues = numpy.array(tissues)
gene_IDs = numpy.array(gene_IDs)
gene_names = numpy.array(gene_names)
expression = numpy.array(expression, dtype = float)

#Without specifying that the values are floaters, numpy will automatically place the values in the array as strings, which
#are difficult to work with numerically. As such we need to specify they are floaters (not integers, which do not contain 
#decimals)

'''
print(tissues)
print(gene_IDs)
print(gene_names)
print(expression)
'''

#Part 3:
#Instructed to skip

#Part 4:
mean = numpy.mean(expression[:10], axis = 1)
# print("The mean expression value for the first ten genes across tissue types is:")
# print(mean)

#Question from Part 4 is not applicable, as Part 3 was skipped.

#Part 5:
total_mean = numpy.mean(expression, axis = 1)
#print("The mean expression value for all genes in the dataset across tissue types is:")
#print(total_mean)

total_median = numpy.median(expression, axis = 1)
#print("The median expression value for all genes in the dataset across tissue types is:")
#print(total_median)

#The median expression values appear to be smaller than the mean expression values. This means that the data is positively 
#skewed (or skewed to the right).

#Part 6:
log = numpy.log2(expression+1)
median_log = numpy.median(log, axis = 1)
print("Median of log values:")
print(median_log)

mean_log = numpy.mean(log, axis = 1)
print("Mean of log values:")
print(mean_log)

#The mean and median values of the log are skewed in the same direction as described in the previous question and with similar 
#degrees of difference between the corresponding mean/median. 

# print(log)

#Part 7:
log_copy = numpy.copy(log)
sorted = numpy.sort(log_copy, axis = 1)
diff_array = sorted[:,-1] - sorted[:,-2]
#print(numpy.sort(diff_array))

#Part 8:
print(len(numpy.sort(diff_array[diff_array > 10])))
#There are 33 genes where the difference between the highest and second highest tissue expression values
# is greater than 10 