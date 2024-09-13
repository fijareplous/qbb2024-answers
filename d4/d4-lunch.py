#!/usr/bin/env python

import sys

import numpy

#Questions 1 - 4

# get gene-tissue file name
filename = sys.argv[1]
# open file as "fs"
fs = open(filename, mode='r')
# create a dictionary to hold samples for gene-tissue pairs
relevant_samples = {}
# step through file (named fs)
for line in fs:
    # Split each line into "fields"
    fields = line.rstrip("\n").split("\t")
    # Create key from gene and tissue
    key = (fields[0], fields[2])
    # Initialize dict from key with list to hold samples
    relevant_samples[key] = []
#Close the file (good habit)
fs.close()


# Get metadata file name
filename = sys.argv[2]
# open file as "fs"
fs = open(filename, mode='r')
# Skip first line
fs.readline()
# Create dictionary to hold samples for tissue name
tissue_samples = {}
# step through file
for line in fs:
    # Split line into fields
    fields = line.rstrip("\n").split("\t")
    # Create key from gene and tissue
    key = fields[6]
    value = fields[0]
    # Initialize dict from key with list to hold samples and 
    # set default so we don’t need to check if a tissue already added to the dictionary
    tissue_samples.setdefault(key, [])
    tissue_samples[key].append(value)
fs.close()

# get metadata file name
filename = sys.argv[3]
# open file
fs = open(filename, mode='r')
# Skip first 2 lines
fs.readline()
fs.readline()
#Want to show column names and split by tabs
header = fs.readline().rstrip("\n").split("\t")
#"Filter" to everything after column 2
header = header[2:]

#Questions 4 and 5
#Dictionary to pull out the columns for all tissue types
tissue_columns={}
for tissue, samples in tissue_samples.items(): # .items gives keys,values (tissues,sampleID)
    #Same explanation as above for .setdefault()
    tissue_columns.setdefault(tissue, [])
    for sample in samples:
        if sample in header:
            position=header.index(sample) # index tells you where in the 'key' do you see 'value'
            tissue_columns[tissue].append(position) # append allows you to add to list
# print(header)

#To find which tissue types have the most samples:
#Setting base max value of 0 so we have something to append
maxValue=0
#Iterate through tissue_samples
for tissue, samples in tissue_samples.items():
    #If the number of samples is greater than the current max value...
    if len(samples) > maxValue: 
        #...update the max value
        maxValue=len(sample)
        #To be able to ID which type it is
        maxValue(key) = tissue 

#Muscle skeletal has the most number of samples

#Need an insanely high minimum value that we can work backwards from to ID tissue type(s) with least # samples
minValue=20000000
#This code is the same as the code above but reversed to find the minimum
for tissue, samples in tissue_samples.items():
    if len(samples) < minValue:
        minValue=len(sample)
        minValue(key) = tissue #if the sample is smaller than the min value, it becomes the new min value
#print(minValueKey) 
#The leukemia cell line has the least number of samples


#See R file for Question 7 (d4-lunch-Q7)