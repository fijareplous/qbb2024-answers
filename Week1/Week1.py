#!/usr/bin/env python

import numpy
import scipy

#Question 1.2 (3x)

#Total length of genome:
genome_size = 1000000
#Read length:
read_length = 100 
#Coverage:
coverage = 3
# 1Mbp * 3x coverage / 100bp reads
num_reads = int((genome_size * coverage)/read_length)


#Use an array to keep track of the coverage at each position in the genome
genome_coverage = numpy.zeros(genome_size, int)


for i in range(num_reads):
  startpos = numpy.random.randint(0, genome_size - read_length + 1)
  endpos = startpos + read_length
  genome_coverage[startpos:endpos] += 1

#Save coverages as .txt for R plotting
numpy.savetxt("genome_coverage_3x.txt", genome_coverage)

#Get the range of coverages observed
maxcoverage = max(genome_coverage)
xs = list(range(0, maxcoverage+1))

#Get the poisson pmf at each of these
poisson_estimates = scipy.stats.poisson.pmf(xs, coverage)

#Get normal pdf at each of these (i.e. the density between each adjacent pair of points)
normal_estimates = scipy.stats.norm.pdf(xs, numpy.mean(genome_coverage), numpy.std(genome_coverage))


#Question 1.5 (10x)

#Total length of genome:
genome_size = 1000000
#Read length:
read_length = 100 
#Coverage:
coverage = 10
# 1Mbp * 10x coverage / 100bp reads
num_reads = int((genome_size * coverage)/read_length)


#Use an array to keep track of the coverage at each position in the genome
genome_coverage = numpy.zeros(genome_size, int)


for i in range(num_reads):
  startpos = numpy.random.randint(0, genome_size - read_length + 1)
  endpos = startpos + read_length
  genome_coverage[startpos:endpos] += 1

#Save coverages as .txt for R plotting
numpy.savetxt("genome_coverage_10x.txt", genome_coverage)

#Get the range of coverages observed
maxcoverage = max(genome_coverage)
xs = list(range(0, maxcoverage+1))

#Get the poisson pmf at each of these
poisson_estimates = scipy.stats.poisson.pmf(xs, coverage)

#Get normal pdf at each of these (i.e. the density between each adjacent pair of points)
normal_estimates = scipy.stats.norm.pdf(xs, numpy.mean(genome_coverage), numpy.std(genome_coverage))


#Question 1.6

#Total length of genome:
genome_size = 1000000
#Read length:
read_length = 100 
#Coverage:
coverage = 30
# 1Mbp * 30x coverage / 100bp reads
num_reads = int((genome_size * coverage)/read_length)


#Use an array to keep track of the coverage at each position in the genome
genome_coverage = numpy.zeros(genome_size, int)


for i in range(num_reads):
  startpos = numpy.random.randint(0, genome_size - read_length + 1)
  endpos = startpos + read_length
  genome_coverage[startpos:endpos] += 1

#Save coverages as .txt for R plotting
numpy.savetxt("genome_coverage_30x.txt", genome_coverage)

#Get the range of coverages observed
maxcoverage = max(genome_coverage)
xs = list(range(0, maxcoverage+1))

#Get the poisson pmf at each of these
poisson_estimates = scipy.stats.poisson.pmf(xs, coverage)

#Get normal pdf at each of these (i.e. the density between each adjacent pair of points)
normal_estimates = scipy.stats.norm.pdf(xs, numpy.mean(genome_coverage), numpy.std(genome_coverage))


#Question 2.1 

#Saved the following as edges.txt, as I started 2.1 after txt had been provided:

'''
digraph{
TTC -> TCA
TCT -> CTT
TGA -> GAT
CAT -> ATT
ATT -> TTC
ATT -> TTG
CTT -> TTA
GAT -> ATT
TAT -> ATT
TCA -> CAT
TTA -> TAT
TTG -> TGA
TTC -> TCT
ATT -> TTT
}
'''


