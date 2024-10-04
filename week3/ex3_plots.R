library(tidyverse)
library(ggplot2)

allelefreq <- read.delim("~cmdb/qbb2024-answers/week3/AF.txt")
depth <- read.delim("~cmdb/qbb2024-answers/week3/DP.txt")

#Allele frequency graph
ggplot(data=allelefreq, mapping = aes(x = allele_frequency)) + geom_histogram(bins = 11, 
  fill = "purple", color = "black") +
  aes(y = after_stat(count)/sum(after_stat(count))) + 
  labs(title = "Spectrum of Allele Frequencies Among Variants",y = "Frequency of Variants", 
  x = "Allele (SNP) Frequency") 

ggsave("~cmdb/qbb2024-answers/week3/AF_plot.png")

#Depth of reads graph
ggplot(data=depth, mapping = aes(x = read_depth)) + geom_histogram(bins=21, 
  fill = "purple", color = "black") +
  xlim(0,20) + aes(y = after_stat(count)/sum(after_stat(count))) + 
  labs(title = "Spectrum of Read Depths Across Variants",y = "Frequency of Variants", 
  x = "Read Depths") 

ggsave("~cmdb/qbb2024-answers/week3/DP_plot.png")


#Question 3.1
#This graph looks as expected. The majority of alleles have a frequency of about 50%, although
#many are more or less common. This makes sense given that many SNPs are not located in "important"
#regions that could have serious consequences and thus will have a frequency of around 50%. Also, I 
#I would say this is a normal/Gaussian distribution. 

#Question 3.2
#This graph also looks as expected. It was previously estimated that there was roughly 4x coverage
#of the genome, and looking at the graph, the average appears to be approximately 4x. Also, I would 
#say this graph represents a Poisson distribution. 
