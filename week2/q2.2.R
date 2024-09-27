library(tidyverse)
library(ggthemes)

data <- read.delim("~/qbb2024-answers/week2/snp_counts.txt")


ggplot(data = data, mapping = aes(x = MAF, y = log2(Enrichment), color = Feature)) + 
  geom_line() + xlab("Minor Allele Frequency (MAF)") + ylab("log2 SNP Enrichment)") + 
  scale_color_discrete(labels = c("cCREs", "Exons", "Introns", "Other"))