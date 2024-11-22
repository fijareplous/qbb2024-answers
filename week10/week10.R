library(dplyr)
library(ggplot2)
library(tidyverse)

#Set working directory
setwd("~/qbb2024-answers/week10")

#Import file
nuclei_signal <- read.table("nuclei_signal.txt")

#Add column names
colnames(nuclei_signal) <- c("Gene","nascentRNA","PCNA_signal", "log2ratio")

#Format as data frame
nuclei_signal <- as.data.frame(nuclei_signal)

#Plot nascent RNA signal for each gene
ggplot(data = nuclei_signal, mapping = aes(x = Gene, y = nascentRNA)) + 
  geom_violin() + ggtitle("Nascent RNA Signal for Gene Knockdowns in HeLa Cells") +
  ylab("Nascent RNA fluorescence signal") + xlab("siRNA Knockdown")
ggsave("nRNA_violin.png")
 
#Plot PCNA signal for each gene
ggplot(data = nuclei_signal, mapping = aes(x = Gene, y = PCNA_signal)) + 
  geom_violin() + ggtitle("PCNA Signal for Gene Knockdowns in HeLa Cells") +
  ylab("PCNA Signal") + xlab("siRNA Knockdown")
ggsave("PCNA_violin.png")

#Plot log2 ratio for each gene
ggplot(data = nuclei_signal, mapping = aes(x = Gene, y = log2ratio)) + 
  geom_violin() + ggtitle("Log2 Ratio of nascent RNA Signal to PCNA Signal for Gene
  Knockdowns in HeLa Cells") + ylab("Log2(nRNA signal to PCNA signal") + 
  xlab("siRNA Knockdown")
ggsave("log2Ratio_violin.png")