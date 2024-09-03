#Q1
library(tidyverse)
gtex <- read_delim("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

#Q2
glimpse(gtex)

#Q3
gtexNew <- filter(gtex,SMGEBTCHT == "TruSeq.v1")

#Q4
ggplot(data=gtexNew, mapping = aes(x = SMTSD)) + geom_bar() + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) +
  xlab("Tissue Type") + ylab("Number of Subjects") + 
  ggtitle("Number of Samples of Each Tissue Type")

#Q5
ggplot(data=gtexNew, mapping=aes(x=SMRIN)) + geom_bar() + xlab("RIN Number") +
  ylab("Frequency") + ggtitle("Distribution of RIN Numbers Across Samples")
#The data is unimodal, skewed right. The average appears to be around 8.

#Q6
ggplot(data=gtexNew, mapping=aes(x=SMTSD, y= SMRIN)) + geom_violin() + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) +
  xlab("Tissue Type") + ylab("RIN") + 
  ggtitle("Distribution of RIN Numbers Across Tissue Types")
#I chose a violin plot to be able to show the relative abundances of RINs 
#across various tissue types
#Fibroblasts and cancerous cells showed greater average RINs than other tissues
#and with less variation. This is due to the fact that samples were taken directly 
#from a cell culture rather than from deceased individuals.

#Q7
#SMGNSDTC = number of genes detected
ggplot(data=gtexNew, mapping=aes(x=SMTSD, y=SMGNSDTC)) + geom_violin() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) +
  xlab("Tissue Type") + ylab("Number of Genes") + 
  ggtitle("Number of Genes Detected Per Sample Across Tissue Types")
#I again chose a violin plot to visualize the number of genes across various tissue
#types and their relative abundances among samples. Testis averaged the hightest 
#in terms of number of genes, which is supported by the fact that testicular cells
#use transcriptional scanning to modulate the rate of mutation within the germline
#(Xia et al., 2021) 

#Q8
ggplot(data=gtexNew, mapping=aes(x=SMTSISCH, y=SMRIN)) + 
  geom_point(size = 0.5, alpha = .5) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) +
  xlab("Ischemic Time") + ylab("RIN") +
  ggtitle("Relationship Between Ischemic Time and RIN") +
  facet_wrap(. ~ SMTSD) + geom_smooth(method = "lm")
#The data tends to show that as the ischemic time increases, the RIN number decreases.
#In other words, as the amount of time after the blood supply has been reduced
#increases, the quality of the RNA decreases, which makes intuitive sense. Some
#tissues see a degrade in RNA faster than others (e.g. fibroblasts retain quality
#RNA longer than the esophagus)

#Q9
ggplot(data=gtexNew, mapping = aes(x=SMTSISCH, y=SMRIN)) + 
  geom_point(mapping=aes(color = SMATSSCR), size = 0.5, alpha = .5) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) +
  xlab("Ischemic Time") + ylab("RIN") +
  ggtitle("Relationship Between Ischemic Time, RIN, and Autolysis Scores Across Tissue Types") + 
  geom_smooth(method = "lm") +
  facet_wrap(. ~ SMTSD) 
#I don't see much of a relationship, although the data is very dense and small 
#on the screen. Colon cells, however, seem to have higher autolysis scores compared
#to other cell types.