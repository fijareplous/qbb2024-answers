library(tidyverse)
library(broom)
library(DESeq2)

#Part 1

#Set WD
setwd("~/qbb2024-answers/week7/")

#Load data
counts_df <- read_delim("gtex_whole_blood_counts_downsample.txt")

#Load metadata
metadata_df <- read_delim("gtex_metadata_downsample.txt")

#Make gene name the row names in the data
counts_df <- column_to_rownames(counts_df, var = "GENE_NAME")

#Make the subject ID the row names in the metadata
metadata_df <- column_to_rownames(metadata_df, var = "SUBJECT_ID")

#Sanity check to make sure that the row names of the metadata match col names in data
colnames(counts_df) == rownames(metadata_df)
table(colnames(counts_df) == rownames(metadata_df))
#They all came back as "true" which is lovely

#Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                              colData = metadata_df,
                              design = ~ DTHHRDY + SEX + AGE)

#Apply a variance stabilizing transformation to give high expressed genes more weight
vsd <- vst(dds)

#PCA plots
png(filename="AGE_PCA_Plot.png")
plotPCA(vsd, intgroup = "AGE")
dev.off()
png(filename="SEX_PCA_Plot.png")
plotPCA(vsd, intgroup = "SEX")
dev.off()
png(filename="DTHHRDY_PCA_Plot.png")
plotPCA(vsd, intgroup = "DTHHRDY")
dev.off()

#The first PC explains 48% of the variance, while the second explains 7% of the variance
#DTHHRDY (cause of death) appears to be associated with the larger of the two PCAs (48%),
#while AGE seems to be associated with the second largest PCA (7%), as there appears to be
#a weak amount of clustering between different age groups.


#Part 2

#Extract the VST expression matrix and bind it to the metadata
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble()
vsd_df <- bind_cols(metadata_df, vsd_df)

#Test for differential expression of the gene WASH7P
wash7p <- lm(formula = WASH7P ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()

#Considering the p value for sex-differentiated expression is only 2.792437e-01, no—
#there is not significant evidence.

#Test for differential expression of the gene SLC25A47
slc <- lm(formula = SLC25A47 ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()

#Yes, there is sufficient evidence to show sex-differentiated expression, as the p value
#is 2.569926e-02, which is less than the commonly used threshold of 0.05. 

#Fit the regression model using DESeq2
dds <- DESeq(dds)

# Differential expression results for sex
res_sex <- results(dds, name = "SEX_male_vs_female")  %>%
  as_tibble(rownames = "GENE_NAME")

#Remove NA values from results
res_sex <- res_sex %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

#Filter results to show genes that are differentially expressed at an FDR of 10%
res_sex <- res_sex %>%
  filter(padj < 0.1) %>%
  arrange(padj)

dim(res_sex)[1]
#There are 262 genes that fall into this category. 

#Open file with gene mapping information
locations = read.delim("gene_locations.txt")

#Merge gene locations with data; 
genes_by_loc = left_join(locations,res_sex,by="GENE_NAME") %>% arrange(padj)

#The X and Y chromosomes understandably encode the genes that are most strongly 
#upregulated in males versus females, respectively. This makes sense given that 
#the X and Y chromosomes are sex chromosomes and thus contribute significantly to 
#sex-based phenotypes.There are more male upregulated genes, as female was set
#as the reference.


WASH7_DESeq = genes_by_loc %>% filter(GENE_NAME == "WASH7")
SLC_DESeq = genes_by_loc %>% filter(GENE_NAME == "SLC25A47")

#****

#Differential expression results for cause of death
res_death <- results(dds, name = "DTHHRDY_ventilator_case_vs_fast_death_of_natural_causes")  %>%
  as_tibble(rownames = "GENE_NAME")

#Remove NA values from results
res_death <- res_death %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

#Filter results to show genes that are differentially expressed at an FDR of 10%
res_death <- res_death %>%
  filter(padj < 0.1) %>%
  arrange(padj)

dim(res_death)[1]
#There are 16069 genes that fall into this category. 

#Given that cause of death was previously determined to contribute to 48% of the 
#variance, it is understandable that there would be such a high number of genes 
#differentially expressed based on type of death and a comparitively lower number 
#of genes differentially expressed based on sex.


#Part 3

#Generate volcano plot, significant at a 10% FDR & abs(log2FoldChange) greater than 1 
ggplot(data = res_sex, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = (abs(log2FoldChange) > 1 & -log10(padj)>1))) +
  geom_text(data = res_sex %>% filter(abs(log2FoldChange) > 2 & -log10(padj) > 10),
  aes(x = log2FoldChange, y = -log10(padj) + 5, label = GENE_NAME), size = 3,) +
  theme_bw() + theme(legend.position = "none") +
  scale_color_manual(values = c("darkgray", "coral")) +
  labs(y = expression(-log[10]("p-adj")), x = expression(log[2]("fold change")))

ggsave("volcano_plot_sex.png")

