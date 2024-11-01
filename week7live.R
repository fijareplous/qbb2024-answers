library(tidyverse)
library(broom)
library(DESeq2)

#Pasilla master splicing regulator, mutated it to see what genes impacted in 
#Drosophila 

# set your working directory to where your data and output will be stored
setwd("~/qbb2024-answers/week7")

# load the gene expression counts
counts_df <- read_delim("pasilla_gene_counts.tsv")

# move the gene_id column to rownames, so that the contents of the
# tibble is entirely numeric
counts_df <- column_to_rownames(counts_df, var = "gene_id")

# look at first five rows
counts_df[1:5,]

# load the metadata
metadata_df <- read_delim("pasilla_metadata.csv")

#Want the metadata to match in terms of data column names, so move the sample 
#IDs from the first column to rownames
metadata_df <- column_to_rownames(metadata_df, var = "SAMPLE_ID")


# check that the columns of the counts are identical and in the same order as the
# rows of the metadata
colnames(counts_df) == rownames(metadata_df)
table(colnames(counts_df) == rownames(metadata_df))

# create a DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                              colData = metadata_df,
                              design = ~ condition + type)
#Design specifies the regression formula 

# apply VST normalization
vsd <- vst(dds)

# apply and plot principal components
plotPCA(vsd, intgroup = "condition")
#PC1 contributes to more than half of the variance (treated vs untreated)
plotPCA(vsd, intgroup = "type")
#PC2 is the type of sequencing, the fact that it's explaining 30% of variance is concerning
plotPCA(vsd, intgroup = c("condition", "type"))

# if I wanted "untreated" to be the reference level
dds$condition <- relevel(dds$condition, ref = "untreated")

# use DESeq2 to perform differential expression analysis across all genes
dds <- DESeq(dds)

pasilla_res <- results(dds, name = "condition_treated_vs_untreated") %>%
  as_tibble(rownames = "GENE_ID")

#how many genes (rows) differentially expressed at a false discovery rate of <10%
pasilla_res%>%
  filter(padj < 0.1) %>%
  nrow()
#1330
#False positive set really low will mean that you will still miss stuff so don't want
#super low

# create a volcano plot

ggplot(data = pasilla_res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = (abs(log2FoldChange) > 2 & pvalue < 1e-20))) +
  geom_text(data = pasilla_res %>% filter(abs(log2FoldChange) > 2 & pvalue < 1e-50),
            aes(x = log2FoldChange, y = -log10(pvalue) + 5, label = GENE_ID), size = 3,) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("darkgray", "coral")) +
  labs(y = expression(-log[10]("p-value")), x = expression(log[2]("fold change")))
