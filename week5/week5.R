library(DESeq2)
library(vsn)
library(matrixStats)
library(readr)
library(tibble)
library(hexbin)
library(ggfortify)
library(tidyverse)


#Part 3.1

#Load file salmon.merged.gene_counts.tsv
salmon = readr::read_tsv("~/qbb2024-answers/week5/salmon.merged.gene_counts.tsv")

#Use the column gene_name as row names 
salmon = column_to_rownames(salmon, var = "gene_name")

#Remove column gene_id
salmon = salmon %>% dplyr::select(-gene_id)

#Convert data from each column from numeric to integer
salmon = salmon %>% mutate_if(is.numeric, as.integer)

#Keep only rows with at least 100 reads
salmon = salmon[rowSums(salmon) > 100,]

#Select only the “narrow region” samples
narrow = salmon %>% select("A1_Rep1":"P2-4_Rep3")



#Part 3.2

#Create metadata tibble with two columns (sample names and replicate number)
narrow_metadata = tibble(tissue=as.factor(c("A1", "A1", "A1", "A2_3", "A2_3", 
    "A2_3", "Cu", "Cu", "Cu", "LFC_Fe", "LFC_Fe", "Fe", "LFC_Fe", "Fe", "Fe", 
    "P1", "P1", "P1", "P2_4", "P2_4", "P2_4")),
    rep=as.factor(c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 1, 3, 2, 3, 1, 2, 3, 1, 2, 3)))

#create a DESeq2 object using data and  metadata tibble
narrowdata = DESeqDataSetFromMatrix(countData=as.matrix(narrow), colData=narrow_metadata, design=~tissue)

#Correct for batch-effects via variance stabilizing transformation (vst)
narrowVstdata = vst(narrowdata)

#Plot the data!
meanSdPlot(assay(narrowVstdata))



#Part 3.3

#Create PCA data
narrowPcaData = plotPCA(narrowVstdata,intgroup=c("rep","tissue"), returnData=TRUE)

#Plot PCA data using ggplot (and save)
ggplot(narrowPcaData, aes(PC1, PC2, color=tissue, shape=rep)) + geom_point(size=5) +
  ggtitle("PCA Plot of Drosophila Midgut Tissues, Corrected Via Vst")
ggsave("~/qbb2024-answers/week5/week5_PCA_graph.png")
#Graph looks... not great. More specifically, LFC_Fe 3 and Fe 1 should be swapped
#To fix this, I manually switched them in the metadata


#Part 3.4

#Convert vst-corrected data into matrix
narrowVstdata = as.matrix(assay(narrowVstdata))

#Average each set of replicates
combined = narrowVstdata[,seq(1, 21, 3)]
combined = combined + narrowVstdata[,seq(2, 21, 3)]
combined = combined + narrowVstdata[,seq(3, 21, 3)]
combined = combined / 3
sds = rowSds(combined)

#Filter out genes with sds lower than 1
filt = rowSds(combined) > 1
narrowVstdata = narrowVstdata[filt,]



#Part 3.5

# Set seed for reproducibility
set.seed(42)

#Cluster genes using k-means 
k=kmeans(narrowVstdata, centers=12)$cluster

#Get ordering of samples to put them into clusters
ordering = order(k)

#Reorder genes
k = k[ordering]

# Plot heatmap of expressions & clusters; save
png("~/qbb2024-answers/week5/week5_heatmap.png")
heatmap(narrowVstdata[ordering,],Rowv=NA,Colv=NA,RowSideColors = RColorBrewer::brewer.pal(12,"Paired")[k])
dev.off()


#Part 3.6

#Get genes from first cluster
genes = rownames(narrowVstdata[k == 1,])

#Save these to text file to import into Panther
write.table(genes, "~/qbb2024-answers/week5/week5_cluster_genes.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)


