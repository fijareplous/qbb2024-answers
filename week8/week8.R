BiocManager::install("zellkonverter")
library(zellkonverter)
library(scuttle)
library(scater)
library(scran)

#Set working directory
setwd("~/qbb2024-answers/week8")


#Part 1


#Create SingleCellExperiment object
gut <- readH5AD("v2_fca_biohub_gut_10x_raw.h5ad")

#Change assay name to "counts"
assayNames(gut) <- "counts"

#Normalize counts
gut <- logNormCounts(gut)

#Inspect gut object
gut
#Q1: There are 13,407 genes being quantified across 11,788 cells. X_pca, X_tsne, and
#X_umap are the three dimension reduction datasets present.

#Look at metadata
as.data.frame(colData(gut))
colnames(colData(gut))
#Q2: There are 39 column names
#I think age, sex, and n_counts could all be interesting, as they provide insight
#into commonly expressed genes among different age and sex groups


#Part 2

#Sum expression of each gene across cells
genecounts <- rowSums(assay(gut))

#Summarize genecounts
summary(genecounts)
head( sort(genecounts, decreasing=TRUE ))
#Q3: the mean count is 3,185, and the median is 254, meaning that the data is heavily
#skewed right (centered at lower gene counts). The three genes with the highest 
#expression levels are lncRNA:Hsromega, pre-rRNA:CR45845, and lncRNA:roX1. All three of
#these genes are non-coding RNAs. 

#Create cellcounts vector
cellcounts <- colSums(assay(gut))

#Plot cellcounts as histogram
hist(cellcounts)

#Find mean from summary
summary(cellcounts)
#Q4a: The mean cell counts is 3,622. Counts that are higher (>10,000) likely have much
#higher basal expression comparatively, perhaps due to cell function.

#Create celldetected vector
celldetected <- colSums(assay(gut)>0)

#Plot celldetected as histogram
hist(celldetected)

#Find mean from summary
summary(celldetected)
#Q4b: The mean number of genes detected per cell is 1,059. Given that there are 13,407
#genes total, the calculated mean represents approximately 7.9% of the total genes. 

#Create vector of mitochondrial gene names
mito <- grep("^mt:", rownames(gut), value = TRUE)

#Create DataFrame for mito
df <- perCellQCMetrics(gut, subsets = list(Mito = mito))

#Convert to data.frame
df <- as.data.frame(df)

#Confirm mean and sum match previous calculations
summary(df)

#Add metrics to cell metadata
colData(gut) <- cbind( colData(gut), df )

#Plot percent of reads from mitochondria
png('week8_mito_plot.png')
plotColData(gut, y = "subsets_Mito_percent", x = "broad_annotation") +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
#Q5: Epithelial cells seem to have the highest proportion of mitochondrial gene
#expression. This may be due to the fact that epithelial cells have a high abundance, 
#as they cover the skin and other organs, as well as the fact that they have a rapid
#turnover rate. 


#Part 3

#Vector for cells of interest relating to epithelial cells
coi <- colData(gut)$broad_annotation == "epithelial cell"

#New SingleCellExperiment object for coi
epi <- gut[, coi]

#Plot epi according to X_umap; color by annotation
png('week8_Xumap_plot.png')
plotReducedDim(epi, "X_umap", colour_by="annotation")
dev.off()

#List with pairwise comparisons between all annotation categories
marker.info <- scoreMarkers(epi, colData(epi)$annotation)

#ID top marker genes in the anterior midgut 
chosen <- marker.info[["enterocyte of anterior adult midgut epithelium"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4])
#Q6b: The top 6 marker genes in the anterior midgut are Mal-A6, Men-b, vnd, betaTry,
#Mal-A1, and Nhe2. At least 4 of the markers have been identified
#as being involved in glucose/carbohydrate digestion, which means that the anterior
#midgut appears to specialize in metabolizing carbohydrates

#Notes on marker function (not to be graded, just for my own sanity lol):
#1) Carbohydrate metabolic process
#2) malate dehydrogenase (decarboxylating) (NADP+) activity. Involved in glucose homeostasis.
#3) transcription factor in neuroectoderm patterning, formation and specification of ventral neuroblasts
#4) putative digestive enzyme with predicted serine-type endopeptidase activity, midgut 
#5) carbohydrate metabolic process
#6) Na[+]/H[+] hydrogen exchanger 2, increases intracellular pH

#Plot expression of Mal-A6 (top marker gene)
png("week8_Mal-A6_plot.png")
plotExpression(gut, "Mal-A6", x="annotation") +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

#Now to repeat for somatic precursor cells:

#Vector for cells of interest relating to somatic precursor cells
spc_coi <- colData(gut)$broad_annotation == "somatic precursor cell"

#New SingleCellExperiment object for coi (somatic precursor cells)
spc <- gut[, spc_all]

#List with pairwise comparisons between all annotation categories
spc_marker <- scoreMarkers(spc, colData(spc)$annotation)

#ID top marker genes in the somatic precursor cells
chosen_spc <- marker.info[["intestinal stem cell"]]
ordered_spc <- selected_spc[order(selected_spc$mean.AUC, decreasing=TRUE),]
head(ordered_spc[,1:4])
#The top marker genes are hdc, kek5, N, zfh2, Tet, and Dl. These genes appear to play
#important roles in signaling, transcription, and development, which make sense for them to be 
#identified in stem cells.

#Notes on marker function (not to be graded, just for my own sanity lol):
#1) Histidine decarboxylase
#2) transmembrane protein involved in BMP signaling regulation
#3) Essential signaling protein, involved in development
#4) putative transcription factor, cell fate and apoptosis during development
#5) demethylates DNA methylated on the 6th position of adenine (N(6)-methyladenosine) DNA
#6) No hit/wouldn't load

#Create vector with names of the top six genes of interest
spc_goi <- rownames(ordered_spc)[1:6]

#Plot expression of top six marker genes across cell types
png("week8_6_genes_plot.png")
plotExpression(spc, features = spc_goi, x = "annotation") +
  theme(axis.text.x=element_text(angle=90)) 
dev.off()
#Enteroblasts and intestinal stem cells have more similar expression based on the graph.
#Dl appears to be most specific for intestinal stem cells. 
