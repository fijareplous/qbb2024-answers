library(dplyr)
library(tidyverse)
library(broom)

relevant_expr <- read_tsv("~/Downloads/dicts_expr.tsv")

mutated <- relevant_expr %>% dplyr::mutate(Tissue_Gene=paste0(Tissue, " ",
   GeneID)) %>% mutate(log_mutant = log2(Expr + 1))

ggplot(data=mutated, mapping=aes(x=Tissue_Gene, y=log_mutant)) + coord_flip() +
  geom_violin() + xlab("Tissue, Gene Pair") + ylab("Expression Value")

'
I am slightly surprised by the high degree of variability in terms of gene expression
across tissues, although these results do make sense. For instance, there is more 
gene expression variability seen in the stomach and small intestines, which makes
sense given that these tissues are likely responding to environmental and changing
their expression patterns accordingly. The testis, on the other hand, would be 
considered highly specialized and thus see less variation, which is reflected in
the graph.
'