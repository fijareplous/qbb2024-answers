# Day 2 Lunch Answers

## Q1
Used the following code to visualize the data:

`head hg38-gene-metadata-feature.tsv | less -S`

Cut everything except for the gene_biotype (7th column)

`cut -f 7 hg38-gene-metadata-feature.tsv | head`

Got just the unique biotypes
`cut -f 7 hg38-gene-metadata-feature.tsv | sort | uniq`

Need to add a tally to determine how many of each biotype exist
`cut -f 7 hg38-gene-metadata-feature.tsv | sort | uniq -c`

There are 19618 protein coding genes. I'd like to learn more about unitary pseudogenes because I've never heard of them before and am intrigued by what they are. 

## Q2
Used the same format from the previous question but added sort after uniq to sort by the number of entries

`cut -f1 hg38-gene-metadata-go.tsv | sort | uniq -c | sort | tail -n1`

The ensemble gene ID ENSG00000168036 has the most go IDs (273). 

Saved file with rows just corresponding to ENSG00000168036
`grep "ENSG00000168036" hg38-gene-metadata-go.tsv > 168036.tsv`

Sorted by name_1006 (3rd column)
`sort -k 3 168036.tsv`

**** add what I think it does

## Q1 - Gencode
To just get the IG genes, not the pseudogenes, and determine how many in each chromosome:
`grep -w "IG_._gene" genes.gtf | cut -f 1 | uniq -c | sort -n`
Results:
   1 chr21
   6 chr16
  16 chr15
  48 chr22
  52 chr2
  91 chr14

To examine the distribution of pseudogenes in the same manner:
`grep -w "IG_._pseudogene" genes.gtf | cut -f 1 | uniq -c | sort -n`
Results:
   1 chr1
   1 chr10
   1 chr18
   1 chr8
   5 chr9
   6 chr15
   8 chr16
  45 chr2
  48 chr22
  83 chr14

The IG pseudogenes are more spread out across chromosomes that IG genes, with four chromosomes having one IG pseudogene but no IG genes. Chromosomes 2, 22, and 14 all have a large number of IG genes and pseudogenes.

## Q2 - Gencode
The reason why grep "pseudogene" genes.gtf wouldn't work to identify just the pseudogenes is because "pseudogene" appears in the tag section, as well as the gene_type column. One way to fix this would be to change it to the following:
`grep "gene_type \".*pseudogene\"; gene_name" genes.gtf | less -S`
The backslashes are to make sure that the quotations are taken as part of the string. Including gene_type in the beginning ensures that only instances where "pseudogene" appears in the gene_type category are selected. .* is used to search for 0+ characters (rather that 1+ with just *) in between the quotation mark and "pseudogene."


