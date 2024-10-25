# Question 1.1 
    conda activate qb24
    fastqc

## Click through the different metrics. Do any stand out as problematic? 
Both the per base sequence content and the per sequence GC content appear wonky. The per base sequence graph should have percentages of approximately 25% for each of the four bases, but the early positions in the graph range from approximately 17% to 42%. The per sequence GC content also has a much sharper and taller peak than what we'd expect. Lastly, 10 sequences in sample 1 and 7 sequences in sample 2 show up as being overrepresented. 

## Can you think of a reason why this sample does not match the expected GC content distribution and base content at the beginning of sequences?
The base distribution in the beginning could be chalked up to issues with the polymerase, as the beginning and end of a read is known to not be high quality. In terms of the GC content distribution, it might be that there are systematic issues, such as contamination of the samples. 

# Question 1.2


    Most overrepresented sequnce (27 counts in sample 1 and 25 counts in sample 2):  CTTGACCAAGATGAAACTGTTCGTATTCCTGGCCTTGGCCGTGGCCGCAA

## What is the origin of the most overrepresented sequence in this sample? Does it make sense?
The BLAST results show 100% identity with Drosophila melanogaster serine protease 3, serine protease 1 (SER1), and serine protease 2 (SER2) genes. Given that the samples came from Drosophila midguts and that serine proteases are involved in the digestion by cleaving peptide bonds in proteins, this result makes sense. 


# Question 2

## If you were to reject any samples with the percentage of unique reads less than 45%, how many samples would you keep?
There are 15 samples with percentages of unique reads over 45%. 

## Can you see the blocks of triplicates clearly? Try adjusting the min slider up. Does this suggest anything about consistency between replicates?
Increasing the min value using the slider allows the triplicates to be seen more clearly. This demonstrates consistency between replicates, as the triplicates appear to have similar sample to sample distances, which is as expected.

# Question 3.3

## Examine the PCA plot. Does everything look okay (We wouldn’t ask if it did)?
Yeah I'm going to be honest, I really thought I did something wrong when the graph was produced... Instead of having clear clusters (which shouldn't overlap), the "clusters" are mixed and/or merged on top of one another, so it's difficult to differentiate them (specifically LCF_Fe and Fe).

## What does the PCA plot suggest to you about the data? Why?
Based on the fact that LFC_Fe 3 appears to be in a cluster space Fe 2 and 3, and that Fe 1 appears to be in one with LFC_Fe 1 and 2, it appears as though the labels were switched or there is some error in the data/metadata.

# Question 3.6
## Do the categories of enrichments make sense? Why?
Broadly speaking they do. Because the samples were isolated from Drosophila midguts, we expect to see enrichments specific to the midgut and overall organism function. The enrichments match these predictions, with categories including carbohydrate catabolic processes, midgut development, nervous system development, etc. 