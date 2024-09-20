# Question 1.1: 
## How many 100bp reads are needed to sequence a 1Mbp genome to 3x coverage?
3 x 1 Mbp / 100 bp = 30,000 reads


# Question 1.4:
## In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
48539 / 1,000,000 = .049 or 4.9%

Approximately 4.9% of the genome has not been sequenced


## How well does this match Poisson expectations? How well does the normal distribution fit the data?
This matches the Poisson expectations well, as the line graph for the Poisson distribution lines up very neatly with the histogram of coverages. The normal distribution also matches fairly well, although the peak of the normal distribution is slightly to the right of the histogram peak. 


# Question 1.5:
## In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
128 / 1,000,000 = .000128 or .0128%


## How well does this match Poisson expectations? How well does the normal distribution fit the data?
This matches both the Poisson and normal distributions very well. The normal and Poisson distributions are more closely aligned with the data seen in the histogram as compared to what is seen in the 3x coverage. 


# Question 1.6:
## In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
14 / 1,000,000 = .000014 or .0014%


## How well does this match Poisson expectations? How well does the normal distribution fit the data?
This matches both the Poisson and normal distributions extremely well. Both the normal and Poisson distributions align very closely with the data seen in the histogram.


# Question 2.2-2.4:
    conda create -n graphviz -c conda-forge graphviz
    conda activate graphviz
    dot -Tpng ./edges.txt -o ex2_digraph.png 


# Question 2.5:
## Assume that the maximum number of occurrences of any 3-mer in the actual genome is five. Using your graph from Step 2.4, write one possible genome sequence that would produce these reads.
ATTCTTATTGATTCATTT

# Question 2.6:
## In a few sentences, what would it take to accurately reconstruct the sequence of the genome? 
To reconstruct the genome, you'd need read lengths that are significantly larger to ensure uniqueness between reads, as well as many more reads to be able to cover the entire genome. Analyzing significantly longer contigs would ensure that repeated regions do not interfere with accurate assembly, and more contigs would help to ensure that more of the genome is covered. 


