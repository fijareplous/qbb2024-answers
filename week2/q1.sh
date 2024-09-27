#!/bin/bash

#Question 1.4
bedtools sort -i genes.bed > genes_sorted.bed
bedtools merge -i genes_sorted.bed > genes_chr1.bed

bedtools sort -i exons.bed > exons_sorted.bed
bedtools merge -i exons_sorted.bed > exons_chr1.bed

bedtools sort -i cCREs.bed > cCREs_sorted.bed
bedtools merge -i cCREs_sorted.bed > cCREs_chr1.bed

#Question 1.5
bedtools subtract -a genes_chr1.bed -b exons_chr1.bed > introns_chr1.bed

#Question 1.6
bedtools subtract -a genome_chr1.bed -b exons_chr1.bed > genome_minus_exons.bed
bedtools subtract -a genome_minus_exons.bed -b introns_chr1.bed > genome_no_ons.bed
bedtools subtract -a genome_no_ons.bed -b cCREs_chr1.bed > other_chr1.b