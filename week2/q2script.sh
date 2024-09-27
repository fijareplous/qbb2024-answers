#!/bin/bash

#Header
echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt
#Go through each MAF file:
for MAF in 0.1 0.2 0.3 0.4 0.5
do 
    MAF_file=chr1_snps_${MAF}.bed
    #Coverage of whole chromosome:
    bedtools coverage -a genome_chr1.bed -b ${MAF_file} > total_coverage${MAF}.txt
    #Sum of SNPs:
    SNP_sum=$(awk '{s+=$4}END{print s}' total_coverage${MAF}.txt)
    #Sum of bases:
    base_sum=$(awk '{s+=$6}END{print s}' total_coverage${MAF}.txt) 
    #Background density:
    background=$(bc -l -e "${SNP_sum} / ${base_sum}")
    for feature in exons introns cCREs other 
    do 
        #Go through feature files:
        file=${feature}_chr1.bed
        #Check coverage by MAF and feature file:
        bedtools coverage -a ${file} -b ${MAF_file} > coverage_${MAF}_${feature}.txt
        #SNP sum for given feature:
        SNP_coverage=$(awk '{s+=$4}END{print s}' coverage_${MAF}_${feature}.txt)
        #Number of total bases in feature:
        feature_total=$(awk '{s+=$6}END{print s}' coverage_${MAF}_${feature}.txt)
        #Ratio of SNPs to feature background (density):
        density=$(bc -l -e "${SNP_coverage} / ${feature_total}")
        #Ratio of SNPs to genome background:
        enrichment=$(bc -l -e "${density} / ${background}")
        #Add results to text file:
        echo -e "${MAF}\t${feature}\t${enrichment}" >> snp_counts.txt
    done
done

    

