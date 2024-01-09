# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 09:44:22 2024
@author: joylove5

"""
### Dependecies 
import pandas as pd
import numpy as np
############################### User Inputs ###################################
"""
INSTRUCTIONS: Everything you (the user) write in this section should be put in quotes EXECPT the coverage_minimum
"""
### Inputs-
wd = 
clone_tables_list = 
coverage_minimum =  # Suggest 25

############################# Haplotype-Awareness ############################

def Make_haplotype_aware(clone_table, clone_name, path, minimum_coverage):
### Taking in Variants  
    variant_df = pd.read_csv(clone_table, sep = "\t", comment = "#", header=None)
    variant_df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", clone_name]
    variant_df["ID"] = variant_df["CHROM"] + "_" + variant_df["POS"].astype(str)
    print(variant_df)
### Building Clone Variant Table
# New Dataframe    
    new_df = pd.DataFrame()
    new_df["ID"] = variant_df["ID"]
    new_df["S288c_Parent_Allele"] = variant_df["REF"]
    new_df["SK1_Parent_Allele"] = variant_df["ALT"]
    new_df[clone_name + "_Variant_Allele"] = variant_df["ALT"]
# Build coverage and frequency collumns
    co_cov= []
    for word in variant_df[clone_name]:
        a = word.strip().split(',')[1]
        co_cov.append(a)
    count = []
    coverage = []
    for num in co_cov:
        first_num = num.strip().split(':')[0]
        count.append(int(first_num))
        second_num = num.strip().split(':')[1]
        coverage.append(int(second_num))
    new_df[clone_name + "_Variant_Coverage"] = coverage
    frequency = []
    for i in range(len(count)):
        try: 
            frequency.append((count[i] / coverage[i]) * 100)
        except ZeroDivisionError:
            frequency.append(0)
    frequency = np.round(frequency, decimals = 2)
    new_df[clone_name + "_Variant_Frequency"] = frequency

### Genotype Determination
    zygosity = []
    for percent in new_df[clone_name + "_Variant_Frequency"]:
        if percent > 89:
            zygosity.append("Homozygous for SK1")
        elif percent > 20 and percent < 80:
            zygosity.append("Heterozygous")
        elif percent < 11:
            zygosity.append("Homozygous for S288c")
        else:
            zygosity.append(None)
     
    for index, cov in enumerate(new_df[clone_name + "_Variant_Coverage"]): 
        if cov <= coverage_minimum:
            zygosity[index] = "No Call Due to Low Coverage"
    new_df["Genotype"] = zygosity

### Making Output CSV file
    new_df.to_csv(path + clone_name + "_Haplotype_Aware_Table(VCF).csv")

def main():
    ### Intergation Loop
    file = open(clone_tables_list, "r+")
    for line in file:
        line = line.strip('\n')
        variant_table = line.strip()
        clone_name_key = line.strip().split('_')[0]
        output_table = Make_haplotype_aware(variant_table, clone_name_key, wd, coverage_minimum)
    
    return output_table

if __name__ == "__main__":
    main()
    