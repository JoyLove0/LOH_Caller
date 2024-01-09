# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 15:31:38 2024
@author: joylove5
"""
### Dependecies 
import pandas as pd
import numpy as np
############################### User Inputs ###################################
"""
INSTRUCTIONS: Everything you (the user) write in this section should be put in quotes EXECPT the coverage_minimum
"""
### Inputs
wd = 
clone_tables_list = 
coverage_minimum =  # Suggest 25

############################# Haplotype-Awareness ############################
def Make_haplotype_aware(clone_table, clone_name, path, minimum_coverage):
### Taking in Variants  
    variants = path + clone_table
    variant_df = pd.read_csv(variants)

### Building Clone Variant Table
# Making marker ID column
    for column in variant_df:
        variant_df["ID"] = variant_df["Chromosome"] + "_" + variant_df["Region"].astype(str)
# New Dataframe    
    new_df = pd.DataFrame()
    new_df["ID"] = variant_df["ID"]
    new_df["S288c_Parent_Allele"] = variant_df["Reference"]
    new_df["SK1_Parent_Allele"] = variant_df["Allele"]
    new_df[clone_name + "_Variant_Allele"] = variant_df["Allele"]
    new_df[clone_name + "_Variant_Coverage"] = variant_df["Coverage"]
    frequency = []
    frequency = (variant_df["Count"]/variant_df["Coverage"]) * 100
    frequency = np.round(frequency, decimals = 2)
    new_df[clone_name + "_Variant_Frequency"] = frequency
    
### Zygosity Determination
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
    new_df.to_csv(path + clone_name + "_Haplotype_Aware_Table.csv")
    
def main():
    ### Genotyping Loop
    file = open(clone_tables_list, "r+")
    for line in file:
        variant_table = line.strip()
        clone_name_key = line.strip().split('_')[0]
        output_table = Make_haplotype_aware(variant_table, clone_name_key, wd, coverage_minimum)
    
    return output_table

if __name__ == "__main__":
    main()