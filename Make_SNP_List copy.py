# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18, 2023
Uploaded on Dec 21, 2023
@author: Joy Love
"""
### Dependecies 
import pandas as pd
############################### User Inputs ###################################
"""
INSTRUCTIONS: Everything you (the user) write in this section MUST be put in quotes
"""
### Inputs
wd = 
diploid_vcf = 
SK1_vcf = 
unreliables_csv = 
date = #yearmonthdate ex: 20231221

### Outputs
SNP_List_Name = 

############################# SNP List Generation ############################

def data_prep(filename):
    """
    Modifies the vcf files so that they are editable in the rest of the code.

    Parameters:
    ----------
    - filename : The names of the vcf files.

    Returns:
    -------
    - df : Dataframe with neccesary information from the vcf files.

    """
    name = filename.split(".vcf")[0]
    df = pd.read_csv(filename, sep = "\t", comment = "#")
    df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", name]
    df["ID"] = df["CHROM"] + "_" + df["POS"].astype(str)
    return df

def semi_join(data_Set1, data_Set2, joiner):
    """
    Combines two dataset on a common collumn in the table. The semi-join means 
    that one table wiil be match to another and anything that does not match the
    first table is discarded. 
    
    Parameters:
    ----------
    - data_Set1 : The dataset that must be matched.
    - data_Set2 : Another dataset. We will only take the values out of this table that match the first dataset
    - joiner : The column name that the dataset will be matched on. Must be present in both tables. 
        
    Returns:
    -------
    A new dataset containing all of the values in dataset 2 that match dataset 1.
    """
    # Set Up   
    semi = data_Set1.merge(data_Set2,on=joiner)
    data_Set1[joiner].isin(data_Set2[joiner])
    semi=data_Set1.merge(data_Set2,on=joiner)
    # Left Semi-join
    new_semi = data_Set1[data_Set1[joiner].isin(semi[joiner])]
    pd.DataFrame(new_semi)

    # Cleaning 
    #new_semi.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "HetSNP_List_Variants"]
    return(new_semi)

def filter_positions(df1, unreliable_ranges):
    """
    Filter positions in the dataframed that contain position in both the diploid
    and in Sk1 that fall within the ranges specified in unreliable_df.

    Parameters:
    -------
    - match_df: DataFrame containing positions to be filtered.
    - unreliable_df: DataFrame containing unreliable ranges.

    Returns:
    -------
    - filtered_df: DataFrame with positions filtered based on unreliable ranges.
    """

    # Create a boolean mask for positions to exclude
    exclude_mask = False
    for index, row in unreliable_ranges.iterrows():
        chromosome = row['Chromosome']
        start_position = row['Region_Start']
        end_position = row['Region_End']
        
        # Create a mask for positions falling within the unreliable range
        mask = (df1['CHROM'] == f'chr{chromosome}') & (df1['POS'].between(start_position, end_position))
        
        # Update the exclude mask with the current mask
        exclude_mask |= mask

    # Filter dip_df using the exclude mask
    filtered_df = df1[~exclude_mask]

    return filtered_df

def add_header(_filtered_df, date):
    """
    Adds a header to a VCF file and appends the content of _filtered_df to it.

    Parameters:
    ----------
    - _filtered_df : DataFrame containing the filtered positions.
    - date : A string containing the date for the VCF file header. This comes
    from the user input

    Returns:
    -------
    None as it writes to a vcf file.
    """
    
    header = """##fileformat=VCFv4.1
##fileDate=""" + date + """
##source=CLC Genomics Workbench 23.0.4 build 20230515014922
##reference=/CLC_Ddrive/S288c reference genome/Reference (Genome)
##contig=<ID=chr1,length=230218>
##contig=<ID=chr2,length=813184>
##contig=<ID=chr3,length=316620>
##contig=<ID=chr4,length=1531933>
##contig=<ID=chr5,length=576874>
##contig=<ID=chr6,length=270161>
##contig=<ID=chr7,length=1090940>
##contig=<ID=chr8,length=562643>
##contig=<ID=chr9,length=439888>
##contig=<ID=chr10,length=745751>
##contig=<ID=chr11,length=666816>
##contig=<ID=chr12,length=1078177>
##contig=<ID=chr13,length=924431>
##contig=<ID=chr14,length=784333>
##contig=<ID=chr15,length=1091291>
##contig=<ID=chr16,length=948066>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total number of filtered reads used to call variants at this position in this sample.">
##FORMAT=<ID=CLCAD2,Number=R,Type=Integer,Description="Allele depth, number of filtered reads supporting the alleles, ordered as listed in REF and ALT.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HetSNP_Variants
"""

    output_VCF = wd + SNP_List_Name
    with open(output_VCF, 'a') as vcf:
        vcf.write(header)

    # Append the content of _filtered_df without header
    final_ouput = _filtered_df.to_csv(output_VCF, sep="\t", mode='a', header=False, index=False)
    return final_ouput

def main():
    """
    First, this function matches the SK1's SNPs to the Parent Diploid's SNP. 
    This will then be used to generate a confident HetSNP List, to identiy these
    known variants in the Diploid clones. 

    Returns:
    -------
    The Confident HetSNP List.

    """
    # Data Prep
    dip_df = data_prep(diploid_vcf)
    sk1_df = data_prep(SK1_vcf)
    
    # Semijoin and cleaning
    match_df = semi_join(dip_df, sk1_df, "ID")
    
    # Filrering out unreliable coordinates
    unreliable_df = pd.read_csv(unreliables_csv)
    filtered_df = filter_positions(match_df, unreliable_df)
    
    # Adding header and export csv
    output_SNP_VCF = add_header(filtered_df, date)
    return output_SNP_VCF
    #Checks
    #print("Diploid Dataframe: ")
    #print(dip_df)
    #print("Sk1 Dataframe: ")
    #print(sk1_df)
    
    #print("Matched DataFrame:")
    #print(match_df)
    
    #print("Unreliable DataFrame:")
    #print(unreliable_df)
    
    #print("Filtered DataFrame:")
    #print(filtered_df)
    
    #print("HetSNP List VCF")
    #print(output_SNP_VCF) #should be none. Check your working directory for the SNP list.
    #print("Should be none. Check your working directory for the SNP list.")
if __name__ == "__main__":
    main()
    
