# LOH_Caller
## DESCRIPTION

A haplotype aware and cost-effective LOH calling solution, leveraging Python and CLC Genomics Workbench.

PROJECT OVERVIEW: Gradualism, the prevailing model for mutation accumulation and evolution, states that independent mutations occurring over extended periods contribute to genetic diversity. However, alternative modes of mutation accumulation, such as punctuated bursts, have also been proposed. In our study on Saccharomyces cerevisiae (baker's yeast), we identified punctuated bursts, specifically in the accumulation of loss-of-heterozygosity (LOH) tracts. LOH tracts occur when heterozygous genomic regions (e.g., Aa) become homozygous for either one of the two available alleles (i.e., AA or aa). Accurate and efficient identification of LOH event from whole genome sequencing (WGS) data is crucial for understanding these mutation bursts. Currently, our group uses Nexus Copy Number to detect LOH tracts. However, this program has limitations: (1) Nexus Copy Number lacks haplotype awareness, making it challenging to determine which allele an LOH clone tract became homozygous for and (2) its license is costly. The goal of our study is to develop a haplotype aware and cost-effective LOH calling solution, leveraging CLC Genomics Workbench 23.0.4 and Python 3.11.4 64-bit. 

## METHODOLOGY
![methods copy](https://github.com/JoyLove0/LOH_Finder/assets/108104001/20b85501-5da6-4908-a9d5-f001fba95c62)

# USAGE

## STEP 1: Running "Make_SNP_List.py"

This script should be ran after you have generates the variant list for the the diploid parent and the non reference parent (SK1) utilizing CLC Genomic Worknech. The workflow that does this process is called "Making Inputs for HetSNP List Generation."

Before running, there is user input section that needs to be filled out. Note that everything you (the user) write in this section MUST be put in quotes.

| Inputs          | Descriptions  |
| --------------  | ------------- |
| wd              | This is the working directory. This will contain your input VCF Files. |
| diploid_vcf     | Name of the vcf file that has the variant for the diploid. |
| SK1_vcf         | Name of the vcf file that has the variant for the the non reference parent (SK1). | 
| unreliables_csv | Name of the csv file containing where we cannot call LOH reliabile. This was created by going through the mapping files of SK1, S288c, and the SK1 X S288c diploid when mapped to S288c and finding areas of abnorally high and low coverage. This is included in the files provided. |   
| date | The date you run the code. It should be written as yearmonthdate. For example 20231221. This will be used to create the header in the resulting VCF file. |
| SNP_List_Name   | The name you wish to use for the resulting HetSNP VCF file. |

On the command line, the command is as follows:

`python Make_SNP_List.py`

Or run script in Python IDE (written in Spyder IDE 5.4.3).

