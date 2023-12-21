# LOH_Finder
## DESCRIPTION

A haplotype aware and cost-effective LOH calling solution, leveraging Python and CLC Genomics Workbench.

PROJECT OVERVIEW: Gradualism, the prevailing model for mutation accumulation and evolution, states that independent mutations occurring over extended periods contribute to genetic diversity. However, alternative modes of mutation accumulation, such as punctuated bursts, have also been proposed. In our study on Saccharomyces cerevisiae (baker's yeast), we identified punctuated bursts, specifically in the accumulation of loss-of-heterozygosity (LOH) tracts. LOH tracts occur when heterozygous genomic regions (e.g., Aa) become homozygous for either one of the two available alleles (i.e., AA or aa). Accurate and efficient identification of LOH event from whole genome sequencing (WGS) data is crucial for understanding these mutation bursts. Currently, our group uses Nexus Copy Number to detect LOH tracts. However, this program has limitations: (1) Nexus Copy Number lacks haplotype awareness, making it challenging to determine which allele an LOH clone tract became homozygous for and (2) its license is costly. The goal of our study is to develop a haplotype aware and cost-effective LOH calling solution, leveraging Python and CLC Genomics Workbench.

## METHODOLOGY: 
![methods copy](https://github.com/JoyLove0/LOH_Finder/assets/108104001/20b85501-5da6-4908-a9d5-f001fba95c62)

### Usage 

## Step 1: Running "Make SNP List.py"

Before running, there is user input section that need to be filled out. Note that everything you (the user) write in this section MUST be put in quotes.

