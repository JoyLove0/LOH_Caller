# LOH_Caller
## DESCRIPTION

This tool is a haplotype aware and cost-effective LOH calling solution, leveraging CLC Genomics Workbench 23.0.4 and Python 3.11.4 64-bit. The Python portion of this tool has both a command line interface (CLI) and a graphical user interface (GUI). The usage for both of these will be detailed below. 

Loss of Heterozygosity (LOH) is a class of chromosomal structural mutation that has important implications in the development of numerous genetic diseases, including cancer. LOH tracts occur when heterozygous genomic regions become homozygous for either one of the two available alleles. While using Saccharomyces cerevisiae as a model organism, we have identified LOH tracts in a mutation accumulation experiment. Accurate and efficient identification of LOH events from short-read whole genome sequencing (WGS) data is crucial for understanding these mutations.

## METHODOLOGY

First, we created a complete list of heterozygous single nucleotide polymorphisms (HetSNPs) in the parent strain. We validated this list by mapping reads from the fully heterozygous parent diploid to one of the haploid parents (S288c, reference genome). We then derived a list of variants shared by the two sets and removed variants present at repetitive regions of the genome to arrive at a high-confidence HetSNP list. This enabled us to generate a custom variant tract in CLC. Subsequently, we mapped WGS reads from the known LOH clones above to the reference, creating a variant list for each LOH clone. These mapping variants were interrogated against the high-confidence HetSNP variant tract. CLC provided allelic identity and frequency outputs at all HetSNP positions, which were used in a custom Python code to make haplotype-aware phased genotyping calls. These calls were validated against known LOH tracts initially identified in Nexus. Our initial results show that this approach can reveal the full complexity of LOH tracts, including discontinuities in gene conversion segments near the recombination endpoints.

![methods copy](https://github.com/JoyLove0/LOH_Finder/assets/108104001/20b85501-5da6-4908-a9d5-f001fba95c62)

# USAGE

## STEP 1 

## Running CLC Genomic Workbench Worflow called "Making Inputs for HetSNP List Generation 

![workflow snapshot](https://github.com/JoyLove0/LOH_Caller/assets/108104001/ba0c8a88-781f-4789-baf7-7d2ed7c3039c)

After hitting run, input the location of the raw reads from the diploid parent and the non-reference parent (SK1). The outputs of this workflow are the variant lists for both parents in VCF format.

## Python: Running "Make_SNP_List.py"

This step should be ran after you have generates the variant list for the diploid parent and the non-reference parent (SK1) utilizing CLC Genomic Worknech. The workflow that does this process is called "Making Inputs for HetSNP List Generation."

## Parameters
| Inputs          | Descriptions  |
| --------------  | ------------- |
| wd              | This is the working directory. This will contain your input VCF Files. |
| diploid_vcf     | Name of the vcf file that has the variant for the diploid. |
| SK1_vcf         | Name of the vcf file that has the variant for the the non reference parent (SK1). | 
| unreliables_csv | Name of the csv file containing where we cannot call LOH reliabile. This was created by going through the mapping files of SK1, S288c, and the SK1 X S288c diploid when mapped to S288c and finding areas of abnorally high and low coverage. This is included in the files provided. |   
| date | The date you run the code. It should be written as yearmonthdate. For example 20231221. This will be used to create the header in the resulting VCF file. |
| SNP_List_Name   | The name you wish to use for the resulting HetSNP VCF file. |

## CMI Interface
Before running, there is user input section that needs to be filled out. Note that everything you (the user) write in this section MUST be put in quotes.

On the command line, the command is as follows:

`python Make_SNP_List.py`

Or run script in Python IDE (written in Spyder IDE 5.4.3).

## STEP 2

## Running CLC Genomic Workbench Worflow called "Evaluation LOH Clone Variant"

![worflow snapshot 2](https://github.com/JoyLove0/LOH_Caller/assets/108104001/a793ad79-6d64-45d2-8539-8b66b18fdde9)

After hitting run, select the location of the raw sequencing data of the LOH clones.

## STEP 3 
## Python: Running "Make_Haplotype_Aware.py"

This script should be run after you have the "Evaluation LOH Clone Variant" worflow in CLC Genomic Workbench. It is reccomended that this code is ran in a new directory.

## Parameters
| Inputs            | Descriptions   |
| ----------------  | -------------- |
| wd                | This is the working directory. This will contain your input VCF/CSV Files. |
| clone_tables_list | This is the name of the txt file containing the list of clone table names. If you used the lines above this will be "Clone_Variant_Tables.txt". |
| coverage_minimum  | This must be a number with no quotesYou get to set the minimum coverage needed to reliable call region in your LOH tracks heterozygous or homozygous. If a position falls underneath this coverage, it will be marked "No Call Due to Low Coverage". Suggest to start with 25 and make this number higher to make this input more strict. | 

There are two version of this script that you will use depending on the file type of the LOH clone variant list. If you use the code ending in "_for_CSV.py" you will download each variant table from CLC. You must take out the collums in the right hand side so that only Chromosome, Region, Refernce, Allele, Count, and Covergae remain. If you are using code ending in "_for_VCF.py" the saves VCF files must be in the same directory of where the code itself is. 

NOTE: You must make a txt file containing all of the names of your LOH variant files. To do this, you can run the following in command line.

`ls *.csv > Clone_Variant_Tables.txt` (when you are using the code made for csv files)

OR 

`ls *.vcf > Clone_Variant_Tables.txt` (when you are using the code made for csv files)

Regardless, check the Clone_Variant_Tables.txt for correctness. An example file is provided. 

# Dependecies

This tool requires CLC Genomic Workbench. The specific tools used in the "Making Inputs for HetSNP Lisy" are Map Reads to Reference, Basic Variant Detection, Filter on Custom Criteria, and Export in VCF. The tools used in "Evaluating LOH Clone Variants" are Map Reads to Reference, Identify Known Mutations from Mapping, and Export as VCF. 

For more information on each these tools, see below:

| Program       | Manuals and/or Useful Links   |
| ------------- | ------------- |
| All Manuals | https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/current/index.php?manual=Introduction_CLC_Genomics_Workbench.html |
| Map Reads to Reference | https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/900/index.php?manual=Map_Reads_Reference.html |
| Basic Variant Detection | https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/900/index.php?manual=Basic_Variant_Detection.html |
| Filter on Custom Criteria | https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/current/index.php?manual=Filter_on_Custom_Criteria.html |
| Export in VCF | https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/current/index.php?manual=Export_in_VCF_format.html |
| Identify Known Mutations from Mapping | https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/950/index.php?manual=Identify_Known_Mutations_from_Sample_Mappings.html |
