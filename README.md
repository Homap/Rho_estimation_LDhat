# VCF to LDhat output

usage: vcf_to_LDhat_genotype.py [-h] vcf bed Nsnps chr Nwin  

Description: Take a gzipped VCF file as an input, the region and the number of  
SNPs to extract from VCF and outputs the LDhat inputs, the sites and locs files.  

**sites file**  
4 10 2  
>GenotypeA  
122110?000  
>GenotypeB  
1111201100  
>GenotypeC  
011111?112  
>GenotypeD  
2112210100  
  
**locs file**  
10 1200 L  
1 57 180 187 223 250 438 509 878 1034  
  
For genotype/unphased data, the convention used is 0 and 1 for the two  
homozygotes, and 2 for heterozygotes. ? is used for missing data. The first line of the  
sites file is as follows: 4: number of sequences, 10: number of SNPs, 2: unphased  
The first line of the locs file is as follows: 10: number of SNPs, 1200: total  
length of the sequence, L: using cross-over model (the other model is the gene  
coversion model)

It then processes the Ldhat outputs to produce plots and dataframes needed for downstream  
analysis.  

positional arguments:  
  vcf         gzipped VCF
  bed         scaffold length file in bed format
  Nsnps       number of SNPs
  chr         chromosome or scaffold
  Nwin        number of windows

optional arguments:
  -h, --help  show this help message and exit