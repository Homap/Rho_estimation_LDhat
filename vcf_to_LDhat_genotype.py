#!/usr/bin/python3
import numpy as np
import textwrap
import argparse
from argparse import ArgumentParser, HelpFormatter


class RawFormatter(HelpFormatter):
	def _fill_text(self, text, width, indent):
		return "\n".join([textwrap.fill(line, width) for line in textwrap.indent(textwrap.dedent(text), indent).splitlines()])

program_descripton = f'''
	Title: vcf_to_genotype_fasta.py
	Version: v1.0
	Date: 2 Nov 2020
	Author: Homa Papoli Yazdi
	
	Description: Take a gzipped VCF file as an input, the region and the number of SNPs 
	to extract from VCF and outputs the LDhat inputs, the sites and locs files.

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

	For genotype/unphased data, the convention used is 0 and 1 for the two homozygotes,
	and 2 for heterozygotes. ? is used for missing data. The first line of the sites
	file is as follows: 4: number of sequences, 10: number of SNPs, 2: unphased 
	The first line of the locs file is as follows: 10: number of SNPs, 1200: total
	length of the sequence, L: using cross-over model (the other model is the gene
	coversion model)
	
	USAGE: 
	Make the program executable by:
	chmod u+x vcf_to_LDhat_genotype.py
	RUN:
	./vcf_to_LDhat_genotype.py file.vcf.gz sites.txt locs.txt superscaffold36 50
	'''

parser = ArgumentParser(description=program_descripton, formatter_class=RawFormatter)
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('--vcf', metavar='-v', required=True, nargs='+', type=argparse.FileType('r'), help='gzipped VCF')
requiredNamed.add_argument('--Nsnps', metavar='-n', required=True, type=int, help='number of SNPs')
requiredNamed.add_argument('--chr', metavar='-c', required=True, type=str, help='chromosome or scaffold')
requiredNamed.add_argument('--sites', metavar='-s', required=True, type=str, help='sites output file')
requiredNamed.add_argument('--locs', metavar='-l', required=True, type=str, help='locs output file')

args = parser.parse_args(['-h'])

# #**********************************************************************************
# # Open the vcf file and the output files
# #**********************************************************************************
# vcf_file = gzip.open(sys.argv[1], "rb")
# sites_out = open(sys.argv[2], "w")
# locs_out = open(sys.argv[3], "w")
# region = sys.argv[4]
# num_snps = int(sys.argv[5])

# #**********************************************************************************
# # Parse the genotype field
# genotype_matrix = []
# sample_names = []
# geno_position = []
# for line in vcf_file:
# 	line = line.rstrip().split()
# 	if line[0].startswith("#CHROM"):
# 		sample_names = line[9:]
# 	if not line[0].startswith("#"):
# 		if line[0] == "superscaffold36":
# 			geno_position.append(line[1])
# 			genotype_fields = line[9:]
# 			genotype_matrix.append(genotype_fields)

# for index, item in enumerate(genotype_matrix):
# 	if index < 50: #1000: # Take the first thousand SNPs
# 		#print(index)
# 		for index_ind, item_ind in enumerate(item):
# 			genotype = genotype_matrix[index][index_ind].split(":")[0]
# 			if genotype == "0/1":
# 				genotype_matrix[index][index_ind] = '2'
# 			elif genotype == "1/0":
# 				genotype_matrix[index][index_ind] = '2'
# 			elif genotype == "0/0":
# 				genotype_matrix[index][index_ind] = '0'
# 			elif genotype == "1/1":
# 				genotype_matrix[index][index_ind] = '1'
# 			elif genotype == "./.":
# 				genotype_matrix[index][index_ind] = '?'
# 			else:
# 				raise Exception("Genotype field contains incorrect notation!")
# 	else:
# 		break

# genotype_np_array = np.array([np.array(xi) for xi in genotype_matrix[0:50]])

# #print(genotype_matrix)
# #print(genotype_np_array)

# gen_array = genotype_np_array.transpose()

# gen_array_seq = np.apply_along_axis(lambda row: row.astype('|S1').tostring(), 
# 					axis=1,arr=gen_array)

# def chunks(s, n):
# 	for start in range(0, len(s), n):
# 		yield s[start:start+n]

# sites_out.write(str(gen_array.shape[0])+" "+str(gen_array.shape[1])+" "+"2"+"\n")
# for index, ind in enumerate(gen_array_seq):
# 	sites_out.write(">"+sample_names[index]+"\n")
# 	for seq in chunks(gen_array_seq[index], 50):
# 		sites_out.write(seq+"\n")

# # print(geno_position[0:1000])

# relative_geno = [str(int(pos) - int(geno_position[0]) + 1) for pos in geno_position]
# # print(relative_geno)

# L = int(geno_position[50]) - int(geno_position[0]) + 1

# # print(geno_position[1000] , geno_position[0] , L)

# locs_out.write(str(gen_array.shape[1])+" "+str(L)+" "+ "L"+"\n"+" ".join(relative_geno[0:50])+"\n")


# vcf_file.close()
# sites_out.close()
# locs_out.close()









