#!/usr/bin/python
import sys
import gzip
import numpy as np

'''
Title: vcf_to_genotype_fasta.py
Date: 2 Nov 2020
Author: Homa Papoli Yazdi

# Description:

# Produce the sites and locs file for the LDhat input


# 4 10 2
# >GenotypeA
# 122110?000
# >GenotypeB
# 1111201100


For genotype/unphased data, the convention used is 0 and 1 for the two homozygotes,
and 2 for heterozygotes.

Usage: 
'''

# Read a vcf file
vcf_file = gzip.open(sys.argv[1], "rb")
sites_out = open(sys.argv[2], "w")
locs_out = open(sys.argv[3], "w")

# Parse the genotype field
genotype_matrix = []
sample_names = []
geno_position = []
for line in vcf_file:
	line = line.rstrip().split()
	if line[0].startswith("#CHROM"):
		sample_names = line[9:]
	if not line[0].startswith("#"):
		if line[0] == "superscaffold36":
			geno_position.append(line[1])
			genotype_fields = line[9:]
			genotype_matrix.append(genotype_fields)

#print(sample_names)
# print(geno_position)
# print(geno_position[0])
# print(geno_position[1999])
# 4522739 - 3525136 + 1 = 997604
#print(genotype_matrix)

for index, item in enumerate(genotype_matrix):
	if index < 1000: # Take the first thousand SNPs
		#print(index)
		for index_ind, item_ind in enumerate(item):
			genotype = genotype_matrix[index][index_ind].split(":")[0]
			if genotype == "0/1":
				genotype_matrix[index][index_ind] = '2'
			elif genotype == "1/0":
				genotype_matrix[index][index_ind] = '2'
			elif genotype == "0/0":
				genotype_matrix[index][index_ind] = '0'
			elif genotype == "1/1":
				genotype_matrix[index][index_ind] = '1'
			elif genotype == "./.":
				genotype_matrix[index][index_ind] = '?'
			else:
				raise Exception("Genotype field contains incorrect notation!")
	else:
		break

genotype_np_array = np.array([np.array(xi) for xi in genotype_matrix[0:1000]])

#print(genotype_matrix)
#print(genotype_np_array)

gen_array = genotype_np_array.transpose()

gen_array_seq = np.apply_along_axis(lambda row: row.astype('|S1').tostring(), 
					axis=1,arr=gen_array)

def chunks(s, n):
	for start in range(0, len(s), n):
		yield s[start:start+n]

sites_out.write(str(gen_array.shape[0])+" "+str(gen_array.shape[1])+" "+"2"+"\n")
for index, ind in enumerate(gen_array_seq):
	sites_out.write(">"+sample_names[index]+"\n")
	for seq in chunks(gen_array_seq[index], 1000):
		sites_out.write(seq+"\n")

# print(geno_position[0:1000])

relative_geno = [str(int(pos) - int(geno_position[0]) + 1) for pos in geno_position]
# print(relative_geno)

L = int(geno_position[1000]) - int(geno_position[0]) + 1

# print(geno_position[1000] , geno_position[0] , L)

locs_out.write(str(gen_array.shape[1])+" "+str(L)+" "+ "L"+"\n"+" ".join(relative_geno[0:1000])+"\n")


vcf_file.close()
sites_out.close()
locs_out.close()









