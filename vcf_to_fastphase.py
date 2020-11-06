#!/usr/bin/python
import sys
import gzip
import numpy as np

"""
Produce input file for fastPHASE input
no.individuals
no.SNPsites
P pos(1) pos(2) ... pos(no.SNPsites) <optional line>
SSS...SSS <optional line>
ID (1)
genotypes(1-a)
genotypes(1-b)
ID (2)
genotypes(2-a)
genotypes(2-b)
.
.
.
ID (no.individuals)
genotypes(no.individuals-a)
genotypes(no.individuals-b)

Example:
3
4
# id 1
1a11
0t01
# id 2
1t11
0a00
# id 3
?a01
?t10
"""

# Read a vcf file
vcf_file = gzip.open(sys.argv[1], "rb")
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

#print(genotype_matrix)
#print(sample_names)
#print(geno_position)


for index, item in enumerate(genotype_matrix):
	for index_ind, item_ind in enumerate(item):
		genotype_matrix[index][index_ind] = genotype_matrix[index][index_ind].split(":")[0]

genotype_np_array = np.array([np.array(xi) for xi in genotype_matrix])
#print(genotype_np_array)
#print(genotype_np_array.shape)

transposed_matrix = genotype_np_array.transpose()

#print(transposed_matrix)
#print(transposed_matrix.shape)

allele_ind = []
for ind in transposed_matrix:
	allele1 = []
	allele2 = []
	for geno in ind:
		allele1.append(geno.split("/")[0])
		allele2.append(geno.split("/")[1])
	allele_ind.append([allele1, allele2])

allele_ind_array = np.array([np.array(xi) for xi in allele_ind])
#print(allele_ind_array)

print(str(transposed_matrix.shape[0])+"\n"+str(transposed_matrix.shape[1]))
for ind_index, ind_item in enumerate(allele_ind):
	print("#"+" "+sample_names[ind_index])
	for hap in ind_item:
		print("".join(hap))


vcf_file.close()
