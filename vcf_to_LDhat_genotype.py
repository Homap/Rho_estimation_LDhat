#!/usr/bin/python3
import numpy as np
import gzip
import io
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
	'''

parser = ArgumentParser(description=program_descripton, formatter_class=RawFormatter)
#requiredNamed = parser.add_argument_group('required named arguments')
parser.add_argument('vcf', help='gzipped VCF')
parser.add_argument('Nsnps', help='number of SNPs', type=int)
parser.add_argument('chr', help='chromosome or scaffold', type=str)
parser.add_argument('sites', help='sites output file', type=str)
parser.add_argument('locs', help='locs output file', type=str)

args = parser.parse_args()

#**********************************************************************************
def chunks(s, n):
	for start in range(0, len(s), n):
		yield s[start:start+n]
#**********************************************************************************
# Parse the genotype field
genotype_matrix = []
sample_names = []
geno_position = []
with io.TextIOWrapper(gzip.open(args.vcf, 'r')) as vcf_file:
	for line in vcf_file:
		line = line.rstrip().split()
		if line[0].startswith("#CHROM"):
			sample_names = line[9:]
		if not line[0].startswith("#"):
			if line[0] == args.chr:
				geno_position.append(line[1])
				genotype_fields = line[9:]
				genotype_matrix.append(genotype_fields)

for index, item in enumerate(genotype_matrix):
	if index < args.Nsnps: 
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

genotype_np_array = np.array([np.array(xi) for xi in genotype_matrix[0:args.Nsnps]])


gen_array = genotype_np_array.transpose()

gen_array_seq = np.apply_along_axis(lambda row: row.astype('|S1').tostring(), axis=1,arr=gen_array)



with open(args.sites, 'w') as sites_out:
	sites_out.write(str(gen_array.shape[0])+" "+str(gen_array.shape[1])+" "+"2"+"\n")
	for index, ind in enumerate(gen_array_seq):
		sites_out.write(">"+sample_names[index]+"\n")
		for seq in chunks(gen_array_seq[index], args.Nsnps):
			sites_out.write(seq+"\n")

relative_geno = [str(int(pos) - int(geno_position[0]) + 1) for pos in geno_position]

L = int(geno_position[args.Nsnps]) - int(geno_position[0]) + 1
with open(args.locs, 'w') as locs_out:
	locs_out.write(str(gen_array.shape[1])+" "+str(L)+" "+ "L"+"\n"+" ".join(relative_geno[0:args.Nsnps])+"\n")










