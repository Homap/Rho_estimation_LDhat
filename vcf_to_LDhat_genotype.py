#!/usr/bin/python3
import numpy as np
import gzip
import io
import textwrap
import argparse
from argparse import ArgumentParser, HelpFormatter
import subprocess
from subprocess import Popen, PIPE, STDOUT
import pandas as pd
import matplotlib.pyplot as plt
import pexpect
import sys
from matplotlib.animation import FuncAnimation
# import seaborn as sns

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
parser.add_argument('vcf', help='gzipped VCF')
parser.add_argument('bed', help='scaffold length file in bed format')
parser.add_argument('Nsnps', help='number of SNPs', type=int)
parser.add_argument('chr', help='chromosome or scaffold', type=str)
parser.add_argument('Nwin', help='number of windows', type=int)

args = parser.parse_args()

#**********************************************************************************
def slidingWindow(sequence_l, winSize, step):
	""" Returns a generator that will iterate through
	the defined chunks of input sequence. Input
	sequence must be iterable."""

	# Verify the inputs
	if not ((type(winSize) == type(0)) and (type(step) == type(0))):
		raise Exception("**ERROR** type(winSize) and type(step) must be int.")
	if step > winSize:
		raise Exception("**ERROR** step cannot be larger than winSize.")
	if winSize > sequence_l:
		pass
	# Pre-compute number of chunks to emit
	numOfChunks = ((int(sequence_l-winSize)/step))+1
	numOfChunks = int(numOfChunks) 

	for i in range(0, numOfChunks*step, step):
		yield i,i+winSize

def chunks(s, n):
	for start in range(0, len(s), n):
		yield s[start:start+n]
#**********************************************************************************
# Read scaffold length file into a dictionary
scaf_len = {}
with open(args.bed, 'r') as lenfile:
	next(lenfile)
	for line in lenfile:
		line = line.rstrip().split()
		scaf_len[line[0]] = int(line[2]) - int(line[1])
# print(scaf_len)
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

# relative_geno = [str(int(pos) - int(geno_position[0]) + 1) for pos in geno_position]
# print(relative_geno)
# Create list of windows
windows = slidingWindow(len(genotype_matrix), args.Nsnps, args.Nsnps-5)
# for i in windows: print(i)

interval_counter = 0
for interval_index, interval in enumerate(windows):
	# print(interval)
	if interval_index < args.Nwin:
		start = list(interval)[0]
		end = list(interval)[1]
		# print(start, end)
		for index, item in enumerate(genotype_matrix):
			index = index + start
			# print("INDEX", index)
			if start <= index <= end: 
				for index_ind, item_ind in enumerate(item):
					genotype = genotype_matrix[index][index_ind].split(":")[0]
					if genotype == "0/1" or genotype == "2":
						genotype_matrix[index][index_ind] = '2'
					elif genotype == "1/0" or genotype == "2":
						genotype_matrix[index][index_ind] = '2'
					elif genotype == "0/0" or genotype == "0":
						genotype_matrix[index][index_ind] = '0'
					elif genotype == "1/1" or genotype == "1":
						genotype_matrix[index][index_ind] = '1'
					elif genotype == "./.":
						genotype_matrix[index][index_ind] = '?'
					else:
						raise Exception("Genotype field contains incorrect notation!")
			else:
				break

		genotype_np_array = np.array([np.array(xi) for xi in genotype_matrix[start:end]])
		# print("1", genotype_np_array)

		gen_array = genotype_np_array.transpose()
		# print("2", gen_array)

		gen_array_seq = np.apply_along_axis(lambda row: row.astype('|S1').tostring().decode('utf-8'), axis=1,arr=gen_array)
		# print("3", gen_array_seq)

		interval_counter += 1
		sites = args.chr + ":" + str(interval_counter) + ".sites.txt"
		# print(outname)
		with open(sites, 'w') as sites_out:
			sites_out.write(str(gen_array.shape[0])+" "+str(gen_array.shape[1])+" "+"2"+"\n")
			for index, ind in enumerate(gen_array_seq):
				sites_out.write(">"+sample_names[index]+"\n")
				for seq in chunks(gen_array_seq[index], args.Nsnps):
					sites_out.write(seq+"\n")

		L = int(geno_position[end]) - int(geno_position[start]) + 1
		# print(L)
		locs = args.chr + ":" + str(interval_counter) + ".locs.txt"
		with open(locs, 'w') as locs_out:
			new_coord = [str(int(coord) - int(geno_position[start]) + 1) for coord in geno_position[start:end]]
			locs_out.write(str(gen_array.shape[1])+" "+str(L)+" "+ "L"+"\n"+"\n".join(new_coord)+"\n")

		original_pos = args.chr + ":" + str(interval_counter) + ".pos.txt"
		with open(original_pos, 'w') as pos_out:
			pos_out.write(str(gen_array.shape[1])+" "+str(L)+" "+ "L"+"\n"+"\n".join(geno_position[start:end])+"\n")

		print("Running LDhat")
		ldhat_out = args.chr + ":" + str(interval_counter) + ".ldhat."
		# cmd = ["./pairwise", "-seq", sites, "-loc", locs, "-lk", "ostrich_genotypenew_lk.txt", "-prefix", ldhat_out]
		# # print(cmd)
		# p = subprocess.Popen(cmd).communicate()
		ldhat_cmd = pexpect.spawn('./pairwise', ["-seq", sites, "-loc", locs, "-lk", "ostrich_genotypenew_lk.txt", "-prefix", ldhat_out])
		# # print(sshAP)
		# sshAP.logfile = sys.stdout#.buffer

		# Add a timeout in case the script fails
		ldhat_cmd.timeout = 60
		print("running grid size question")
		ldhat_cmd.expect('Do you wish to change grid over which to estimate likelihoods')
		ldhat_cmd.sendline('0')

		print("running sliding window question")
		ldhat_cmd.expect('Do you wish to carry out a sliding windows analysis')
		ldhat_cmd.sendline('0')

		print("running rmin full table")
		ldhat_cmd.expect('Full table')
		ldhat_cmd.sendline('2')

		print("running 4Ner method of moment question")
		ldhat_cmd.expect('Estimate 4Ner by moment method')
		ldhat_cmd.sendline('1')

		print("Plotting Composite-likelihood as a function of 4Ner")

		composite_out = args.chr + ":" + str(interval_counter) + ".ldhat." + "outfile.txt"
		composite_png = args.chr + "_" + str(interval_counter) + ".png"
		outfile = pd.read_table(composite_out, \
		skip_blank_lines=True, skipinitialspace=True, sep='\s+',\
		skiprows=lambda x: x in [0, 1, 2, 3, 4, 5])

		plt.figure()
		plt.plot(outfile['4Ner(region)'], outfile['Pairwise'])
		plt.xlabel('4Ner (region)')
		plt.ylabel('Composite-likelihood')
		plt.savefig(composite_png)

		f = open(composite_out, "r")
		lk = ""
		for line in f:
			if line.startswith("Maximum"):
				line = line.rstrip()
				lk = lk + line		
		print("chr"+"\t"+"SNP1_pos"+"\t"+"SNP2_pos"+"\t"+"rho"+"\t"+"Lk")
		print(args.chr, geno_position[start], geno_position[end], lk.split()[4], lk.split()[8])

		rmin_out = args.chr + ":" + str(interval_counter) + ".ldhat." + "rmin.txt"
		rmin_png = args.chr + "_" + str(interval_counter) + ".png"
		pairwise_rmin_out = args.chr + "_" + str(interval_counter) + ".pairwise.rmin.txt"
		rmin_file = open(rmin_out, "r")
		rmin = rmin_file.readlines()[5:]
		# Create a list of lower diagonal
		lower_diagonal = [0.000000]
		for index, line in enumerate(rmin):
			if index == 0:
				pass
			else:
				line = line.rstrip().split()
				if index == int(line[0].replace(":", "")) - 1:
					new_list = line[1:index+1]
					new_list.append('0')
					for element in new_list:
						lower_diagonal.append(float(element))
		# Build the whole matrix
		n = args.Nsnps - 1
		mat = np.zeros((n,n)) # Initialize nxn matrix 
		# Find lower left indices of a triangular nxn matrix
		tril = np.tril_indices(n)
		# Find upper right indices of a triangular nxn matrix 
		triu = np.triu_indices(n, 1)
		mat[tril] = lower_diagonal
		# Make the matrix symmetric
		mat[triu] = mat.T[triu]
		# ax = sns.heatmap(mat)
		# fig = ax.get_figure()
		# fig.savefig(rmin_png)
		with open(pairwise_rmin_out, "w") as pairwise_rmin_out_file:
			pairwise_rmin_out_file.write("chr"+"\t"+"snp1"+"\t"+"snp2"+"\t"+"\t"+"snp1_pos"+"\t"+"snp2_pos"+"\t"+"rho"+"\n")
			for snppair in zip(triu[0], triu[1]):
				snp1 = "SNP." + str(snppair[0] + 1)
				snp2 = "SNP." + str(snppair[1] + 1)
				snp1_original_pos = geno_position[snppair[0]]
				snp2_original_pos = geno_position[snppair[1]]
				pairwise_rmin_out_file.write(args.chr+"\t"+snp1+"\t"+snp2+"\t"+snp1_original_pos+"\t"+snp2_original_pos+"\t"+str(mat[snppair])+"\n")
				geno_position[start:end]

		# window_out = args.chr + ":" + str(interval_counter) + ".ldhat." + "window_out.txt"
		# window_png = args.chr + ":" + str(interval_counter) + ".window.png"
		# window_txt = args.chr + ":" + str(interval_counter) + ".originalPOS.window.txt"
		# with open(window_out, 'r') as w_file:
		# 	w = pd.read_table(w_file, \
		# 		skip_blank_lines=True, sep ="\t", skipinitialspace=True, \
		# 		skiprows=lambda x: x in [0, 1, 2])

		# 	w_no_na = w.dropna()

		# 	mid_position = (w_no_na["SNP_L"].astype(float)+w_no_na["SNP_R"])/2
		# 	mid = mid_position/1000
		# 	rho = w_no_na["4Ner/bp/kb"]

		# 	plt.figure()
		# 	plt.plot(mid, rho*1000)
		# 	plt.ylabel("4Ner/kb")
		# 	plt.xlabel("Position (kb)")
		# 	plt.savefig(window_png)

		# 	# Create a dictionary with new and original coordinates
		# 	new_coord_to_original = dict(zip(new_coord, geno_position))

		# 	original_snp_l = [new_coord_to_original[str(int(float(item)))] for index, item in enumerate(w_no_na["SNP_L"])]
		# 	original_snp_r = [new_coord_to_original[str(int(float(item)))] for index, item in enumerate(w_no_na["SNP_R"])]

		# 	w_no_na.insert(1, "SNP_L_POS", original_snp_l, True)
		# 	w_no_na.insert(3, "SNP_R_POS", original_snp_r, True)

		# 	# Saving dataframe as CSV
		# 	#with open(window_txt, "w") as df_out:
		# 	w_no_na.to_csv(window_txt, sep="\t", index = False)


# Create animation of all likelihood plots
# fig = plt.figure()
# def animate(interval_counter):
# 	composite_out = "superscaffold36" + ":" + str(interval_counter) + ".ldhat." + "outfile.txt"
# 	f = open(composite_out, "r")
# 	lk = ""
# 	for line in f:
# 		if line.startswith("Maximum"):
# 			line = line.rstrip()
# 			lk = lk + line
# 	outfile = pd.read_table(composite_out, \
# 	skip_blank_lines=True, skipinitialspace=True, sep='\s+',\
# 	skiprows=lambda x: x in [0, 1, 2, 3, 4, 5])
# 	x = outfile['4Ner(region)']
# 	y = outfile['Pairwise']
# 	plt.cla()
# 	im = plt.plot(x, y)
# 	plt.xlabel('4Ner (region)')
# 	plt.ylabel('Composite-likelihood')
# 	plt.title(lk)
# 	return im


# ani = matplotlib.animation.FuncAnimation(fig, animate, range(1, 21), repeat=False, interval = 1000)
# # Set up formatting for the movie files
# Writer = matplotlib.animation.FFMpegWriter(fps=30, codec='libx264')
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("ld_images.mp4")

# plt.show()


