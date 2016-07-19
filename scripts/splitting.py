from Bio import Seq
from Bio import SeqIO
from Bio.Seq import Seq
import os
import os.path
import sys
import argparse

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--in_file', help='input file [fasta]')
	parser.add_argument('--out_dir', help='output directory')
	parser.add_argument('--len_spl', type=int, help='length of files')
	return parser.parse_args()

tmp = parse_args()
path_in = tmp.in_file
path_out = tmp.out_dir
len_spl = tmp.len_spl

data = []
for i in SeqIO.parse(path_in, format = "fasta"):
	if (len(str(i.seq)) >= 300):
		data.append(i)

num_spl = (len(data) + len_spl - 1) // len_spl
for i in range(num_spl):
	arr = data[i * len_spl : min(len(data), (i + 1) * len_spl)]
	SeqIO.write(arr, open(os.path.join(path_out, str(i) + ".fasta"), "w"), "fasta")
