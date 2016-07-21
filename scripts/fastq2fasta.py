from os import listdir
from Bio import SeqIO
import os.path
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_dir', help='input directory')
    parser.add_argument('--out_dir', help='output directory')
    return parser.parse_args()


tmp = parse_args()
in_directory = tmp.in_dir
out_directory = tmp.out_dir
files = listdir(in_directory)
files = filter(lambda x: x.endswith('.fastq'), files)
for file in  files:
    outfile = open(os.path.join(out_directory,file[0:file.rfind(".")]+".fasta"),"w")
    arr = []
    for record in SeqIO.parse(os.path.join(in_directory,file), format="fastq"):
        arr.append(record)
    SeqIO.write(arr, outfile, "fasta")
    outfile.close()
