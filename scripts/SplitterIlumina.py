import argparse
from Bio import SeqIO
import os.path
from os import mkdir

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in1', help='input file')
    parser.add_argument('--in2', help='input file')
    parser.add_argument('--out_dir', help='output directory')
    parser.add_argument('--n', help='output directory')
    return parser.parse_args()

n = tmp.n
tmp = parse_args()
in1 = tmp.in1
in2 = tmp.in2
out_directory = tmp.out_dir
arr = []
strin = 1
num = 0
filenum = 0
try:
    mkdir(os.path.join(out_directory,"R1"))
except OSError:
    print("file")
try:
    mkdir(os.path.join(out_directory,"R2"))
except OSError:
    print("file")
infile = open(in1,"r")
outfile = open(os.path.join(out_directory,"R1","R1_"+str(filenum)+".fasta"),"w")
while strin:
    if not (num % (4*n)) and (num!=0):
        filenum += 1
        outfile.close()
        outfile = open(os.path.join(out_directory,"R1","R1_"+str(filenum)+".fasta"),"w")
    strin = infile.readline()
    outfile.write(strin)
    num += 1
outfile.close()
infile.close()
strin = 1
num = 0
filenum = 0
infile = open(in2,"r")
outfile = open(os.path.join(out_directory,"R2","R2_"+str(filenum)+".fasta"),"w")
while strin:
    if not (num % (4*n)) and (num!=0):
        filenum += 1
        outfile.close()
        outfile = open(os.path.join(out_directory,"R2","R2_"+str(filenum)+".fasta"),"w")
    strin = infile.readline()
    outfile.write(strin)
    num += 1
outfile.close()
infile.close()
print(filenum)
