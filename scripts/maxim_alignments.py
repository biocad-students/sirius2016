from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', help='input file [csv]')
    parser.add_argument('--out_file', help='output file [fasta]')
    return parser.parse_args()

tmp = parse_args()
file_in = tmp.in_file
file_out = tmp.out_file

t1 = []
l1 = []
t2 = []
l2 = []
t3 = []
l3 = []
t4 = []
for x in open(file_in, "r").readlines():
    if (x.count(",") == 0):
        continue
    arr = x.strip().split(",")
    if (arr[-1] == "1"):
        arr[1] = str(Seq(arr[1]).reverse_complement())
    m = arr[1]
    a1, a2, a3, a4 = arr[6].split(), arr[7].split(), arr[8].split(), arr[9].split()
    s = [int(a1[0]) * 3, int(a1[1]) * 3 + 3, int(a2[0]) * 3, int(a2[1]) * 3 + 3, int(a3[0]) * 3, int(a3[1]) * 3 + 3, int(a4[0]) * 3, int(a4[1]) * 3 + 3]
    t1.append(SeqRecord(id=arr[0], description=arr[0], seq=Seq(m[s[0]:s[1]])))
    l1.append(SeqRecord(id=arr[0], description=arr[0], seq=Seq(m[s[1]:s[2]])))
    t2.append(SeqRecord(id=arr[0], description=arr[0], seq=Seq(m[s[2]:s[3]])))
    l2.append(SeqRecord(id=arr[0], description=arr[0], seq=Seq(m[s[3]:s[4]])))
    t3.append(SeqRecord(id=arr[0], description=arr[0], seq=Seq(m[s[4]:s[5]])))
    l3.append(SeqRecord(id=arr[0], description=arr[0], seq=Seq(m[s[5]:s[6]])))
    t4.append(SeqRecord(id=arr[0], description=arr[0], seq=Seq(m[s[6]:s[7]])))

SeqIO.write(l1, file_out + "1", format = "fasta")
SeqIO.write(l2, file_out + "2", format = "fasta")
SeqIO.write(l3, file_out + "3", format = "fasta")

clustalomega_cline = ClustalOmegaCommandline(infile=file_out+"1", outfile=file_out+"10", verbose=True, auto=True, force=True)
clustalomega_cline()
clustalomega_cline = ClustalOmegaCommandline(infile=file_out+"2", outfile=file_out+"20", verbose=True, auto=True, force=True)
clustalomega_cline()
clustalomega_cline = ClustalOmegaCommandline(infile=file_out+"3", outfile=file_out+"30", verbose=True, auto=True, force=True)
clustalomega_cline()

l1 = []
l2 = []
l3 = []
for i in SeqIO.parse(file_out+"10", format="fasta"):
    l1.append(i)
for i in SeqIO.parse(file_out+"20", format="fasta"):
    l2.append(i)
for i in SeqIO.parse(file_out+"30", format="fasta"):
    l3.append(i)
out_str = []
for i in range(len(l1)):
    out_str.append(t1[i] + l1[i] + t2[i] + l2[i] + t3[i] + l3[i] + t4[i])

SeqIO.write(out_str, file_out, format="fasta")  
