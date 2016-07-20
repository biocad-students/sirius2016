import sklearn
import numpy
import Bio
from sklearn.cluster import DBSCAN
from sklearn import metrics
import numpy as np
from sklearn.cluster import dbscan
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio import AlignIO, SeqIO
from Bio.Alphabet import generic_rna
from Bio.Alphabet import generic_dna

a = open("good.VH.csv", "r") #сюда вводить имя файла Макса
out = open("out_file.fasta",'w') #а сюда имя фаста файла
b =  a.readlines()
for i in range(len(b)):
    b[i] = b[i].rstrip()
fast = []
print(b)
for i in b:
    i  = i.split(',')
    obrez = i[-2]
    reverse = i[-1]
    posl = i[1]
    for l in range(int(obrez)):
        posl = posl[1:]
    if reverse:
        posl = Seq(posl)
        posl.reverse_complement()
        posl = str(posl)
    codd = Seq(posl[2:],generic_dna)
    posl = codd.translate()
    fr1 = list(i[6].split())[0]
    fr4 = list(i[9].split())[1]
    posl = posl[int(fr1):int(fr4)]
    fast.append(posl)
count = 0
for i in fast:
    out.write('>' + str(count)+'\n')
    out.write(str(fast[count]) + '\n')
    count += 1
a.close()
out.close()