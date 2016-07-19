from Bio.Seq import reverse_complement
import numpy as np
from multiprocessing import Pool
from sys import argv
import os.path
import argparse

GAP_PENALTY = 5
COINCIDENCE_SCORE = 4


#Функция возврашает кол-во идентичных элементов не считая гэпы
def similarity(x, y):
    score = 0
    for t in range(len(x)):
        if (x[t] == y[t]) and (x[t]!="-") and (y[t]!="-"):
            score +=1
    return score

#я не знаю что это
def f(x, y):
    if x==y:
        return COINCIDENCE_SCORE
    return -COINCIDENCE_SCORE

#полуглобальное выравнивание Максима
def sga(s1, s2, s3, s4):
    n, m = int(len(s1)), int(len(s2))
    d = np.zeros((n + 1, m + 1))
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d[i][j] = max(d[i][j-1] - GAP_PENALTY, d[i-1][j] - GAP_PENALTY, d[i-1][j-1] + f(s1[i-1], s2[j-1]))
    ans = 0
    x, y = 0, 0
    for i in range(m + 1):
        if (d[n][i] >= d[x][y]):
            x, y = n, i
    for i in range(n + 1):
        if (d[i][m] >= d[x][y]):
            x, y = i, m
    ans1 = s1[x:] + "-" * (m - y)
    ans3 = s3[x:] + "-" * (m - y)
    ans2 = s2[y:] + "-" * (n - x)
    ans4 = s4[y:] + "-" * (n - x)
    while (x * y != 0):
        arr = d[x][y-1] - GAP_PENALTY, d[x-1][y] - GAP_PENALTY, d[x][y] + f(s1[x-1], s2[y-1])
        scr = arr.index(max(arr))
        if (scr == 2):
            ans1 = s1[x-1] + ans1
            ans2 = s2[y-1] + ans2
            ans3 = s3[x-1] + ans3
            ans4 = s4[y-1] + ans4
            x, y = x - 1, y - 1
        elif (scr == 0):
            ans2 = s2[y-1] + ans2
            ans1 = "-" + ans1
            ans4 = s4[y - 1] + ans4
            ans3 = "-" + ans3
            x, y = x, y-1
        else:
            ans2 = "-" + ans2
            ans1 = s1[x-1] + ans1
            ans4 = "-" + ans4
            ans3 = s3[x-1] + ans3
            x, y = x-1, y
    ans1 = s1[:x] + "-" * y + ans1
    ans2 = s2[:y] + "-" * x + ans2
    ans3 = s3[:x] + "-" * y + ans3
    ans4 = s4[:y] + "-" * x + ans4
    return ans1, ans2, ans3, ans4 #

#функция на которую передаеться массив в котором массив с массивами последовательностей, номер, выходной файл

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in1', help='first input file')
    parser.add_argument('--in2', help='second input file')
    parser.add_argument('--out_dir', help='output directory')
    return parser.parse_args()



def pairarr(arr):
    print('File_' + str(arr[1]) + "_open")
    outfile = open(os.path.join(arr[2], "SmallPair_" + str(arr[1]) + ".fasta"), "at")
    for mas in arr[0]:
        seq1 = mas[2]
        seq2 = mas[3]
        desc1 = mas[6]
        desc2 = mas[7]
        seq1 = seq1[0:-1]
        seq2 = seq2[0:-1]
        desc1 = desc1[0:-1]
        desc2 = desc2[0:-1]
        seq2 = reverse_complement(seq2)
        desc2 = desc2[::-1]
        seq1, seq2, desc1, desc2 = sga(seq1, seq2, desc1, desc2)
        seq = str()
        desc = str()
        if similarity(seq1, seq2) > 50:
            for n in range(len(seq1)):
                if seq1[n] == seq2[n]:
                    seq += seq1[n]
                    desc += max([desc1[n], desc2[n]])
                elif seq1[n] == "-":
                    seq += seq2[n]
                    desc += desc2[n]
                elif seq2[n] == "-":
                    seq += seq1[n]
                    desc += desc1[n]
                else:
                    if desc1[n] > desc2[n]:
                        seq += seq1[n]
                        desc += desc1[n]
                    else:
                        seq += seq2[n]
                        desc += desc2[n]
            outfile.write(mas[0])
            outfile.write(seq + "\n")
            outfile.write(mas[4])
            outfile.write(desc + "\n")
        del seq
        del desc
    outfile.close()
    print('File_' + str(arr[1]) + "_Close")

#основная функция которая составляет массивы для функции и параллелит функции
if __name__ == '__main__':
    # print("infile No1:")
    # in1 = input()
    # print("infile No2:")
    # in2 = input()
    # print("OutFolder:")
    # out = input()


    tmp = parse_args()
    in1 = tmp.in1
    in2 = tmp.in2
    out = tmp.out_dir


    infile1 = open(in1, "r")
    infile2 = open(in2, "r")
    strin = 1
    nu = 0
    num = 0
    a = []
    arr = []
    array = []
    while strin:
        if not(num % 1100) and (num != 0):
            arr.append([array, nu, out])
            array = []
            nu +=1
        del a
        a = []
        for i in range(4):
            strin = infile1.readline()
            a.append(strin)
            strin = infile2.readline()
            a.append(strin)
        array.append(a)
        num += 1
    infile1.close()
    infile2.close()
    arr.append([array, nu, out])
    p1 = Pool()
    p1.map(pairarr, arr)
    p1.close()
