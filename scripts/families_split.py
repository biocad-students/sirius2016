from Bio.Seq import Seq
import argparse
import os.path
from time import sleep

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', help='input file [csv]')
    parser.add_argument('--out_dir', help='output directory')
    parser.add_argument('--is_good', type=int, help='1 if file contains bad seq, 0 if it doesn`t')
    return parser.parse_args()


def main():
    tmp = parse_args()
    file_in = open(tmp.in_file, 'r')  # входной файл
    number = 0
    flag = bool(tmp.is_good)
    for record in file_in:
        b=0
        number += 1
        score = [0] * 10
        name = []
        l = 0
        h = 0
        split = record.split(',')
        if int(split[-1]):
            dna = Seq(split[1])
            sequnse = dna.reverse_complement()
        else:
            sequnse = split[1]
        sequnse_out = sequnse[int(split[-2]):]
        if flag:
            s = split[6].split(' ')
            s2 = split[9].split(' ')
            c1 = int(s[0]) * 3
            c2 = int(s2[1]) * 3 + 3
            sequnse = sequnse_out[c1:c2]
        for i in split[2:6]:
            if i != '.':
                name.append(i)
                # чуть-чуть велосипеда/индусского кода
        for i in name:
            b+=1
            if int(i[7]) == 3:
                if i[3] == 'l':
                    l += 2
                else:
                    h += 2
                    score[int(i[9])] += 2
            if int(i[7]) == 2:
                if i[3] == 'l':
                    l += 1
                else:
                    score[int(i[9])] += 1
                    h += 1
            if int(i[7]) == 1:
                if i[3] == 'l':
                    l += 1.5
                else:
                    h += 1.5
                    score[int(i[9])] += 1.5
        print(score)
        if l > h:
            file_out = open(os.path.join(tmp.out_dir, 'lama.csv'),'a')
            file_out.write(record)
            file_out.close()
        else:
            file_name = score.index(max(score))
            file_name = 'human_' + str(file_name) + '.csv'
            file_out = open(os.path.join(tmp.out_dir, file_name),'a')
            file_out.write(record)
            file_out.close()

if __name__ == "__main__":
    main()
