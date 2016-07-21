from Bio import SeqIO
import numpy as np
import os, sys
from Bio.Seq import Seq
import argparse

GAP_PENALTY = 5
COINCIDENCE_SCORE = 4

# вспомогательная функция для sga
def similarity(x, y):
    start = 0
    while (y[start] == "-"):
        start += 1
    finish = len(x) - 1
    while (y[finish] == "-"):
        finish -= 1
    score = 0
    for i in range(start, finish + 1):
        score += (x[i] != y[i])
    return [(1 - (score / (finish - start + 1))), start, finish]

# вспомогательная функция для sga
def f(x, y):
    if x == y:
        return COINCIDENCE_SCORE
    return -COINCIDENCE_SCORE

# вполуглобальное выравнивание
def sga(s1, s2):
    n, m = int(len(s1)), int(len(s2))
    d = np.zeros((n + 1, m + 1))
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d[i][j] = max(d[i][j - 1] - GAP_PENALTY, d[i - 1][j] - GAP_PENALTY,
                          d[i - 1][j - 1] + f(s1[i - 1], s2[j - 1]))
    ans = 0
    x, y = 0, 0
    for i in range(m + 1):
        if (d[n][i] >= d[x][y]):
            x, y = n, i
    for i in range(n + 1):
        if (d[i][m] >= d[x][y]):
            x, y = i, m
    ans1 = s1[x:] + "-" * (m - y)
    ans2 = s2[y:] + "-" * (n - x)
    while (x * y != 0):
        arr = [d[x][y - 1] - GAP_PENALTY, d[x - 1][y] - GAP_PENALTY, d[x][y] + f(s1[x - 1], s2[y - 1])]
        scr = arr.index(max(arr))
        if (scr == 2):
            ans1 = s1[x - 1] + ans1
            ans2 = s2[y - 1] + ans2
            x, y = x - 1, y - 1
        elif (scr == 0):
            ans2 = s2[y - 1] + ans2
            ans1 = "-" + ans1
            x, y = x, y - 1
        else:
            ans2 = "-" + ans2
            ans1 = s1[x - 1] + ans1
            x, y = x - 1, y
    ans1 = s1[:x] + "-" * y + ans1
    ans2 = s2[:y] + "-" * x + ans2
    return similarity(ans1, ans2)

# словарь, для перевода нуклеотидов в аминокислоты
map = {
    "TC-":"S",
    "CT-":"L",
    "CC-":"P",
    "CG-":"R",
    "AC-":"T",
    "GT-":"V",
    "GC-":"A",
    "GG-":"G",
    "TT-":"X",
    "TA-":"X",
    "TG-":"X",
    "CA-":"X",
    "AT-":"X",
    "AA-":"X",
    "AG-":"X",
    "GA-":"X",
}

# перевод нуклеотидов в аминокислоты
def tr(cur_seq):
    if len(cur_seq) % 3 == 1:
        cur_seq = cur_seq[:-1]
    if len(cur_seq) % 3 == 2:
        cur_seq += '-'

    for i in range(0, len(cur_seq), 3):
        count = 0
        for j in cur_seq[i:i + 3]:
            if j == '-':
                count += 1
        if (count == 1 and cur_seq[i + 2] == '-'):
            ans += map[cur_seq[i:i + 3]]
        elif count == 0:
            ans += Seq(cur_seq[i:i + 3]).translate()
        else:
            ans += 'X'
    return ans

# работа с консолью
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_in', help='first input file')  # файл с bad
    parser.add_argument('--dir_out_b', help='first (bad) output directory')  # папка для записи bad
    parser.add_argument('--dir_out_g', help='second (good) output directory')  # папка для записи good
    parser.add_argument('--file_in_fr', help='frame with VXes')  # папка с FR

    return parser.parse_args()

# словарь для хранения FR
map_fr = {}

# создание словаря для хранения FR

def get_map(arg):
    for i in ["VH_FR", "VL_FR", "VK_FR"]:
        for j in range(1, 5):
            map_fr[str(i + str(j) + ".fasta")] = []

    path_fr = arg
    dirs_fr = os.listdir(path_fr)

    for file in dirs_fr: #бежим по папке со всеми FR(VH, VK, VL)
        file_in_fr = os.path.join(path_fr, file)

        path_fr_loc = file_in_fr
        dirs_fr_loc = os.listdir(path_fr_loc)

        for file_loc in dirs_fr_loc: # бежим по папке (VH, VK, VL) там 4 FR
            file_parse = SeqIO.parse(os.path.join(path_fr_loc, file_loc), "fasta")
            st = file + '_' + file_loc
            for i in file_parse:
                map_fr[st].append(i)

def main():
    args = parse_args()

    get_map(args.file_in_fr)

    path = args.file_in
    dirs = os.listdir(path)

    for file in filter(lambda x: x.contain("bad"), dirs):
        file_in = open(os.path.join(path, file), 'r')
        file_out_b = open(os.path.join(args.dir_out_b, file), "at")
        file_out_g = open(os.path.join(args.dir_out_g, file), "at")

        for str_ in file_in:  # перебор всех последовательностей
            cur = str_.split(',')
            seq = cur[1]

            if int(cur[11]) == 1:  # нужно развернуть
                seq = str(Seq(seq).reverse_complement())

            for q in range(6, 10):  # старт и стоп перевести из аминокислотных в нуклеотидные
                if cur[q] != '.':
                    cur[q] = cur[q].split()
                    cur[q][0] = int(cur[q][0])
                    cur[q][1] = int(cur[q][1])
                    cur[q][0] *= 3
                    cur[q][1] *= 3
                    cur[q][1] += 2

            seq = seq[int(cur[10]):]

            count_fr = 0
            for i in range(2, 6):  # сколько найдено FR
                if cur[i] != '.':
                    count_fr += 1

            if count_fr == 3:
                if cur[5] == '.':  # если нет FR4
                    stop_ = int(cur[8][1]) + 1 - int(cur[10])  # где заканчивается FR3
                    cur_seq = seq[:stop_] + '-' + seq[
                                                  stop_:]  # последовательность для перевода в аминокислоты с одним GAP

                    cur_aa = tr(cur_seq)  # перевод в аминокислоты

                    name = 'V' + cur[2][1]  # определение VH/VL/VH

                    file_in_fr = map_fr[name + "_FR4.fasta"]  # все последовательности FR4

                    fl = '0'  # '0' - еще не найден FR4, что-то другое - имя найденного FR4

                    for fr in file_in_fr:  # пробегаем по FR4
                        tmp = sga(cur_aa, fr.seq)  # выравнивание FR4 на последовательность
                        scor = tmp[0]
                        if scor > 0.75:  # выравнивание хорошее
                            for_wr = str_.split(',')  # нужно заменить поля str_: вставить имя FR4, заменить последовательность на новую и убрать флаг на сдвиг и переворот последовательности
                            for_wr[1] = str(cur_seq)
                            for_wr[11] = '0'
                            for_wr[10] = '0'
                            for_wr[5] = fr.name
                            for_wr[9] = str(tmp[1]) + " " + str(tmp[2])
                            file_out_g.write(','.join(for_wr) + "\n")
                            file_out_g.flush()
                            fl = fr.name
                            break

                    if fl == '0':  # не нашли FR4 с одним GAP
                        cur_seq = seq[0:stop_] + '--' + seq[stop_:]  # последовательность для перевода в аминокислоты с двумя GAP

                        cur_aa = tr(cur_seq)  # перевод в аминокислоты

                        for fr in file_in_fr:  # пробегаем по FR4

                            tmp = sga(cur_aa, fr.seq)  # выравнивание FR4 на последовательность
                            scor = tmp[0]
                            if scor > 0.75:  # выравнивание хорошее
                                for_wr = str_.split(
                                    ',')  # нужно заменить поля str_: вставить имя FR4, заменить последовательность на новую и убрать флаг на сдвиг и переворот последовательности
                                for_wr[1] = str(cur_seq)
                                for_wr[11] = '0'
                                for_wr[10] = '0'
                                for_wr[5] = fr.name
                                for_wr[9] = str(tmp[1]) + " " + str(tmp[2])
                                file_out_g.write(','.join(for_wr) + "\n")
                                file_out_g.flush()
                                fl = fr.name
                                break

                        if fl == '0':  # если так и не нашли FR4
                            file_out_b.write(str_)
                            file_out_b.flush()
                elif cur[2] == '.':  # если нет FR1
                    start_ = int(cur[7][0]) - int(cur[10]) - 1  # где начинается FR2
                    cur_seq = seq[2:start_] + '--' + seq[
                                                     start_:]  # последовательность для перевода в аминокислоты с одним GAP

                    cur_aa = tr(cur_seq)  # перевод в аминокислоты

                    name = 'V' + cur[5][1]  # определение VH/VL/VH

                    file_in_fr = map_fr[name + "_FR1.fasta"]  # все последовательности FR1

                    fl = '0'  # '0' - еще не найден FR4, что-то другое - имя найденного FR1

                    for fr in file_in_fr:  # пробегаем по FR1
                        tmp = sga(cur_aa, fr.seq)  # выравнивание FR1 на последовательность
                        scor = tmp[0]
                        if scor > 0.75:  # выравнивание хорошее
                            for_wr = str_.split(
                                ',')  # нужно заменить поля str_: вставить имя FR4, заменить последовательность на новую и убрать флаг на сдвиг и переворот последовательности
                            for_wr[1] = str(cur_seq)
                            for_wr[11] = '0'
                            for_wr[10] = '0'
                            for_wr[2] = fr.name
                            for_wr[6] = str(tmp[1]) + " " + str(tmp[2])
                            file_out_g.write(','.join(for_wr) + "\n")
                            file_out_g.flush()
                            fl = fr.name
                            break

                    if fl == '0':  # так и не нашли FR1 с одним GAP
                        cur_seq = seq[1:start_] + '-' + seq[start_:]  # последовательность для перевода в аминокислоты с двумя GAP

                        cur_aa = tr(cur_seq)  # перевод в аминокислоты

                        for fr in file_in_fr:  # пробегаем по FR1
                            tmp = sga(cur_aa, fr.seq)  # выравнивание FR1 на последовательность
                            scor = tmp[0]
                            if scor > 0.75:  # выравнивание хорошее
                                for_wr = str_.split(',')  # нужно заменить поля str_: вставить имя FR4, заменить последовательность на новую и убрать флаг на сдвиг и переворот последовательности
                                for_wr[1] = str(cur_seq)
                                for_wr[11] = '0'
                                for_wr[10] = '0'
                                for_wr[2] = fr.name
                                for_wr[6] = str(tmp[1]) + " " + str(tmp[2])
                                file_out_g.write(','.join(for_wr) + "\n")
                                file_out_g.flush()
                                fl = fr.name
                                break

                        if fl == '0':  # так и не нашли FR1
                            file_out_b.write(str_)
                            file_out_b.flush()
            else:
                file_out_b.write(str_)
                file_out_b.flush()


if __name__ == "__main__":
    main()
