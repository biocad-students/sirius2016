import argparse
import os.path

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', help='input file [csv]')
    parser.add_argument('--out_dir', help='output directory')
    return parser.parse_args()

def main():
    tmp = parse_args()
    file_in = open(tmp.in_file, 'r')
    for record in file_in:
        b=0
        score = [0] * 10
        name = []
        l = 0
        h = 0
        split = record.split(',')
        for i in split[2:6]:
            if i != '.':
                name.append(i)
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
