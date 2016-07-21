import argparse
import sys

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', help='input file [csv]')
    parser.add_argument('--out_file', help='output file [fasta]')
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    if not (args.in_file and args.out_file):
        print(parser.usage)
        sys.exit(0)

    file_in = args.in_file
    file_out = args.out_file
    
    with open(file_out, "wt") as fd:
        with open(file_in, "rt") as lines:
            for line in lines:
                if line.find(",") == -1:
                    continue
                seqid, seq, fr1, fr2, fr3, fr4, fr1i, fr2i, fr3i, fr4i, shift, rflag = line.strip().split(',')
                if int(rflag):
                    seq = seq[::-1].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()
                seq = seq[int(shift):]
                fr1_start = int(fr1i.split()[0]) * 3 if fr1i != '.' else 0
                fr4_end   = int(fr4i.split()[1]) * 3 if fr4i != '.' else len(seq)
                fd.write(">%s\n" % seqid)
                fd.write("%s\n" % seq[fr1_start:fr4_end])

if __name__ == "__main__":
    main()
