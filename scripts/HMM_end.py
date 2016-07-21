from Bio import SeqIO
from Bio import Seq
import argparse
import os.path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_in_good', help='input file [fasta]')
    parser.add_argument('--file_in_bad', help='input file [fasta]')
    parser.add_argument('--dir_out', help='output directory [fasta]')
    return parser.parse_args()


def viterbi(obs, states, start_p, trans_p, emit_p):
    V = [{}]

    for st in states:  # filling base
        try:
            V[0][st] = {"prob": start_p[st] * emit_p[st][obs[0]], "prev": None}
        except:
            V[0][st] = {"prob": 0, "prev": None}
    for t in range(1, len(obs)):  # filling dinamic
        V.append({})

        for st in states:
            max_tr_prob = 0

            for prev_st in states:
                try:
                    prob = V[t-1][prev_st]["prob"] * trans_p[prev_st][st]
                    max_tr_prob = max(prob, max_tr_prob)
                except:
                    pass

            for prev_st in states:
                try:
                    prob = V[t-1][prev_st]["prob"] * trans_p[prev_st][st]
                except:
                    prob = 0
                if prob == max_tr_prob:
                    try:
                        max_prob = max_tr_prob * emit_p[st][obs[t]]
                    except:
                        max_prob = 0
                    V[t][st] = {"prob": max_prob, "prev": prev_st}
                    break
            opt = []
    max_prob = max(value["prob"] for value in V[-1].values())
    previous = None

    for st, data in V[-1].items():  # regen path
        if data["prob"] == max_prob:
            opt.append(st)
            previous = st
            break

    for t in range(len(V) - 2, -1, -1):  # filling path
        opt.insert(0, V[t + 1][previous]["prev"])
        previous = V[t + 1][previous]["prev"]
    return (V, opt, max_prob)


def hamdist(str1, str2):
    diffs = 0

    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1

    return diffs


def main():
    tmp = parse_args()
    file_in_good = open(tmp.in_file_good, 'r')
    file_in_bad = open(tmp.in_file_bad, 'r')
    bads = []
    n = len(good)
    LEN = 3
    for i in SeqIO.parse(tmp.in_file_bad, format="fasta"):
        bads.append(str(i.seq))
    number = len(bads)
    profile1 = tmp.in_file_good
    in_file = tmp.in_file_bad
    out_file = os.path.join(tmp.dir_out, 'outfile.fasta')
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file,
                                                 profile1=profile1,
                                                 outfile=out_file,
                                                 verbose=True, auto=True)
    clustalomega_cline()
    arr = []
    for i in SeqIO.parse(os.path.join(tmp.dir_out, 'outfile.fasta'),
                         format="fasta"):
        arr.append(str(i.seq))
    good = arr[:len(arr) - number]  # train
    for l in range(number):
        bad = arr[len(arr) - number + i]  # for each bad from chosed bads
        i = 0
        while i != len(bad):  # removing inserts
            flag = True
            for j in range(n):
                if good[j][i] != '-':
                    flag = False
            if bad[i] != '-' and flag:
                bad = bad[:i] + bad[i + 1:]
                for j in range(n):
                    good[j] = good[j][:i] + good[j][i + 1:]
                i -= 1
            i += 1
        emit_p = dict()  # prob of transition from visible to hidden states
        trans_p = dict()  # prob of transition from hidden to hidden states
        start_p = dict()  # prob of starting states
        states = []  # hidden states
        statesset = set()  # done
        obs = []  # visible states
        start_p['0' + bad[:LEN]] = 1
        last = [bad[:LEN]]
        emit_p['0' + bad[:LEN]] = {'0' + bad[:LEN]: 1}
        for k in range(len(bad) - LEN):
            s = bad[k:k + LEN]  # filling obs
            obs.append(str(k) + s)
            newlast = []

            for i in range(len(last)):  # filling trans_p, emit_p and states
                if bad[k + LEN] == '-':
                    arr = []
                    summ = 0
                    s1 = last[i] + bad[k + LEN]
                    s2 = str(k) + last[i]
                    for j in range(n):
                        dist = hamdist(s1, good[j][k:k + LEN + 1])
                        arr.append((good[j][k + LEN], dist))
                        summ += dist
                    trans_p[s2] = dict()
                    states.append(s2)
                    statesset.add(s2)

                    for j in range(len(arr)):
                        s3 = str(k + 1) + last[i][1:] + arr[j][0]
                        s5 = str(k + 1) + bad[k + 1:k + LEN + 1]
                        p = 2 / n - arr[j][1] / summ
                        if s3 not in trans_p[s2]:
                            trans_p[s2][s3] = p
                            emit_p[s3] = dict()
                            emit_p[s3][s5] = p
                            if k + 1 == len(bad) - LEN and s3 not in statesset:
                                states.append(s3)
                                statesset.add(s3)
                            if '-' not in last[i][1:] + arr[j][0]:
                                newlast.append(last[i][1:] + arr[j][0])
                        else:
                            s4 = str(k + 1) + last[i][1:] + arr[j][0]
                            trans_p[s2][s4] += p
                            emit_p[s4][s5] += p
                else:
                    s6 = str(k + 1) + last[i][1:] + bad[k + LEN]
                    trans_p[str(k) + last[i]] = dict()
                    trans_p[str(k) + last[i]][s6] = 1
                    states.append(str(k) + last[i])
                    statesset.add(str(k) + last[i])
                    emit_p[s6] = dict()
                    emit_p[s6][str(k + 1) + bad[k + 1:k + LEN + 1]] = 1
                    if k + 1 == len(bad) - LEN and s6 not in statesset:
                        states.append(s6)
                        statesset.add(s6)
                    if '-' not in last[i][1:] + bad[k + LEN]:
                        newlast.append(last[i][1:] + bad[k + LEN])
            last = newlast

        obs.append(str(k + 1) + bad[k + 1:k + LEN + 1])
        V, opt, max_prob = viterbi(obs, states, start_p, trans_p, emit_p)
        ans = [opt[0][1:]]
        for i in opt[1:]:
            ans.append(i[-1])
        file_out = open(os.path.join(tmp.dir_out, str(l) + '_bad.fasta'), 'a')
        file_out.write(''.join(ans))
        file_out.close()
    file_in_bad.close()
    file_in_good.close()
if __name__ == "__main__":
    main()
