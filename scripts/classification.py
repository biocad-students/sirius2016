import numpy as np
import os
import os.path
from Bio import Seq
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import time
import argparse

GAP_PENALTY = 5
COINCIDENCE_SCORE = 4

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
	return [(1 - (score/(finish - start + 1))), start, finish]

def f(x, y):
	if x==y:
		return COINCIDENCE_SCORE
	return -COINCIDENCE_SCORE

def sga(s1, s2):
	n, m = int(len(s1)), int(len(s2))
	d = np.zeros((n + 1, m + 1))
	for i in range(1, n + 1):
		for j in range(1, m + 1):
			d[i][j] = max(d[i][j-1] - GAP_PENALTY, d[i-1][j] - GAP_PENALTY, d[i - 1][j - 1] + f(s1[i-1], s2[j-1]))
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
		arr = [d[x][y-1] - GAP_PENALTY, d[x-1][y] - GAP_PENALTY, d[x][y] + f(s1[x-1], s2[y-1])]
		scr = arr.index(max(arr))
		if (scr == 2):
			ans1 = s1[x-1] + ans1
			ans2 = s2[y-1] + ans2
			x, y = x - 1, y - 1
		elif (scr == 0):
			ans2 = s2[y-1] + ans2
			ans1 = "-" + ans1
			x, y = x, y-1
		else:
			ans2 = "-" + ans2
			ans1 = s1[x-1] + ans1
			x, y = x-1, y
	ans1 = s1[:x] + "-" * y + ans1
	ans2 = s2[:y] + "-" * x + ans2
	return similarity(ans1, ans2)

def to_csv_format(a):
	record, table_names, table_scores, table_params, frame, direction = a
	return (",".join([record.id, str(record.seq), ",".join(table_names), ",".join(table_params), str(frame), str(direction)])) + "\n"

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--in_file', help='input file [fasta]')
	parser.add_argument('--out_dir', help='output directory')
	parser.add_argument('--path_germline', help='directory of germlines')
	parser.add_argument('--is_heavy', type=bool, help='1 if file contains HC, 0 if it contains LC')
	return parser.parse_args()


# In[5]:
def main():
	tmp = parse_args()
	path_file = tmp.in_file
	path_out = tmp.out_dir
	path_germline =  tmp.path_germline
	is_vh = tmp.is_heavy



	data = []
	for i in SeqIO.parse(path_file, format = "fasta"):
		if (len(str(i.seq)) >= 300):
			data.append(i)

	s = time.time()
	if (is_vh):
		classes = ["VH"]
	else:
		classes = ["VL", "VK"]

	frames_classes = dict()
	for type_ig in classes:
		frames_classes[type_ig] = [[] for i in range(5)]
		for fr in range(1, 5):
			for j in SeqIO.parse(os.path.join(path_germline, os.path.join(type_ig, "FR" + str(fr) + ".fasta")), format = "fasta"):
				frames_classes[type_ig][fr].append(j)

	good = {x:open(os.path.join(path_out, "good." + x + ".csv"), "at") for x in ["VH", "VL", "VK"]}
	bad = {x:open(os.path.join(path_out, "bad." + x + ".csv"), "at") for x in ["VH", "VL", "VK"]}
	trash = open(os.path.join(path_germline, "trash.csv"), "at")
	for record in data:
		processing_record = record
		found_good = False
		for direction in range(2):
			for frame in range(3):
				if (found_good):
					break
				current_seq = str(Seq(str(processing_record.seq)[frame:]).translate())
				for type_ig in classes:
					current_score = 0 
					table_scores = [0] * 5
					table_names = ["."] * 5
					table_params = ["."] * 5
					for ig_frame in range(1, 5):
						max_score = [0]
						ig_name_max = "" 
						for frame_type in frames_classes[type_ig][ig_frame]:
							score = sga(current_seq, str(frame_type.seq))
							if (score[0] > max_score[0]):
								max_score = score
								ig_name_max = frame_type.id
						if (max_score[0] > 0.75):
							table_scores[ig_frame] = 1
							table_names[ig_frame] = ig_name_max
							table_params[ig_frame] = str(max_score[1]) + " " + str(max_score[2])
						else:
							table_scores[ig_frame] = 0
					if (sum(table_scores) == 4):
						good[type_ig].write(to_csv_format((record, table_names[1:], table_scores[1:], table_params[1:], frame, direction)))
						found_good = True
					elif (sum(table_scores) > 1):
						bad[type_ig].write(to_csv_format((record, table_names[1:], table_scores[1:], table_params[1:], frame, direction)))
						found_good = True			  
			processing_record = processing_record.reverse_complement()
		if (not found_good):
			trash.write(record.id + "\n")

	trash.close()
	for i in good.values():
		i.close()
	for i in bad.values():
		i.close()
	print(time.time() - s)

if __name__ == "__main__":
	main()
