#!/usr/bin/env python

import sys

ref_fasta = sys.argv[1]
path_fasta = sys.argv[2]
# result to stdout

k = 51
max_diagonal_distance = 10
min_diagonal_matches = 20

def revcomp(s):
	comp = {"A": "T", "T": "A", "C": "G", "G": "C", "a": "T", "t": "A", "c": "G", "g": "C"}
	return "".join(comp[c] for c in s[::-1])

ref_seq = ""
with open(ref_fasta) as f:
	for l in f:
		if l[0] == ">": continue
		ref_seq += l.strip()

path_seq = ""
path_name = ""
with open(path_fasta) as f:
	for l in f:
		if l[0] == ">":
			path_name = l[1:].strip()
		else:
			path_seq += l.strip()

kmer_poses = {}
for i in range(0, len(ref_seq)-k):
	kmer = ref_seq[i:i+k]
	if kmer in kmer_poses:
		kmer_poses[kmer] = None
	else:
		kmer_poses[kmer] = i;
rc_ref_seq = revcomp(ref_seq)
for i in range(0, len(rc_ref_seq)-k):
	kmer = rc_ref_seq[i:i+k]
	if kmer in kmer_poses:
		kmer_poses[kmer] = None

fw_matches = []
for i in range(0, len(path_seq)-k):
	kmer = path_seq[i:i+k]
	if kmer not in kmer_poses: continue
	pos = kmer_poses[kmer]
	if pos is None: continue
	fw_matches.append((pos, i))

bw_matches = []
rc_path_seq = revcomp(path_seq)
for i in range(0, len(rc_path_seq)-k):
	kmer = rc_path_seq[i:i+k]
	if kmer not in kmer_poses: continue
	pos = kmer_poses[kmer]
	if pos is None: continue
	bw_matches.append((pos, i))

sys.stderr.write("fw matches: " + str(len(fw_matches)) + "\n")
sys.stderr.write("bw matches: " + str(len(bw_matches)) + "\n")
if len(bw_matches) > len(fw_matches):
	path_seq = revcomp(path_seq)
	path_name += "_revcomp"
	matches = bw_matches
else:
	matches = fw_matches

assert len(matches) > min_diagonal_matches

matches.sort(key=lambda x: x[0]-x[1])

min_refpos = len(ref_seq)
rotate = 0
for i in range(0, len(matches)-min_diagonal_matches):
	startdiagonal = matches[i][0] - matches[i][1]
	enddiagonal = matches[i+min_diagonal_matches-1][0] - matches[i+min_diagonal_matches-1][1]
	if enddiagonal - startdiagonal >= max_diagonal_distance: continue
	min_refpos_here = len(ref_seq)
	rotate_here = 0
	for j in range(i, i+min_diagonal_matches):
		if matches[i][0] < min_refpos_here:
			min_refpos_here = matches[i][0]
			rotate_here = matches[i][0] - matches[i][1]
	if min_refpos_here < min_refpos:
		min_refpos = min_refpos_here
		rotate = rotate_here

assert min_refpos < len(ref_seq)

path_name += "_rotate" + str(rotate)
path_seq = path_seq[-rotate:]+path_seq[:-rotate]

print(">" + path_name)
print(path_seq)
