#!/usr/bin/env python

import sys

graph_file = sys.argv[1]
read_paths_file = sys.argv[2]
heavy_path_file = sys.argv[3]
min_coverage = int(sys.argv[4])
# output to stdout

def revnode(n):
	return ("<" if n[0] == ">" else ">") + n[1:]

def revcomp(s):
	comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
	return "".join(comp[c] for c in s[::-1])

def canonpath(path):
	fwstr = "".join(path)
	bwstr = "".join(revnode(n) for n in path[::-1])
	if bwstr < fwstr: return tuple(revnode(n) for n in path[::-1])
	return tuple(path)

def parse_heavy_path(pathstr):
	parts = pathstr.split("_")
	heavy_path = parts[0].replace("w", "").replace("f", ">").replace("b", "<").replace(">", "\t>").replace("<", "\t<").strip().split("\t")
	heavy_path_nodes = set(n[1:] for n in heavy_path)
	heavy_path_edges = set()
	for i in range(1, len(heavy_path)):
		heavy_path_edges.add(canonpath(heavy_path[i-1:i+1]))
	heavy_path_edges.add(canonpath([heavy_path[-1], heavy_path[0]]))
	heavy_path_flipped = False
	heavy_path_rotate = 0
	for part in parts[1:]:
		if part == "revcomp":
			heavy_path_flipped = True
		elif part[0:7] == "rotate-":
			heavy_path_rotate = int(part[7:])
	return (heavy_path, heavy_path_nodes, heavy_path_edges, heavy_path_flipped, heavy_path_rotate)

def get_path_sequence(nodeseqs, edge_overlaps, path):
	result = ""
	for i in range(0, len(path)):
		overlap = 0
		if i > 0: overlap = edge_overlaps[canonpath([path[i-1], path[i]])]
		add_seq = nodeseqs[path[i][1:]]
		if path[i][0] == "<": add_seq = revcomp(add_seq)
		add_seq = add_seq[overlap:]
		result += add_seq
	return result

def get_heavy_path_start_pos(nodeseqs, edge_overlaps, heavy_path, start_index):
	pos = 0
	for i in range(0, start_index):
		pos += len(nodeseqs[heavy_path[i][1:]]) - edge_overlaps[canonpath([heavy_path[i-1], heavy_path[i]])]
	pos -= edge_overlaps[canonpath([heavy_path[start_index-1], heavy_path[start_index]])]
	return pos

def get_bubble_type(nodeseqs, edge_overlaps, heavy_path, bubble_path):
	heavy_path_start = None
	heavy_path_end = None
	heavy_path_fw = None
	for i in range(0, len(heavy_path)):
		if heavy_path[i][1:] == bubble_path[0][1:]:
			heavy_path_start = i
			path_fw = heavy_path[i][0] == bubble_path[0][0]
			assert heavy_path_fw is None or heavy_path_fw == path_fw
			heavy_path_fw = path_fw
		if heavy_path[i][1:] == bubble_path[-1][1:]:
			heavy_path_end = i
			path_fw = heavy_path[i][0] == bubble_path[-1][0]
			assert heavy_path_fw is None or heavy_path_fw == path_fw
			heavy_path_fw = path_fw
	assert heavy_path_fw is not None
	assert heavy_path_start is not None
	assert heavy_path_end is not None
	if not heavy_path_fw:
		bubble_path = [revnode(n) for n in bubble_path[::-1]]
		(heavy_path_start, heavy_path_end) = (heavy_path_end, heavy_path_start)
	if heavy_path_end > heavy_path_start:
		partial_heavy_path = heavy_path[heavy_path_start:heavy_path_end+1]
	else:
		partial_heavy_path = heavy_path[heavy_path_start:] + heavy_path[:heavy_path_end+1]
	ref_pos = get_heavy_path_start_pos(nodeseqs, edge_overlaps, heavy_path, heavy_path_start)
	ref_seq = get_path_sequence(nodeseqs, edge_overlaps, partial_heavy_path)
	bubble_seq = get_path_sequence(nodeseqs, edge_overlaps, bubble_path)
	assert partial_heavy_path[0] == bubble_path[0]
	assert partial_heavy_path[-1] == bubble_path[-1]
	start_match = 0
	while start_match < len(ref_seq) and start_match < len(bubble_seq) and ref_seq[start_match] == bubble_seq[start_match]:
		start_match += 1
	ref_seq = ref_seq[start_match:]
	bubble_seq = bubble_seq[start_match:]
	end_match = 0
	while end_match < len(ref_seq) and end_match < len(bubble_seq) and ref_seq[-end_match-1] == bubble_seq[-end_match-1]:
		end_match += 1
	if end_match > 0:
		ref_seq = ref_seq[:-end_match]
		bubble_seq = bubble_seq[:-end_match]
	assert len(ref_seq) > 0 or len(bubble_seq) > 0
	ref_pos += start_match
	if len(ref_seq) == 1 and len(bubble_seq) == 1: return (ref_pos, "SNP " + ref_seq + "->" + bubble_seq)
	if len(ref_seq) == 0: return (ref_pos, "INSERTION " + str(len(bubble_seq)) + "bp")
	if len(bubble_seq) == 0: return (ref_pos, "DELETION " + str(len(ref_seq)) + "bp")
	return (ref_pos, "COMPLEX ref:" + ref_seq + " vs alt:" + bubble_seq + " " + str(max(len(ref_seq), len(bubble_seq))) + "bp")

nodeseqs = {}
edges = {}
edge_overlaps = {}
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split("\t")
		if parts[0] == "S":
			nodeseqs[parts[1]] = parts[2]
		elif parts[0] == "L":
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			overlap = int(parts[5][:-1])
			if fromnode not in edges: edges[fromnode] = set()
			if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
			edges[fromnode].add(tonode)
			edges[revnode(tonode)].add(revnode(fromnode))
			edge_overlaps[canonpath([fromnode, tonode])] = overlap

heavy_path = []
heavy_path_edges = set()
heavy_path_nodes = set()
heavy_path_flipped = False
heavy_path_rotate = 0
heavy_path_seq = ""
with open(heavy_path_file) as f:
	for l in f:
		if l[0] == ">":
			(heavy_path, heavy_path_nodes, heavy_path_edges, heavy_path_flipped, heavy_path_rotate) = parse_heavy_path(l[1:].strip())
		else:
			heavy_path_end = l.strip()

bubble_support = {}
with open(read_paths_file) as f:
	for l in f:
		parts = l.strip().split("\t")
		path = parts[5].replace(">", "\t>").replace("<", "\t<").strip().split("\t")
		last_match = None
		for i in range(0, len(path)):
			if path[i][1:] in heavy_path_nodes:
				if last_match is None:
					last_match = i
					continue
				key = canonpath(path[last_match:i+1])
				last_match = i
				if key in heavy_path_edges: continue
				if key not in bubble_support: bubble_support[key] = 0
				bubble_support[key] += 1

valid_bubbles = []
for bubble in bubble_support:
	if bubble_support[bubble] < min_coverage: continue
	valid_bubbles.append((bubble, bubble_support[bubble]))
valid_bubbles.sort(key=lambda x: -x[1])

for bubble in valid_bubbles:
	(ref_start, bubble_type) = get_bubble_type(nodeseqs, edge_overlaps, heavy_path, bubble[0])
	print("".join(bubble[0]) + "\t" + str(bubble[1]) + "\t" + str(ref_start) + "\t" + bubble_type)
