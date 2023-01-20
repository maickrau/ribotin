#!/usr/bin/env python

import sys

path_name = sys.argv[1]
# gfa from stdin
# sequence to stdout
# no path outputed!

def revnode(n):
	return (">" if n[0] == "<" else "<") + n[1:]

def revcomp(s):
	comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	return "".join(comp[c] for c in s[::-1])

def canon(n1, n2):
	if revnode(n2) + revnode(n1) < n1 + n2: return (revnode(n2), revnode(n1))
	return (n1, n2)

def pathseq(path):
	result = ""
	for i in range(0, len(path)):
		edge = canon(path[i-1], path[i])
		add = ""
		if path[i][0] == ">":
			add = node_seqs[path[i][1:]]
		else:
			add = revcomp(node_seqs[path[i][1:]])
		add = add[edge_overlaps[edge]:]
		result += add
	return result

edges = {}
node_seqs = {}
node_coverages = {}
edge_overlaps = {}
for l in sys.stdin:
	parts = l.strip().split('\t')
	if parts[0] == "S":
		node_seqs[parts[1]] = parts[2]
		node_coverages[parts[1]] = float(parts[3][5:])
	if parts[0] == 'L':
		fromnode = (">" if parts[2] == "+" else "<") + parts[1]
		tonode = (">" if parts[4] == "+" else "<") + parts[3]
		if fromnode not in edges: edges[fromnode] = set()
		edges[fromnode].add(tonode)
		if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
		edges[revnode(tonode)].add(revnode(fromnode))
		overlap = int(parts[5][:-1])
		edge_overlaps[canon(fromnode, tonode)] = overlap

max_node = None
for node in node_seqs:
	assert node in node_coverages
	if max_node is None or node_coverages[node] > node_coverages[max_node]:
		max_node = node

max_path = [">" + max_node]
seen = set()
seen.add(max_node)
while True:
	node = max_path[-1]
	assert node in edges
	max_edge = None
	for edge in edges[node]:
		if max_edge is None or node_coverages[edge[1:]] > node_coverages[max_edge[1:]]:
			max_edge = edge
	if max_edge == max_path[0]: break
	if max_edge[1:] in seen:
		sys.stderr.write(max_edge + "\n")
	assert max_edge[1:] not in seen
	seen.add(max_edge[1:])
	max_path.append(max_edge)

print(">" + "".join(("fw" if node[0] == ">" else "bw") + node[1:] for node in max_path))
print(pathseq(max_path))
