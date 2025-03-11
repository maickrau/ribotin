#include <sstream>
#include <map>
#include <fstream>
#include <set>
#include "ClusterMisc.h"

Node::Node(size_t id, bool forward) :
	value(id + (forward ? 0x8000000000000000 : 0))
{
}

size_t Node::id() const
{
	return value & 0x7FFFFFFFFFFFFFFF;
}

bool Node::forward() const
{
	return (value & 0x8000000000000000) == 0x8000000000000000;
}

bool Node::operator<(const Node& other) const
{
	return value < other.value;
}

bool Node::operator==(const Node& other) const
{
	return value == other.value;
}

bool Node::operator!=(const Node& other) const
{
	return value != other.value;
}

size_t Node::rawValue() const
{
	return value;
}

std::vector<char> getComplement()
{
	std::vector<char> result;
	result.resize(256, 0);
	result['a'] = 'T';
	result['A'] = 'T';
	result['c'] = 'G';
	result['C'] = 'G';
	result['g'] = 'C';
	result['G'] = 'C';
	result['t'] = 'A';
	result['T'] = 'A';
	return result;
}

std::vector<char> complement = getComplement();

std::string revcomp(std::string seq)
{
	std::reverse(seq.begin(), seq.end());
	for (size_t i = 0; i < seq.size(); i++)
	{
		seq[i] = complement[seq[i]];
	}
	return seq;
}

Node reverse(Node node)
{
	return Node { node.id(), !node.forward() };
}

std::string reverse(std::string node)
{
	if (node[0] == '>')
	{
		node[0] = '<';
	}
	else
	{
		assert(node[0] == '<');
		node[0] = '>';
	}
	return node;
}

std::vector<Node> reverse(const std::vector<Node>& path)
{
	std::vector<Node> bw { path.rbegin(), path.rend() };
	for (size_t i = 0; i < bw.size(); i++)
	{
		bw[i] = reverse(bw[i]);
	}
	return bw;
}

std::vector<std::string> reverse(const std::vector<std::string>& path)
{
	std::vector<std::string> bw { path.rbegin(), path.rend() };
	for (size_t i = 0; i < bw.size(); i++)
	{
		bw[i] = reverse(bw[i]);
	}
	return bw;
}

std::pair<std::vector<Node>, std::vector<Node>> canon(const std::vector<Node>& left, const std::vector<Node>& right)
{
	if (left < right) return std::make_pair(left, right);
	return std::make_pair(right, left);
}

std::pair<std::vector<std::string>, std::vector<std::string>> canon(const std::vector<std::string>& left, const std::vector<std::string>& right)
{
	if (left < right) return std::make_pair(left, right);
	return std::make_pair(right, left);
}

std::vector<Node> canon(const std::vector<Node>& path)
{
	auto bw = reverse(path);
	if (bw < path) return bw;
	return path;
}

std::vector<std::string> canon(const std::vector<std::string>& path)
{
	auto bw = reverse(path);
	if (bw < path) return bw;
	return path;
}

std::pair<Node, Node> canon(Node left, Node right)
{
	if (reverse(right) < left) return std::make_pair(reverse(right), reverse(left));
	if (left < reverse(right)) return std::make_pair(left, right);
	if (reverse(left) < right) return std::make_pair(reverse(right), reverse(left));
	return std::make_pair(left, right);
}

std::pair<std::string, std::string> canon(const std::string& left, const std::string& right)
{
	auto revleft = reverse(left);
	auto revright = reverse(right);
	if (revright < left) return std::make_pair(revright, revleft);
	if (left < revright) return std::make_pair(left, right);
	if (revleft < right) return std::make_pair(revright, revleft);
	return std::make_pair(left, right);
}

size_t find(std::vector<size_t>& parent, size_t key)
{
	while (parent.at(key) != parent.at(parent.at(key))) parent[key] = parent[parent[key]];
	return parent[key];
}

void merge(std::vector<size_t>& parent, size_t left, size_t right)
{
	left = find(parent, left);
	right = find(parent, right);
	assert(parent.at(left) == left);
	assert(parent.at(right) == right);
	parent[right] = left;
}

void GfaGraph::loadFromFile(std::string gfaFile)
{
	std::ifstream file { gfaFile };
	std::set<std::tuple<Node, Node, size_t, size_t>> fileedges;
	while (file.good())
	{
		std::string line;
		getline(file, line);
		std::stringstream sstr { line };
		std::string test;
		sstr >> test;
		if (test == "S")
		{
			std::string node, seq, coveragetag;
			sstr >> node >> seq >> coveragetag;
			if (nodeNameToId.count(node) == 0)
			{
				size_t id = nodeNameToId.size();
				nodeNameToId[node] = id;
				nodeSeqs.emplace_back();
				revCompNodeSeqs.emplace_back();
				nodeCoverages.emplace_back(0);
				nodeNames.push_back(node);
			}
			nodeSeqs[nodeNameToId[node]] = seq;
			revCompNodeSeqs[nodeNameToId[node]] = revcomp(seq);
			nodeCoverages[nodeNameToId[node]] = std::stof(coveragetag.substr(5));
		}
		else if (test == "L")
		{
			std::string fromnodename;
			std::string tonodename;
			std::string fromorient, toorient;
			std::string overlapstr, edgecoveragetag;
			sstr >> fromnodename >> fromorient >> tonodename >> toorient >> overlapstr >> edgecoveragetag;
			if (nodeNameToId.count(fromnodename) == 0)
			{
				size_t id = nodeNameToId.size();
				nodeNameToId[fromnodename] = id;
				nodeSeqs.emplace_back();
				revCompNodeSeqs.emplace_back();
				nodeCoverages.emplace_back(0);
				nodeNames.push_back(fromnodename);
			}
			if (nodeNameToId.count(tonodename) == 0)
			{
				size_t id = nodeNameToId.size();
				nodeNameToId[tonodename] = id;
				nodeSeqs.emplace_back();
				revCompNodeSeqs.emplace_back();
				nodeCoverages.emplace_back(0);
				nodeNames.push_back(tonodename);
			}
			Node fromnode, tonode;
			if (fromorient == "+")
			{
				fromnode = Node { nodeNameToId[fromnodename], true };
			}
			else
			{
				fromnode = Node { nodeNameToId[fromnodename], false };
			}
			if (toorient == "+")
			{
				tonode = Node { nodeNameToId[tonodename], true };
			}
			else
			{
				tonode = Node { nodeNameToId[tonodename], false };
			}
			overlapstr.pop_back();
			size_t overlap = std::stoull(overlapstr);
			size_t coverage = std::stoull(edgecoveragetag.substr(5));
			fileedges.emplace(fromnode, tonode, overlap, coverage);
			fileedges.emplace(reverse(tonode), reverse(fromnode), overlap, coverage);
		}
	}
	for (auto t : fileedges)
	{
		edges[std::get<0>(t)].emplace(std::get<1>(t), std::get<2>(t), std::get<3>(t));
	}
	for (size_t i = 0; i < nodeSeqs.size(); i++)
	{
		nodeSeqsTwobit.emplace_back(nodeSeqs[i]);
		revCompNodeSeqsTwobit.emplace_back(revCompNodeSeqs[i]);
	}
}

size_t GfaGraph::numNodes() const
{
	return nodeNameToId.size();
}

Path::Path() :
	nodes(),
	overlaps(),
	leftClip(0),
	rightClip(0)
{
}

std::string Path::getSequence(const std::vector<std::string>& nodeSeqs) const
{
	std::string result;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		std::string add;
		add = nodeSeqs.at(nodes[i].id());
		if (!nodes[i].forward())
		{
			add = revcomp(add);
		}
		if (i > 0) add = add.substr(overlaps[i]);
		result += add;
	}
	result = result.substr(leftClip, result.size() - leftClip - rightClip);
	return result;
}

PathSequenceView::PathSequenceView(const std::vector<Node>& nodes, const std::vector<TwobitString>& nodeSeqs, const std::vector<TwobitString>& revCompNodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges, size_t leftClip, size_t rightClip)
{
	size_t realPos = 0;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		if (nodes[i].forward())
		{
			nodePointers.emplace_back(&(nodeSeqs[nodes[i].id()]));
		}
		else
		{
			nodePointers.emplace_back(&(revCompNodeSeqs[nodes[i].id()]));
		}
		size_t overlap = 0;
		if (i > 0)
		{
			overlap = std::numeric_limits<size_t>::max();
			for (auto edge : edges.at(nodes[i-1]))
			{
				if (std::get<0>(edge) != nodes[i]) continue;
				assert(overlap == std::numeric_limits<size_t>::max());
				overlap = std::get<1>(edge);
			}
			assert(overlap != std::numeric_limits<size_t>::max());
		}
		assert(overlap < nodePointers.back()->size());
		nodeStartOverlap.push_back(overlap);
		nodeStartPoses.push_back(realPos);
		realPos += nodePointers.back()->size() - overlap;
	}
	this->leftClip = leftClip;
	this->rightClip = rightClip;
}

char PathSequenceView::operator[](size_t index) const
{
	if (leftClip > 0) index += leftClip;
	size_t nodeIndex = 0;
	while (nodeIndex+1 < nodeStartPoses.size() && nodeStartPoses[nodeIndex+1] <= index) nodeIndex += 1;
	size_t nodeOffset = index - nodeStartPoses[nodeIndex] + nodeStartOverlap[nodeIndex];
	return (*(nodePointers[nodeIndex])).getChar(nodeOffset);
}

size_t PathSequenceView::size() const
{
	return nodeStartPoses.back() + nodePointers.back()->size() - nodeStartOverlap.back() - leftClip - rightClip;
}

size_t getOverlap(const Node& from, const Node& to, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges)
{
	size_t overlap = std::numeric_limits<size_t>::max();
	for (auto edge : edges.at(from))
	{
		if (std::get<0>(edge) != to) continue;
		assert(overlap == std::numeric_limits<size_t>::max());
		overlap = std::get<1>(edge);
	}
	assert(overlap != std::numeric_limits<size_t>::max());
	return overlap;
}

size_t getPathLength(const std::vector<Node>& nodes, const std::vector<std::string>& nodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges)
{
	size_t result = 0;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		size_t add = nodeSeqs.at(nodes[i].id()).size();
		if (i > 0)
		{
			size_t overlap = getOverlap(nodes[i-1], nodes[i], edges);
			add -= overlap;
		}
		result += add;
	}
	return result;
}

PathSequenceView getSequenceView(const std::vector<Node>& nodes, const std::vector<TwobitString>& nodeSeqs, const std::vector<TwobitString>& revCompNodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges, size_t leftClip, size_t rightClip)
{
	return PathSequenceView(nodes, nodeSeqs, revCompNodeSeqs, edges, leftClip, rightClip);
}
