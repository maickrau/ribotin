#include <limits>
#include <queue>
#include <map>
#include <cstdlib>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <tuple>
#include <cassert>
#include <set>
#include <vector>
#include <algorithm>
#include <sstream>
#include <phmap.h>
#include <thread>
#include <chrono>
#include <filesystem>
#include "ReadExtractor.h"
#include "fastqloader.h"
#include "ClusterHandler.h"
#include "RibotinUtils.h"
#include "WfaHelper.h"
#include "TwobitString.h"
#include "edlib.h"

const size_t minPhaseCoverage = 10;

class Node
{
public:
	Node() = default;
	Node(const Node& other) = default;
	Node(Node&& other) = default;
	Node& operator=(const Node& other) = default;
	Node& operator=(Node&& other) = default;
	Node(size_t id, bool forward) :
		value(id + (forward ? 0x8000000000000000 : 0))
	{
	}
	size_t id() const
	{
		return value & 0x7FFFFFFFFFFFFFFF;
	}
	bool forward() const
	{
		return (value & 0x8000000000000000) == 0x8000000000000000;
	}
	bool operator<(const Node& other) const
	{
		return value < other.value;
	}
	bool operator==(const Node& other) const
	{
		return value == other.value;
	}
	bool operator!=(const Node& other) const
	{
		return value != other.value;
	}
	size_t rawValue() const
	{
		return value;
	}
private:
	size_t value;
};

namespace std
{
	template <> struct hash<std::pair<std::string, std::string>>
	{
		size_t operator()(const std::pair<std::string, std::string>& x) const
		{
			return std::hash<std::string>{}(x.first) ^ std::hash<std::string>{}(x.second);
		}
	};
	template <> struct hash<Node>
	{
		size_t operator()(Node x) const
		{
			return std::hash<size_t>{}(x.rawValue());
		}
	};
	template <> struct hash<std::pair<Node, Node>>
	{
		size_t operator()(const std::pair<Node, Node>& x) const
		{
			return std::hash<Node>{}(x.first) ^ std::hash<Node>{}(x.second);
		}
	};
	template <> struct hash<std::vector<Node>>
	{
		size_t operator()(const std::vector<Node>& x) const
		{
			size_t result = 0;
			for (size_t i = 0; i < x.size(); i++)
			{
				result *= 3;
				result += std::hash<Node>{}(x[i]);
			}
			return result;
		}
	};
	template <> struct hash<std::pair<std::vector<Node>, std::vector<Node>>>
	{
		size_t operator()(const std::pair<std::vector<Node>, std::vector<Node>>& x) const
		{
			return std::hash<std::vector<Node>>{}(x.first) ^ std::hash<std::vector<Node>>{}(x.second);
		}
	};
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

class OntLoop
{
public:
	std::vector<Node> path;
	std::string readName;
	size_t approxStart;
	size_t approxEnd;
	size_t originalReadLength;
};

class MorphConsensus
{
public:
	std::vector<Node> path;
	std::string sequence;
	std::string name;
	std::vector<OntLoop> ontLoops;
	size_t coverage;
};

class ReadPath
{
public:
	std::string readName;
	size_t readLength;
	size_t readStart;
	size_t readEnd;
	std::vector<Node> path;
	size_t pathStartClip;
	size_t pathEndClip;
	bool reverse;
};

class Variant
{
public:
	Variant() = default;
	std::string name;
	std::vector<Node> path;
	std::vector<Node> referencePath;
	size_t coverage;
	size_t referenceCoverage;
	size_t referenceStartPos;
	size_t referenceEndPos;
	std::string variantSeq;
	std::string referenceSeq;
};

class GfaGraph
{
public:
	void loadFromFile(std::string gfaFile)
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
	size_t numNodes() const
	{
		return nodeNameToId.size();
	}
	std::unordered_map<std::string, size_t> nodeNameToId;
	std::vector<std::string> nodeNames;
	std::vector<size_t> nodeCoverages;
	std::vector<std::string> nodeSeqs;
	std::vector<std::string> revCompNodeSeqs;
	std::vector<TwobitString> nodeSeqsTwobit;
	std::vector<TwobitString> revCompNodeSeqsTwobit;
	std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>> edges;
};

class Path
{
public:
	Path() :
		nodes(),
		overlaps(),
		leftClip(0),
		rightClip(0)
	{
	}
	std::vector<Node> nodes;
	std::vector<size_t> overlaps;
	size_t leftClip;
	size_t rightClip;
	std::string getSequence(const std::vector<std::string>& nodeSeqs) const
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
};

class PathSequenceView
{
public:
	PathSequenceView(const std::vector<Node>& nodes, const std::vector<TwobitString>& nodeSeqs, const std::vector<TwobitString>& revCompNodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges, size_t leftClip, size_t rightClip)
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
	char operator[](size_t index) const
	{
		if (leftClip > 0) index += leftClip;
		size_t nodeIndex = 0;
		while (nodeIndex+1 < nodeStartPoses.size() && nodeStartPoses[nodeIndex+1] <= index) nodeIndex += 1;
		size_t nodeOffset = index - nodeStartPoses[nodeIndex] + nodeStartOverlap[nodeIndex];
		return (*(nodePointers[nodeIndex])).getChar(nodeOffset);
	}
	size_t size() const
	{
		return nodeStartPoses.back() + nodePointers.back()->size() - nodeStartOverlap.back() - leftClip - rightClip;
	}
private:
	std::vector<const TwobitString*> nodePointers;
	std::vector<size_t> nodeStartPoses;
	std::vector<size_t> nodeStartOverlap;
	size_t leftClip;
	size_t rightClip;
	friend size_t getMatchLength(const PathSequenceView& left, const PathSequenceView& right, const size_t leftStart, const size_t rightStart);
};

__attribute__((always_inline)) size_t getBitvectorMatchLength(uint64_t leftFirst, uint64_t leftSecond, uint64_t rightFirst, uint64_t rightSecond)
{
	uint64_t mismatch = (leftFirst ^ rightFirst) | (leftSecond ^ rightSecond);
	uint64_t prefixZerosPlusOneOne = (mismatch-1) ^ mismatch;
	uint64_t prefixMatchPlusOne = popcount(prefixZerosPlusOneOne);
	assert(prefixMatchPlusOne >= 1);
	return prefixMatchPlusOne-1;
}

__attribute__((always_inline)) size_t getMatchLength(const TwobitString& left, const TwobitString& right, const size_t leftStart, const size_t rightStart)
{
	size_t result = 0;
	while (true)
	{
		auto leftBits = left.getBitvectorsPossiblyTruncated(leftStart+result);
		auto rightBits = right.getBitvectorsPossiblyTruncated(rightStart+result);
		auto matchHere = getBitvectorMatchLength(leftBits.first, leftBits.second, rightBits.first, rightBits.second);
		size_t maxMatch = std::min(64 - ((leftStart+result) % 64), 64 - ((rightStart+result) % 64));
		matchHere = std::min(matchHere, maxMatch);
		if (matchHere == 0) return result;
		result += matchHere;
		if (leftStart+result >= left.size() || rightStart+result >= right.size())
		{
			result = std::min(result, left.size()-leftStart);
			result = std::min(result, right.size()-rightStart);
			return result;
		}
	}
}

size_t getMatchLength(const PathSequenceView& left, const PathSequenceView& right, const size_t leftStart, const size_t rightStart)
{
	size_t result = 0;
	size_t leftIndex = leftStart + left.leftClip;
	size_t rightIndex = rightStart + right.leftClip;
	size_t leftSize = left.size();
	size_t rightSize = right.size();
	if (leftIndex == leftSize + left.leftClip) return 0;
	if (rightIndex == rightSize + right.leftClip) return 0;
	assert(leftIndex < leftSize + left.leftClip);
	assert(rightIndex < rightSize + right.leftClip);
	size_t leftStartNode = 0;
	size_t rightStartNode = 0;
	while (true)
	{
		while (leftStartNode+1 < left.nodeStartPoses.size() && left.nodeStartPoses[leftStartNode+1] <= leftIndex) leftStartNode += 1;
		while (rightStartNode+1 < right.nodeStartPoses.size() && right.nodeStartPoses[rightStartNode+1] <= rightIndex) rightStartNode += 1;
		size_t leftOffset = leftIndex - left.nodeStartPoses[leftStartNode] + left.nodeStartOverlap[leftStartNode];
		size_t rightOffset = rightIndex - right.nodeStartPoses[rightStartNode] + right.nodeStartOverlap[rightStartNode];
		auto matchHere = getMatchLength(*left.nodePointers[leftStartNode], *right.nodePointers[rightStartNode], leftOffset, rightOffset);
		if (matchHere == 0) return result;
		result += matchHere;
		if (result >= leftSize - leftStart || result >= rightSize - rightStart)
		{
			result = std::min(result, leftSize - leftStart);
			result = std::min(result, rightSize - rightStart);
			return result;
		}
		leftIndex += matchHere;
		rightIndex += matchHere;
	}
}

PathSequenceView getSequenceView(const std::vector<Node>& nodes, const std::vector<TwobitString>& nodeSeqs, const std::vector<TwobitString>& revCompNodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges, size_t leftClip, size_t rightClip)
{
	return PathSequenceView(nodes, nodeSeqs, revCompNodeSeqs, edges, leftClip, rightClip);
}

std::string getSequence(const std::vector<Node>& nodes, const std::vector<std::string>& nodeSeqs, const std::vector<std::string>& revCompNodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges)
{
	std::string result;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		std::string add;
		if (nodes[i].forward())
		{
			add = nodeSeqs.at(nodes[i].id());
		}
		else
		{
			add = revCompNodeSeqs.at(nodes[i].id());
		}
		if (i > 0)
		{
			size_t overlap = std::numeric_limits<size_t>::max();
			for (auto edge : edges.at(nodes[i-1]))
			{
				if (std::get<0>(edge) != nodes[i]) continue;
				assert(overlap == std::numeric_limits<size_t>::max());
				overlap = std::get<1>(edge);
			}
			assert(overlap != std::numeric_limits<size_t>::max());
			add = add.substr(overlap);
		}
		result += add;
	}
	return result;
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

std::string getPathGaf(const Path& path, const GfaGraph& graph)
{
	std::string pathseq = path.getSequence(graph.nodeSeqs);
	std::string result = "heavy_path\t";
	result += std::to_string(pathseq.size());
	result += "\t0\t";
	result += std::to_string(pathseq.size());
	result += "\t+\t";
	for (const auto& node : path.nodes)
	{
		result += (node.forward() ? ">" : "<");
		result += graph.nodeNames.at(node.id());
	}
	result += "\t";
	result += std::to_string(pathseq.size() + path.leftClip + path.rightClip);
	result += "\t";
	result += std::to_string(path.leftClip);
	result += "\t";
	result += std::to_string(pathseq.size() + path.leftClip);
	result += "\t";
	result += std::to_string(pathseq.size());
	result += "\t";
	result += std::to_string(pathseq.size());
	result += "\t60";
	return result;
}

class CoverageComparer
{
public:
	bool operator()(const std::tuple<Node, Node, double>& lhs, const std::tuple<Node, Node, double>& rhs) const
	{
		return std::get<2>(lhs) < std::get<2>(rhs);
	}
};

std::vector<Node> parsePath(const std::string& pathstr, const std::unordered_map<std::string, size_t>& nodeNameToId)
{
	std::vector<Node> result;
	size_t lastStart = 0;
	for (size_t i = 1; i < pathstr.size(); i++)
	{
		if (pathstr[i] != '<' && pathstr[i] != '>') continue;
		result.emplace_back(nodeNameToId.at(pathstr.substr(lastStart+1, i - lastStart-1)), pathstr[lastStart] == '>');
		lastStart = i;
	}
	result.emplace_back(nodeNameToId.at(pathstr.substr(lastStart+1)), pathstr[lastStart] == '>');
	return result;
}

Path readHeavyPath(const GfaGraph& graph, const std::string& heavyPathGafFile)
{
	std::ifstream file { heavyPathGafFile };
	std::string line;
	getline(file, line);
	auto parts = split(line, '\t');
	Path result;
	size_t pathLength = std::stoull(parts[6]);
	size_t pathStart = std::stoull(parts[7]);
	size_t pathEnd = std::stoull(parts[8]);
	result.leftClip = pathStart;
	result.rightClip = pathLength - pathEnd;
	result.nodes = parsePath(parts[5], graph.nodeNameToId);
	for (size_t i = 0; i < result.nodes.size(); i++)
	{
		Node prev;
		if (i > 0)
		{
			prev = result.nodes[i-1];
		}
		else
		{
			prev = result.nodes.back();
			if (prev == result.nodes[0])
			{
				assert(result.nodes.size() >= 3);
				prev = result.nodes[result.nodes.size()-2];
			}
		}
		assert(graph.edges.count(prev) == 1);
		for (auto edge : graph.edges.at(prev))
		{
			if (std::get<0>(edge) != result.nodes[i]) continue;
			result.overlaps.push_back(std::get<1>(edge));
		}
	}
	assert(result.getSequence(graph.nodeSeqs).size() == pathEnd - pathStart);
	return result;
}

void filterOut(std::unordered_set<size_t>& nodes, const std::unordered_set<size_t>& removeThese)
{
	for (auto node : removeThese)
	{
		if (nodes.count(node) == 1) nodes.erase(node);
	}
}

void filterOut(std::unordered_map<Node, std::unordered_set<Node>>& edges, const std::unordered_set<size_t>& removeThese)
{
	for (auto node : removeThese)
	{
		if (edges.count(Node { node, true }) == 1) edges.erase(Node { node, true });
		if (edges.count(Node { node, false }) == 1) edges.erase(Node { node, false });
	}
	for (auto& pair : edges)
	{
		for (auto node : removeThese)
		{
			if (pair.second.count(Node { node, true }) == 1) pair.second.erase(Node { node, true });
			if (pair.second.count(Node { node, false }) == 1) pair.second.erase(Node { node, false });
		}
	}
}

bool isCircular(const std::unordered_set<size_t>& nodes, const std::unordered_map<Node, std::unordered_set<Node>>& edges, Node startNode)
{
	assert(nodes.count(startNode.id()) == 1);
	std::unordered_set<Node> reachable;
	std::vector<Node> stack;
	stack.push_back(startNode);
	while (stack.size() > 0)
	{
		auto top = stack.back();
		assert(nodes.count(top.id()) == 1);
		stack.pop_back();
		if (reachable.count(top) == 1) continue;
		reachable.insert(top);
		if (edges.count(top) == 0) continue;
		for (auto edge : edges.at(top))
		{
			if (nodes.count(edge.id()) == 0) continue;
			stack.push_back(edge);
			if (edge == startNode) return true;
		}
	}
	return false;
}

void filterOutNonCircularParts(std::unordered_set<size_t>& nodes, std::unordered_map<Node, std::unordered_set<Node>>& edges, size_t startNode)
{
	assert(nodes.count(startNode) == 1);
	std::unordered_set<Node> reachable;
	std::vector<Node> stack;
	stack.push_back(Node { startNode, true });
	stack.push_back(Node { startNode, false });
	while (stack.size() > 0)
	{
		auto top = stack.back();
		stack.pop_back();
		if (reachable.count(top) == 1) continue;
		reachable.insert(top);
		if (edges.count(top) == 0) continue;
		for (auto edge : edges.at(top))
		{
			if (nodes.count(edge.id()) == 0) continue;
			stack.push_back(edge);
		}
	}
	std::unordered_set<size_t> notCircular;
	for (auto node : nodes)
	{
		if (reachable.count(Node { node, true }) == 1 && reachable.count(Node { node, false }) == 1)
		{
			continue;
		}
		notCircular.insert(node);
	}
	assert(notCircular.count(startNode) == 0);
	filterOut(nodes, notCircular);
	filterOut(edges, notCircular);
}

class PathRecWideCoverage
{
public:
	PathRecWideCoverage() = default;
	PathRecWideCoverage(const size_t coverage, const size_t length)
	{
		coverages.emplace_back(coverage, length);
	}
	bool operator<(const PathRecWideCoverage& other) const
	{
		for (size_t i = 0; i < coverages.size() && i < other.coverages.size(); i++)
		{
			if (coverages[i].first < other.coverages[i].first) return true;
			if (coverages[i].first > other.coverages[i].first) return false;
			if (coverages[i].second > other.coverages[i].second) return true;
			if (coverages[i].second < other.coverages[i].second) return false;
		}
		if (coverages.size() < other.coverages.size()) return true;
		return false;
	}
	std::vector<std::pair<size_t, size_t>> coverages;
private:
};

PathRecWideCoverage operator+(const PathRecWideCoverage& left, const PathRecWideCoverage& right)
{
	PathRecWideCoverage result;
	size_t leftIndex = 0;
	size_t rightIndex = 0;
	while (leftIndex < left.coverages.size() && rightIndex < right.coverages.size())
	{
		if (left.coverages[leftIndex].first < right.coverages[rightIndex].first)
		{
			result.coverages.emplace_back(left.coverages[leftIndex]);
			leftIndex += 1;
		}
		else if (left.coverages[leftIndex].first > right.coverages[rightIndex].first)
		{
			result.coverages.emplace_back(right.coverages[rightIndex]);
			rightIndex += 1;
		}
		else
		{
			assert(left.coverages[leftIndex].first == right.coverages[rightIndex].first);
			result.coverages.emplace_back(left.coverages[leftIndex]);
			result.coverages.back().second += right.coverages[rightIndex].second;
			leftIndex += 1;
			rightIndex += 1;
		}
	}
	while (leftIndex < left.coverages.size())
	{
		result.coverages.emplace_back(left.coverages[leftIndex]);
		leftIndex += 1;
	}
	while (rightIndex < right.coverages.size())
	{
		result.coverages.emplace_back(right.coverages[rightIndex]);
		rightIndex += 1;
	}
	return result;
}

Path getRecWidestPath(const std::unordered_set<size_t>& coveredNodes, const std::unordered_map<Node, std::unordered_set<Node>>& coveredEdges, const GfaGraph& graph, const Node start)
{
	assert(coveredNodes.count(start.id()) == 1);
	std::unordered_map<Node, PathRecWideCoverage> nodeRecWideCoverage;
	std::unordered_map<Node, Node> predecessor;
	nodeRecWideCoverage[start] = PathRecWideCoverage { graph.nodeCoverages[start.id()], graph.nodeSeqs[start.id()].size() };
	std::vector<Node> checkStack;
	assert(coveredEdges.count(start) == 1);
	for (auto edge : coveredEdges.at(start))
	{
		if (coveredNodes.count(edge.id()) == 0) continue;
		checkStack.push_back(edge);
	}
	while (checkStack.size() > 0)
	{
		auto top = checkStack.back();
		assert(coveredNodes.count(top.id()) == 1);
		checkStack.pop_back();
		if (predecessor.count(top) == 1)
		{
			assert(nodeRecWideCoverage.count(top) == 1);
			continue;
		}
		assert(predecessor.count(top) == 0);
		assert(nodeRecWideCoverage.count(top) == 0 || top == start);
		bool hasAllNeighbors = true;
		Node bestPredecessor;
		bool hasAny = false;
		assert(coveredEdges.count(reverse(top)) == 1);
		for (auto revedge : coveredEdges.at(reverse(top)))
		{
			if (coveredNodes.count(revedge.id()) == 0) continue;
			auto pre = reverse(revedge);
			if (nodeRecWideCoverage.count(pre) == 0)
			{
				hasAllNeighbors = false;
			}
			else
			{
				if (!hasAny)
				{
					bestPredecessor = pre;
					hasAny = true;
				}
				else if (nodeRecWideCoverage.at(bestPredecessor) < nodeRecWideCoverage.at(pre))
				{
					bestPredecessor = pre;
					hasAny = true;
				}
			}
		}
		assert(hasAny);
		if (!hasAllNeighbors) continue;
		if (hasAllNeighbors)
		{
			predecessor[top] = bestPredecessor;
			nodeRecWideCoverage[top] = nodeRecWideCoverage.at(bestPredecessor) + PathRecWideCoverage { graph.nodeCoverages[top.id()], graph.nodeSeqs[top.id()].size() };
			assert(coveredEdges.count(top) == 1);
			for (auto edge : coveredEdges.at(top))
			{
				if (coveredNodes.count(edge.id()) == 0) continue;
				assert(coveredEdges.count(reverse(edge)) == 1);
				assert(coveredEdges.at(reverse(edge)).count(reverse(top)) == 1);
				checkStack.push_back(edge);
			}
		}
	}
	assert(predecessor.count(start) == 1);
	assert(predecessor.count(reverse(start)) == 0);
	Path result;
	Node pos = predecessor.at(start);
	while (pos != start)
	{
		result.nodes.push_back(pos);
		pos = predecessor.at(pos);
	}
	result.nodes.push_back(start);
	std::reverse(result.nodes.begin(), result.nodes.end());
	return result;
}

std::vector<Node> getCoreNodeOrder(const std::unordered_set<size_t>& nodes, const std::unordered_map<Node, std::unordered_set<Node>>& edges, const std::unordered_set<size_t> coreNodes, const size_t startNode)
{
	std::vector<Node> result;
	result.emplace_back(startNode, true);
	while (true)
	{
		std::vector<Node> stack { result.back() };
		Node pos = result.back();
		Node nextPos = pos;
		std::unordered_set<size_t> visited;
		while (stack.size() >= 1)
		{
			auto top = stack.back();
			stack.pop_back();
			if (visited.count(top.id()) == 1) continue;
			visited.insert(top.id());
			assert(nodes.count(top.id()) == 1);
			for (auto edge : edges.at(top))
			{
				if (nodes.count(edge.id()) == 0) continue;
				if (coreNodes.count(edge.id()) == 1)
				{
					nextPos = edge;
					stack.clear();
					break;
				}
				stack.push_back(edge);
			}
		}
		assert(nextPos != pos);
		result.push_back(nextPos);
		pos = nextPos;
		if (pos.id() == startNode) break;
		assert(result.size() <= coreNodes.size());
	}
	assert(result.size() == coreNodes.size()+1);
	assert(result[0].id() == startNode);
	assert(result.back().id() == startNode);
	return result;
}

bool isReachable(std::unordered_set<size_t>& nodes, std::unordered_map<Node, std::unordered_set<Node>>& edges, const Node start, const Node end)
{
	std::unordered_set<size_t> checked;
	std::vector<Node> stack;
	stack.push_back(start);
	while (stack.size() >= 1)
	{
		auto top = stack.back();
		stack.pop_back();
		if (checked.count(top.id()) == 1) continue;
		checked.insert(top.id());
		if (edges.count(top) == 0) continue;
		for (auto edge : edges.at(top))
		{
			if (nodes.count(edge.id()) == 0) continue;
			stack.push_back(edge);
			if (edge == end) return true;
		}
	}
	return false;
}

void filterOutNonThroughGoers(std::unordered_set<size_t>& nodes, const std::unordered_set<size_t>& removables, std::unordered_map<Node, std::unordered_set<Node>>& edges, const Node start, const Node end)
{
	std::unordered_set<Node> visited;
	std::vector<Node> checkStack { start, reverse(end) };
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		assert(nodes.count(top.id()) == 1);
		if (visited.count(top) == 1) continue;
		visited.insert(top);
		if (top == end) continue;
		if (top == reverse(start)) continue;
		if (edges.count(top) == 0) continue;
		for (auto edge : edges.at(top))
		{
			if (nodes.count(edge.id()) == 0) continue;
			checkStack.push_back(edge);
		}
	}
	std::unordered_set<size_t> removeThese;
	for (size_t node : removables)
	{
		if (visited.count(Node { node, true }) == 1 && visited.count(Node { node, false }) == 1) continue;
		removeThese.insert(node);
	}
	filterOut(nodes, removeThese);
	filterOut(edges, removeThese);
	assert(isReachable(nodes, edges, start, end));
}

void filterOutDefinitelyNonConsensusNodes(std::unordered_set<size_t>& nodes, std::unordered_map<Node, std::unordered_set<Node>>& edges, const Node start, const Node end, const GfaGraph& graph)
{
	if (edges.at(start).count(end) == 1) return; // corner case, just don't handle it and hope for no problem
	std::unordered_set<size_t> nodesInSubgraph;
	std::vector<Node> checkStack { start };
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (nodesInSubgraph.count(top.id()) == 1) continue;
		nodesInSubgraph.insert(top.id());
		if (top == end) continue;
		for (auto edge : edges.at(top))
		{
			if (nodes.count(edge.id()) == 0) continue;
			checkStack.push_back(edge);
		}
	}
	assert(nodesInSubgraph.size() >= 3);
	assert(nodesInSubgraph.count(start.id()) == 1);
	assert(nodesInSubgraph.count(end.id()) == 1);
	nodesInSubgraph.erase(start.id());
	nodesInSubgraph.erase(end.id());
	std::vector<std::pair<size_t, size_t>> nodeAndCoverage;
	std::vector<std::tuple<Node, Node, size_t>> edgeAndCoverage;
	for (size_t node : nodesInSubgraph)
	{
		nodeAndCoverage.emplace_back(node, graph.nodeCoverages[node]);
		Node fw { node, true };
		if (graph.edges.count(fw) == 1)
		{
			for (const auto& target : graph.edges.at(fw))
			{
				edgeAndCoverage.emplace_back(fw, std::get<0>(target), std::get<2>(target));
			}
		}
		Node bw { node, false };
		if (graph.edges.count(bw) == 1)
		{
			for (const auto& target : graph.edges.at(bw))
			{
				edgeAndCoverage.emplace_back(bw, std::get<0>(target), std::get<2>(target));
			}
		}
	}
	filterOut(nodes, nodesInSubgraph);
	filterOut(edges, nodesInSubgraph);
	std::sort(nodeAndCoverage.begin(), nodeAndCoverage.end(), [](auto left, auto right) { return left.second > right.second; });
	std::sort(edgeAndCoverage.begin(), edgeAndCoverage.end(), [](auto left, auto right) { return std::get<2>(left) > std::get<2>(right); });
	size_t minCoverage = nodeAndCoverage[0].second;
	size_t nodeIndex = 0;
	size_t edgeIndex = 0;
	while (minCoverage > 0)
	{
		size_t nextCoverage = 0;
		while (nodeIndex < nodeAndCoverage.size() && nodeAndCoverage[nodeIndex].second >= minCoverage)
		{
			nodes.emplace(nodeAndCoverage[nodeIndex].first);
			nodeIndex += 1;
		}
		while (edgeIndex < edgeAndCoverage.size() && std::get<2>(edgeAndCoverage[edgeIndex]) >= minCoverage)
		{
			edges[std::get<0>(edgeAndCoverage[edgeIndex])].emplace(std::get<1>(edgeAndCoverage[edgeIndex]));
			edges[reverse(std::get<1>(edgeAndCoverage[edgeIndex]))].emplace(reverse(std::get<0>(edgeAndCoverage[edgeIndex])));
			edgeIndex += 1;
		}
		if (nodeIndex < nodeAndCoverage.size()) nextCoverage = nodeAndCoverage[nodeIndex].second;
		if (edgeIndex < edgeAndCoverage.size()) nextCoverage = std::max(nextCoverage, std::get<2>(edgeAndCoverage[edgeIndex]));
		assert(nextCoverage < minCoverage);
		minCoverage = nextCoverage;
		if (isReachable(nodes, edges, start, end)) break;
	}
	filterOutNonThroughGoers(nodes, nodesInSubgraph, edges, start, end);
}

struct DistantNode
{
public:
	DistantNode() = default;
	DistantNode(size_t distance, Node node) :
		distance(distance),
		node(node)
	{
	}
	size_t distance;
	Node node;
};

struct DistantNodeComparator
{
public:
	bool operator()(const DistantNode left, const DistantNode right)
	{
		return left.distance > right.distance;
	}
};

size_t splitUndirectedDistance(const GfaGraph& graph, const std::unordered_set<size_t>& nodes, const std::unordered_map<Node, std::unordered_set<Node>>& edges, const size_t middleNode)
{
	std::unordered_map<Node, size_t> distance;
	std::priority_queue<DistantNode, std::vector<DistantNode>, DistantNodeComparator> queue;
	queue.emplace(0, Node { middleNode, true });
	while (queue.size() >= 1)
	{
		auto top = queue.top();
		queue.pop();
		if (distance.count(top.node) == 1)
		{
			assert(distance.at(top.node) <= top.distance);
			continue;
		}
		distance[top.node] = top.distance;
		if (top.node.id() != middleNode)
		{
			queue.emplace(top.distance + graph.nodeSeqs.at(top.node.id()).size(), reverse(top.node));
		}
		if (edges.count(top.node) == 1)
		{
			for (auto edge : edges.at(top.node))
			{
				if (nodes.count(edge.id() == 0)) continue;
				queue.emplace(top.distance, reverse(edge));
			}
		}
	}
	if (distance.count(Node { middleNode, false }) == 0) return std::numeric_limits<size_t>::max();
	return distance.at(Node { middleNode, false });
}

void filterOutDefinitelyNonConsensusNodes(std::unordered_set<size_t>& nodes, std::unordered_map<Node, std::unordered_set<Node>>& edges, const size_t startNode, const GfaGraph& graph, const size_t localResolveLength)
{
	std::unordered_set<size_t> coreNodes;
	std::vector<size_t> testNodes { nodes.begin(), nodes.end() };
	for (size_t node : testNodes)
	{
		if (node == startNode) continue;
		nodes.erase(node);
		if (!isCircular(nodes, edges, Node { startNode, true }))
		{
			size_t splitDistance = splitUndirectedDistance(graph, nodes, edges, node);
			if (splitDistance >= localResolveLength)
			{
				coreNodes.insert(node);
			}
		}
		nodes.insert(node);
	}
	if (coreNodes.size() < 2) return;
	nodes.erase(startNode);
	if (!isCircular(nodes, edges, Node { *coreNodes.begin(), true })) coreNodes.insert(startNode);
	nodes.insert(startNode);
	std::vector<Node> coreNodeOrder = getCoreNodeOrder(nodes, edges, coreNodes, (coreNodes.count(startNode) == 1) ? startNode : *coreNodes.begin());
	for (size_t i = 1; i < coreNodeOrder.size(); i++)
	{
		filterOutDefinitelyNonConsensusNodes(nodes, edges, coreNodeOrder[i-1], coreNodeOrder[i], graph);
	}
}

Path getHeavyPath(const GfaGraph& graph, const size_t localResolveLength)
{
	std::unordered_set<size_t> coveredNodes;
	std::unordered_map<Node, std::unordered_set<Node>> coveredEdges;
	std::vector<std::pair<size_t, size_t>> nodeAndCoverage;
	std::vector<std::tuple<Node, Node, size_t>> edgeAndCoverage;
	for (size_t i = 0; i < graph.nodeCoverages.size(); i++)
	{
		nodeAndCoverage.emplace_back(i, graph.nodeCoverages[i]);
	}
	for (const auto& pair : graph.edges)
	{
		for (const auto& target : pair.second)
		{
			edgeAndCoverage.emplace_back(pair.first, std::get<0>(target), std::get<2>(target));
		}
	}
	std::sort(nodeAndCoverage.begin(), nodeAndCoverage.end(), [](auto left, auto right) { return left.second > right.second; });
	std::sort(edgeAndCoverage.begin(), edgeAndCoverage.end(), [](auto left, auto right) { return std::get<2>(left) > std::get<2>(right); });
	size_t topCoverageNode = nodeAndCoverage[0].first;
	size_t minCoverage = nodeAndCoverage[0].second;
	size_t nodeIndex = 0;
	size_t edgeIndex = 0;
	while (minCoverage > 0)
	{
		size_t nextCoverage = 0;
		while (nodeIndex < nodeAndCoverage.size() && nodeAndCoverage[nodeIndex].second >= minCoverage)
		{
			coveredNodes.emplace(nodeAndCoverage[nodeIndex].first);
			nodeIndex += 1;
		}
		while (edgeIndex < edgeAndCoverage.size() && std::get<2>(edgeAndCoverage[edgeIndex]) >= minCoverage)
		{
			coveredEdges[std::get<0>(edgeAndCoverage[edgeIndex])].emplace(std::get<1>(edgeAndCoverage[edgeIndex]));
			coveredEdges[reverse(std::get<1>(edgeAndCoverage[edgeIndex]))].emplace(reverse(std::get<0>(edgeAndCoverage[edgeIndex])));
			edgeIndex += 1;
		}
		if (nodeIndex < nodeAndCoverage.size()) nextCoverage = nodeAndCoverage[nodeIndex].second;
		if (edgeIndex < edgeAndCoverage.size()) nextCoverage = std::max(nextCoverage, std::get<2>(edgeAndCoverage[edgeIndex]));
		assert(nextCoverage < minCoverage);
		minCoverage = nextCoverage;
		if (isCircular(coveredNodes, coveredEdges, Node { topCoverageNode, true })) break;
	}
	filterOutNonCircularParts(coveredNodes, coveredEdges, topCoverageNode);
	while (true)
	{
		size_t sizeBeforeFilter = coveredNodes.size();
		filterOutDefinitelyNonConsensusNodes(coveredNodes, coveredEdges, topCoverageNode, graph, localResolveLength);
		filterOutNonCircularParts(coveredNodes, coveredEdges, topCoverageNode);
		if (coveredNodes.size() == sizeBeforeFilter) break;
	}
	Path result = getRecWidestPath(coveredNodes, coveredEdges, graph, Node { topCoverageNode, true });
	for (size_t i = 0; i < result.nodes.size(); i++)
	{
		Node prev;
		if (i > 0)
		{
			prev = result.nodes[i-1];
		}
		else
		{
			prev = result.nodes.back();
		}
		assert(graph.edges.count(prev) == 1);
		for (auto edge : graph.edges.at(prev))
		{
			if (std::get<0>(edge) != result.nodes[i]) continue;
			result.overlaps.push_back(std::get<1>(edge));
		}
	}
	assert(result.overlaps.size() == result.nodes.size());
	result.leftClip = 0;
	result.rightClip = result.overlaps[0];
	return result;
}

void writePathSequence(const Path& path, const GfaGraph& graph, std::string outputFile, const std::string& namePrefix)
{
	std::string pathseq = path.getSequence(graph.nodeSeqs);
	std::ofstream file { outputFile };
	file << ">";
	if (namePrefix != "") file << namePrefix << "_";
	file << "heavy_path" << std::endl;
	file << pathseq << std::endl;
}

void writePathGaf(const Path& path, const GfaGraph& graph, std::string outputFile)
{
	std::ofstream file { outputFile };
	file << getPathGaf(path, graph);
}

void runMBG(std::string basePath, std::string readPath, std::string MBGPath, const size_t k, const size_t maxResolveLength, const size_t numThreads)
{
	std::string mbgCommand;
	mbgCommand = MBGPath + " -o " + basePath + "/graph.gfa -t " + std::to_string(numThreads) + " -i " + readPath + " -k " + std::to_string(k) + " -w " + std::to_string(k-30) + " -a 2 -u 3 -r " + std::to_string(maxResolveLength) + " -R 4000 --error-masking=msat --output-sequence-paths " + basePath + "/paths.gaf --only-local-resolve 1> " + basePath + "/mbg_stdout.txt 2> " + basePath + "/mbg_stderr.txt";
	std::cerr << "MBG command:" << std::endl;
	std::cerr << mbgCommand << std::endl;
	int result = system(mbgCommand.c_str());
	if (result != 0)
	{
		std::cerr << "MBG did not run successfully" << std::endl;
		std::abort();
	}
}

std::vector<ReadPath> loadReadPaths(const std::string& filename, const GfaGraph& graph)
{
	std::vector<ReadPath> result;
	std::ifstream file { filename };
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (!file.good()) break;
		auto parts = split(line, '\t');
		ReadPath here;
		std::string readname = parts[0];
		size_t readlength = std::stoull(parts[1]);
		size_t readstart = std::stoull(parts[2]);
		size_t readend = std::stoull(parts[3]);
		std::string pathstr = parts[5];
		size_t pathLength = std::stoull(parts[6]);
		size_t pathStart = std::stoull(parts[7]);
		size_t pathEnd = std::stoull(parts[8]);
		here.readLength = readlength;
		here.readStart = readstart;
		here.readEnd = readend;
		here.pathStartClip = pathStart;
		here.pathEndClip = pathLength - pathEnd;
		here.path = parsePath(pathstr, graph.nodeNameToId);
		assert(here.path.size() >= 1);
		result.emplace_back(here);
	}
	return result;
}

std::vector<Node> getReferenceAllele(const Path& heavyPath, const Node startNode, const Node endNode)
{
	size_t startIndex = heavyPath.nodes.size();
	size_t endIndex = heavyPath.nodes.size();
	for (size_t i = 0; i < heavyPath.nodes.size(); i++)
	{
		if (i == heavyPath.nodes.size()-1 && heavyPath.nodes[0] == heavyPath.nodes.back()) continue;
		if (heavyPath.nodes[i] == startNode) startIndex = i;
		if (heavyPath.nodes[i] == endNode) endIndex = i;
	}
	assert(startIndex < heavyPath.nodes.size());
	assert(endIndex < heavyPath.nodes.size());
	std::vector<Node> result;
	if (endIndex > startIndex)
	{
		result.insert(result.end(), heavyPath.nodes.begin() + startIndex, heavyPath.nodes.begin() + endIndex + 1);
	}
	else
	{
		result.insert(result.end(), heavyPath.nodes.begin() + startIndex, heavyPath.nodes.end());
		assert(result.size() >= 1);
		if (heavyPath.nodes[0] == heavyPath.nodes.back())
		{
			assert(result.back() == heavyPath.nodes[0]);
			result.pop_back();
		}
		result.insert(result.end(), heavyPath.nodes.begin(), heavyPath.nodes.begin() + endIndex + 1);
	}
	return result;
}

bool isSubstring(const std::vector<Node>& small, const std::vector<Node>& big)
{
	if (small.size() > big.size()) return false;
	for (size_t i = 0; i <= big.size() - small.size(); i++)
	{
		bool match = true;
		for (size_t j = 0; j < small.size(); j++)
		{
			if (big[i+j] != small[j])
			{
				match = false;
				break;
			}
		}
		if (match)
		{
			return true;
		}
	}
	return false;
}

size_t getReferenceCoverage(const std::vector<ReadPath>& readPaths, const std::vector<Node>& referenceAllele)
{
	size_t result = 0;
	auto reverseRefAllele = reverse(referenceAllele);
	for (const auto& readPath : readPaths)
	{
		if (isSubstring(referenceAllele, readPath.path))
		{
			result += 1;
			continue;
		}
		if (isSubstring(reverseRefAllele, readPath.path))
		{
			result += 1;
			continue;
		}
	}
	return result;
}

std::vector<Variant> getVariants(const GfaGraph& graph, const Path& heavyPath, const std::vector<ReadPath>& readPaths, const size_t minCoverage)
{
	std::vector<bool> partOfHeavyPath;
	std::vector<bool> nodeOrientation;
	partOfHeavyPath.resize(graph.numNodes(), false);
	nodeOrientation.resize(graph.numNodes(), false);
	for (auto node : heavyPath.nodes)
	{
		partOfHeavyPath[node.id()] = true;
		nodeOrientation[node.id()] = node.forward();
	}
	std::map<std::vector<Node>, size_t> bubbleCoverages;
	bool hasPalindromicVariants = false;
	for (const auto& path : readPaths)
	{
		size_t lastStart = std::numeric_limits<size_t>::max();
		bool lastMatchesReferenceOrientation = true;
		for (size_t i = 0; i < path.path.size(); i++)
		{
			if (!partOfHeavyPath.at(path.path[i].id())) continue;
			bool thisMatchesReferenceOrientation = (path.path[i].forward()) == nodeOrientation.at(path.path[i].id());
			if (lastStart != std::numeric_limits<size_t>::max())
			{
				if (lastMatchesReferenceOrientation == thisMatchesReferenceOrientation)
				{
					std::vector<Node> part { path.path.begin() + lastStart, path.path.begin() + i + 1 };
					part = canon(part);
					bubbleCoverages[part] += 1;
				}
				else
				{
					hasPalindromicVariants = true;
				}
			}
			lastStart = i;
			lastMatchesReferenceOrientation = thisMatchesReferenceOrientation;
		}
	}
	if (hasPalindromicVariants)
	{
		std::cerr << "Note: the genome has palindromic variants." << std::endl;
		std::cerr << "Palindromic variants are ignored and missing from output." << std::endl;
	}
	std::vector<Variant> result;
	for (const auto& pair : bubbleCoverages)
	{
		if (pair.second < minCoverage) continue;
		std::vector<Node> path = pair.first;
		if (nodeOrientation[path[0].id()] != path[0].forward())
		{
			path = reverse(path);
		}
		std::vector<Node> referenceAllele = getReferenceAllele(heavyPath, path[0], path.back());
		if (path == referenceAllele) continue;
		result.emplace_back();
		result.back().path = path;
		result.back().coverage = pair.second;
		result.back().referencePath = referenceAllele;
		result.back().referenceCoverage = getReferenceCoverage(readPaths, result.back().referencePath);
	}
	return result;
}

void writeVariants(const Path& heavyPath, const GfaGraph& graph, const std::vector<Variant>& variants, const std::string& filename)
{
	std::ofstream file { filename };
	for (const auto& variant : variants)
	{
		file << variant.name << "\t";
		for (const auto& node : variant.path)
		{
			file << (node.forward() ? ">" : "<") << graph.nodeNames.at(node.id());
		}
		file << "\t";
		for (const auto& node : variant.referencePath)
		{
			file << (node.forward() ? ">" : "<") << graph.nodeNames.at(node.id());
		}
		file << "\t" << variant.coverage << "\t" << variant.referenceCoverage;
		file << "\t" << variant.variantSeq << "\t" << variant.referenceSeq << std::endl;
	}
}

bool alleleMatchLeft(const std::vector<Variant>& variants, size_t lefti, size_t leftj, size_t righti, size_t rightj)
{
	if (leftj != rightj) return false;
	for (size_t k = 0; k < leftj; k++)
	{
		if (variants[lefti].path[k] != variants[righti].path[k]) return false;
	}
	return true;
}

bool alleleMatchRight(const std::vector<Variant>& variants, size_t lefti, size_t leftj, size_t righti, size_t rightj)
{
	if (variants[lefti].path.size()-leftj != variants[righti].path.size()-rightj) return false;
	for (size_t k = 0; k < variants[lefti].path.size()-leftj; k++)
	{
		if (variants[lefti].path[leftj+k] != variants[righti].path[rightj+k]) return false;
	}
	return true;
}

void writeAlleleGraph(const GfaGraph& fullGraph, const Path& heavyPath, const std::vector<Variant>& variants, std::string alleleGraphFileName)
{
	std::vector<bool> referenceNode;
	referenceNode.resize(fullGraph.numNodes(), false);
	for (size_t i = 0; i < heavyPath.nodes.size(); i++)
	{
		referenceNode[heavyPath.nodes[i].id()] = true;
	}
	std::ofstream out { alleleGraphFileName };
	for (size_t i = 0; i < fullGraph.numNodes(); i++)
	{
		if (!referenceNode[i]) continue;
		out << "S\t" << fullGraph.nodeNames[i] << "\t" << fullGraph.nodeSeqs[i] << "\tll:f:" << fullGraph.nodeCoverages[i] << "\tFC:f:" << (fullGraph.nodeCoverages[i] * fullGraph.nodeSeqs[i].size()) << std::endl;
	}
	for (size_t i = 0; i < heavyPath.nodes.size(); i++)
	{
		if (i == 0 && heavyPath.nodes[0] == heavyPath.nodes.back()) continue;
		Node prev;
		if (i == 0)
		{
			prev = heavyPath.nodes.back();
		}
		else
		{
			prev = heavyPath.nodes[i-1];
		}
		Node curr = heavyPath.nodes[i];
		assert(fullGraph.edges.count(prev) == 1);
		size_t overlap = std::numeric_limits<size_t>::max();
		size_t coverage = std::numeric_limits<size_t>::max();
		for (auto edge : fullGraph.edges.at(prev))
		{
			if (std::get<0>(edge) != curr) continue;
			assert(overlap == std::numeric_limits<size_t>::max());
			overlap = std::get<1>(edge);
			coverage = std::get<2>(edge);
		}
		assert(overlap != std::numeric_limits<size_t>::max());
		out << "L\t" << fullGraph.nodeNames[prev.id()] << "\t" << (prev.forward() ? "+" : "-") << "\t" << fullGraph.nodeNames[curr.id()] << "\t" << (curr.forward() ? "+" : "-") << "\t" << overlap << "M" << "\tec:i:" << coverage << std::endl;
	}
	std::vector<std::vector<std::pair<size_t, size_t>>> parent;
	parent.resize(variants.size());
	std::vector<std::vector<std::pair<size_t, size_t>>> nodePositions;
	nodePositions.resize(fullGraph.numNodes());
	for (size_t i = 0; i < variants.size(); i++)
	{
		parent[i].resize(variants[i].path.size());
		for (size_t j = 0; j < variants[i].path.size(); j++)
		{
			parent[i][j] = std::make_pair(i, j);
			nodePositions[variants[i].path[j].id()].emplace_back(i, j);
		}
	}
	for (size_t i = 0; i < nodePositions.size(); i++)
	{
		if (nodePositions[i].size() < 2) continue;
		for (size_t j = 0; j < nodePositions[i].size(); j++)
		{
			for (size_t k = 0; k < j; k++)
			{
				auto leftPos = nodePositions[i][j];
				auto rightPos = nodePositions[i][k];
				if (alleleMatchLeft(variants, leftPos.first, leftPos.second, rightPos.first, rightPos.second) || alleleMatchRight(variants, leftPos.first, leftPos.second, rightPos.first, rightPos.second))
				{
					while (parent[leftPos.first][leftPos.second] != parent[parent[leftPos.first][leftPos.second].first][parent[leftPos.first][leftPos.second].second])
					{
						parent[leftPos.first][leftPos.second] = parent[parent[leftPos.first][leftPos.second].first][parent[leftPos.first][leftPos.second].second];
					}
					while (parent[rightPos.first][rightPos.second] != parent[parent[rightPos.first][rightPos.second].first][parent[rightPos.first][rightPos.second].second])
					{
						parent[rightPos.first][rightPos.second] = parent[parent[rightPos.first][rightPos.second].first][parent[rightPos.first][rightPos.second].second];
					}
					parent[rightPos.first][rightPos.second] = parent[leftPos.first][leftPos.second];
				}
			}
		}
	}
	std::map<std::pair<size_t, size_t>, size_t> nodeCoverage;
	for (size_t i = 0; i < variants.size(); i++)
	{
		for (size_t j = 1; j < variants[i].path.size()-1; j++)
		{
			auto key = parent[i][j];
			while (parent[key.first][key.second] != key) key = parent[key.first][key.second];
			nodeCoverage[key] += variants[i].coverage;
		}
	}
	std::map<std::pair<std::string, std::string>, size_t> edgeCoverage;
	for (size_t i = 0; i < variants.size(); i++)
	{
		for (size_t j = 1; j < variants[i].path.size(); j++)
		{
			std::string fromname = fullGraph.nodeNames[variants[i].path[j-1].id()];
			if (j != 1)
			{
				auto key = parent[i][j-1];
				while (parent[key.first][key.second] != key) key = parent[key.first][key.second];
				fromname += "_" + std::to_string(key.first) + "_" + std::to_string(key.second);
			}
			std::string toname = fullGraph.nodeNames[variants[i].path[j].id()];
			if (j != variants[i].path.size()-1)
			{
				auto key = parent[i][j];
				while (parent[key.first][key.second] != key) key = parent[key.first][key.second];
				toname += "_" + std::to_string(key.first) + "_" + std::to_string(key.second);
			}
			edgeCoverage[std::make_pair(fromname, toname)] += variants[i].coverage;
		}
	}
	for (size_t i = 0; i < variants.size(); i++)
	{
		assert(variants[i].path.size() >= 2);
		if (variants[i].path.size() == 2)
		{
			Node prev = variants[i].path[0];
			Node curr = variants[i].path[1];
			assert(fullGraph.edges.count(prev) == 1);
			size_t overlap = std::numeric_limits<size_t>::max();
			for (auto edge : fullGraph.edges.at(prev))
			{
				if (std::get<0>(edge) != curr) continue;
				assert(overlap == std::numeric_limits<size_t>::max());
				overlap = std::get<1>(edge);
			}
			assert(overlap != std::numeric_limits<size_t>::max());
			out << "L\t" << fullGraph.nodeNames[prev.id()] << "\t" << (prev.forward() ? "+" : "-") << "\t" << fullGraph.nodeNames[curr.id()] << "\t" << (curr.forward() ? "+" : "-") << "\t" << overlap << "M" << "\tec:i:" << variants[i].coverage << std::endl;
			continue;
		}
		assert(referenceNode[variants[i].path[0].id()]);
		assert(referenceNode[variants[i].path.back().id()]);
		for (size_t j = 1; j < variants[i].path.size()-1; j++)
		{
			std::pair<size_t, size_t> key = parent[i][j];
			while (parent[key.first][key.second] != key)
			{
				key = parent[key.first][key.second];
			}
			if (key.first == i && key.second == j)
			{
				out << "S\t" << fullGraph.nodeNames[variants[i].path[j].id()] << "_" << key.first << "_" << key.second << "\t" << fullGraph.nodeSeqs[variants[i].path[j].id()] << "\tll:f:" << nodeCoverage.at(key) << "\tFC:f:" << (nodeCoverage.at(key) * fullGraph.nodeSeqs[variants[i].path[j].id()].size()) << std::endl;
			}
		}
		{
			Node prev = variants[i].path[0];
			Node curr = variants[i].path[1];
			std::pair<size_t, size_t> key = parent[i][1];
			while (parent[key.first][key.second] != key)
			{
				key = parent[key.first][key.second];
			}
			assert(fullGraph.edges.count(prev) == 1);
			size_t overlap = std::numeric_limits<size_t>::max();
			for (auto edge : fullGraph.edges.at(prev))
			{
				if (std::get<0>(edge) != curr) continue;
				assert(overlap == std::numeric_limits<size_t>::max());
				overlap = std::get<1>(edge);
			}
			assert(overlap != std::numeric_limits<size_t>::max());
			std::string fromname = fullGraph.nodeNames[prev.id()];
			std::string toname = fullGraph.nodeNames[curr.id()] + "_" + std::to_string(key.first) + "_" + std::to_string(key.second);
			out << "L\t" << fromname << "\t" << (prev.forward() ? "+" : "-") << "\t" << toname << "\t" << (curr.forward() ? "+" : "-") << "\t" << overlap << "M" << "\tec:i:" << edgeCoverage.at(std::make_pair(fromname, toname)) << std::endl;
		}
		{
			Node prev = variants[i].path[variants[i].path.size()-2];
			Node curr = variants[i].path[variants[i].path.size()-1];
			std::pair<size_t, size_t> key = parent[i][variants[i].path.size()-2];
			while (parent[key.first][key.second] != key)
			{
				key = parent[key.first][key.second];
			}
			assert(fullGraph.edges.count(prev) == 1);
			size_t overlap = std::numeric_limits<size_t>::max();
			for (auto edge : fullGraph.edges.at(prev))
			{
				if (std::get<0>(edge) != curr) continue;
				assert(overlap == std::numeric_limits<size_t>::max());
				overlap = std::get<1>(edge);
			}
			assert(overlap != std::numeric_limits<size_t>::max());
			std::string fromname = fullGraph.nodeNames[prev.id()] + "_" + std::to_string(key.first) + "_" + std::to_string(key.second);
			std::string toname = fullGraph.nodeNames[curr.id()];
			out << "L\t" << fromname << "\t" << (prev.forward() ? "+" : "-") << "\t" << toname << "\t" << (curr.forward() ? "+" : "-") << "\t" << overlap << "M" << "\tec:i:" << edgeCoverage.at(std::make_pair(fromname, toname)) << std::endl;
		}
		for (size_t j = 2; j < variants[i].path.size()-1; j++)
		{
			Node prev = variants[i].path[j-1];
			Node curr = variants[i].path[j];
			std::pair<size_t, size_t> fromkey = parent[i][j-1];
			while (parent[fromkey.first][fromkey.second] != fromkey)
			{
				fromkey = parent[fromkey.first][fromkey.second];
			}
			std::pair<size_t, size_t> tokey = parent[i][j];
			while (parent[tokey.first][tokey.second] != tokey)
			{
				tokey = parent[tokey.first][tokey.second];
			}
			assert(fullGraph.edges.count(prev) == 1);
			size_t overlap = std::numeric_limits<size_t>::max();
			for (auto edge : fullGraph.edges.at(prev))
			{
				if (std::get<0>(edge) != curr) continue;
				assert(overlap == std::numeric_limits<size_t>::max());
				overlap = std::get<1>(edge);
			}
			assert(overlap != std::numeric_limits<size_t>::max());
			std::string fromname = fullGraph.nodeNames[prev.id()] + "_" + std::to_string(fromkey.first) + "_" + std::to_string(fromkey.second);
			std::string toname = fullGraph.nodeNames[curr.id()] + "_" + std::to_string(tokey.first) + "_" + std::to_string(tokey.second);
			out << "L\t" << fromname << "\t" << (prev.forward() ? "+" : "-") << "\t" << toname << "\t" << (curr.forward() ? "+" : "-") << "\t" << overlap << "M" << "\tec:i:" << edgeCoverage.at(std::make_pair(fromname, toname)) << std::endl;
		}
	}
}

void writeVariantGraph(std::string graphFileName, const GfaGraph& fullGraph, const Path& heavyPath, const std::vector<Variant>& variants, std::string variantGraphFileName)
{
	std::vector<bool> allowedNodes;
	allowedNodes.resize(fullGraph.numNodes(), false);
	std::set<std::pair<Node, Node>> allowedEdges;
	for (size_t i = 0; i < heavyPath.nodes.size(); i++)
	{
		Node prev;
		if (i == 0)
		{
			prev = heavyPath.nodes.back();
		}
		else
		{
			prev = heavyPath.nodes[i-1];
		}
		Node curr = heavyPath.nodes[i];
		allowedEdges.emplace(prev, curr);
		allowedEdges.emplace(reverse(curr), reverse(prev));
		allowedNodes[curr.id()] = true;
	}
	for (const auto& variant : variants)
	{
		for (size_t i = 0; i < variant.path.size(); i++)
		{
			Node prev;
			if (i == 0)
			{
				prev = variant.path.back();
			}
			else
			{
				prev = variant.path[i-1];
			}
			Node curr = variant.path[i];
			allowedEdges.emplace(prev, curr);
			allowedEdges.emplace(reverse(curr), reverse(prev));
			allowedNodes[curr.id()] = true;
		}
	}
	std::ifstream in { graphFileName };
	std::ofstream out { variantGraphFileName };
	while (in.good())
	{
		std::string line;
		getline(in, line);
		std::stringstream sstr { line };
		std::string test;
		sstr >> test;
		if (test == "S")
		{
			std::string nodename;
			sstr >> nodename;
			if (allowedNodes[fullGraph.nodeNameToId.at(nodename)])
			{
				out << line << std::endl;
			}
		}
		else if (test == "L")
		{
			std::string fromnodename, fromorient, tonodename, toorient;
			sstr >> fromnodename >> fromorient >> tonodename >> toorient;
			Node fromnode { fullGraph.nodeNameToId.at(fromnodename), fromorient == "+" };
			Node tonode { fullGraph.nodeNameToId.at(tonodename), toorient == "+" };
			if (allowedEdges.count(std::make_pair(fromnode, tonode)) == 1)
			{
				out << line << std::endl;
			}
		}
	}
}

Path orientPath(const GfaGraph& graph, const Path& rawPath, const std::string& orientReferencePath, size_t k)
{
	Path result = rawPath;
	std::unordered_map<uint64_t, size_t> kmerReferencePosition;
	FastQ::streamFastqFromFile(orientReferencePath, false, [&kmerReferencePosition, k](FastQ& fastq)
	{
		for (size_t i = 0; i < fastq.sequence.size()-k; i++)
		{
			uint64_t hash = std::hash<std::string>{}(fastq.sequence.substr(i, k));
			if (kmerReferencePosition.count(hash) == 1)
			{
				kmerReferencePosition[hash] = std::numeric_limits<size_t>::max();
			}
			else
			{
				kmerReferencePosition[hash] = i;
			}
		}
	});
	{
		std::string pathSeq = result.getSequence(graph.nodeSeqs);
		std::string reverseSeq = revcomp(pathSeq);
		size_t fwMatches = 0;
		size_t bwMatches = 0;
		for (size_t i = 0; i < pathSeq.size() - k; i++)
		{
			uint64_t hash = std::hash<std::string>{}(pathSeq.substr(i, k));
			if (kmerReferencePosition.count(hash) == 1) fwMatches += 1; // no need to check if the hash is unique, it's still valid for orientation
		}
		for (size_t i = 0; i < reverseSeq.size() - k; i++)
		{
			uint64_t hash = std::hash<std::string>{}(reverseSeq.substr(i, k));
			if (kmerReferencePosition.count(hash) == 1) bwMatches += 1; // no need to check if the hash is unique, it's still valid for orientation
		}
		if (fwMatches == 0 && bwMatches == 0)
		{
			std::cerr << "Can't be matched to the reference!" << std::endl;
			return rawPath;
		}
		if (bwMatches > fwMatches)
		{
			std::cerr << "reverse complement" << std::endl;
			std::reverse(result.nodes.begin(), result.nodes.end());
			for (size_t i = 0; i < result.nodes.size(); i++)
			{
				result.nodes[i] = reverse(result.nodes[i]);
			}
			std::reverse(result.overlaps.begin(), result.overlaps.end());
			result.overlaps.insert(result.overlaps.begin(), result.overlaps.back());
			result.overlaps.pop_back();
		}
	}
	std::string pathSeq = result.getSequence(graph.nodeSeqs);
	// todo: this does not work in general due to spurious k-mer hits, but seems to be good enough for now?
	// should instead do chaining between matches and pick leftmost ref pos from them
	std::pair<size_t, size_t> minMatchPos { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
	for (size_t i = 0; i < pathSeq.size() - k; i++)
	{
		uint64_t hash = std::hash<std::string>{}(pathSeq.substr(i, k));
		if (kmerReferencePosition.count(hash) == 0) continue;
		auto refPos = kmerReferencePosition.at(hash);
		if (refPos == std::numeric_limits<size_t>::max()) continue;
		if (refPos < minMatchPos.second)
		{
			minMatchPos.first = i;
			minMatchPos.second = refPos;
		}
	}
	assert(minMatchPos.first != std::numeric_limits<size_t>::max());
	size_t rotatePosition = 0;
	if (minMatchPos.first >= minMatchPos.second)
	{
		rotatePosition = minMatchPos.first - minMatchPos.second;
	}
	else
	{
		assert(pathSeq.size() + minMatchPos.first >= minMatchPos.second);
		rotatePosition = pathSeq.size() + minMatchPos.first - minMatchPos.second;
	}
	size_t rotateNodes = 0;
	size_t pos = 0;
	size_t extraRotate = 0;
	for (size_t i = 1; i < result.nodes.size(); i++)
	{
		pos += graph.nodeSeqs.at(result.nodes[i-1].id()).size();
		pos -= result.overlaps[i];
		if (pos > rotatePosition)
		{
			break;
		}
		rotateNodes = i;
		assert(rotatePosition >= pos);
		extraRotate = rotatePosition - pos;
	}
	assert(rotateNodes != std::numeric_limits<size_t>::max());
	assert(rotateNodes < result.nodes.size());
	std::cerr << "rotate by " << rotateNodes << " nodes and " << extraRotate << " base pairs" << std::endl;
	if (rotateNodes != 0)
	{
		Path result2;
		result2.nodes.insert(result2.nodes.end(), result.nodes.begin() + rotateNodes, result.nodes.end());
		result2.nodes.insert(result2.nodes.end(), result.nodes.begin(), result.nodes.begin() + rotateNodes);
		result2.overlaps.insert(result2.overlaps.end(), result.overlaps.begin() + rotateNodes, result.overlaps.end());
		result2.overlaps.insert(result2.overlaps.end(), result.overlaps.begin(), result.overlaps.begin() + rotateNodes);
		result = result2;
	}
	result.leftClip = extraRotate;
	if (result.overlaps[0] <= result.leftClip)
	{
		result.nodes.push_back(result.nodes[0]);
		result.overlaps.push_back(result.overlaps[0]);
		result.rightClip = graph.nodeSeqs.at(result.nodes.back().id()).size() - extraRotate;
	}
	else
	{
		result.rightClip = result.overlaps[0] - result.leftClip;
	}
	return result;
}

void nameVariants(std::vector<Variant>& variants, const GfaGraph& graph, const Path& heavyPath)
{
	size_t pathLength = heavyPath.getSequence(graph.nodeSeqs).size();
	for (size_t variant = 0; variant < variants.size(); variant++)
	{
		Node startNode = variants[variant].referencePath[0];
		Node endNode = variants[variant].referencePath.back();
		size_t startIndex = std::numeric_limits<size_t>::max();
		size_t endIndex = std::numeric_limits<size_t>::max();
		size_t startPos = std::numeric_limits<size_t>::max();
		size_t endPos = std::numeric_limits<size_t>::max();
		size_t pathPos = 0;
		for (size_t i = 0; i < heavyPath.nodes.size(); i++)
		{
			if (i > 0)
			{
				pathPos += graph.nodeSeqs.at(heavyPath.nodes[i-1].id()).size() - heavyPath.overlaps[i];
			}
			if (i == heavyPath.nodes.size()-1 && heavyPath.nodes[0] == heavyPath.nodes.back()) break;
			if (heavyPath.nodes[i] == startNode)
			{
				startIndex = i;
				startPos = pathPos;
			}
			if (heavyPath.nodes[i] == endNode)
			{
				endIndex = i;
				endPos = pathPos + graph.nodeSeqs.at(heavyPath.nodes[i].id()).size();
			}
		}
		assert(startIndex != std::numeric_limits<size_t>::max());
		assert(endIndex != std::numeric_limits<size_t>::max());
		if (startPos >= heavyPath.leftClip)
		{
			variants[variant].referenceStartPos = startPos - heavyPath.leftClip;
		}
		else
		{
			variants[variant].referenceStartPos = pathLength + startPos - heavyPath.leftClip;
		}
		if (endPos >= heavyPath.leftClip)
		{
			variants[variant].referenceEndPos = endPos - heavyPath.leftClip;
		}
		else
		{
			variants[variant].referenceEndPos = pathLength + endPos - heavyPath.leftClip;
		}
		assert(variants[variant].referenceStartPos < pathLength);
		variants[variant].referenceEndPos = variants[variant].referenceEndPos % pathLength;
		std::string referenceSeq = getSequence(variants[variant].referencePath, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
		std::string variantSeq = getSequence(variants[variant].path, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
		assert(variants[variant].referenceEndPos < variants[variant].referenceStartPos || referenceSeq.size() % pathLength == variants[variant].referenceEndPos - variants[variant].referenceStartPos);
		assert(variants[variant].referenceEndPos > variants[variant].referenceStartPos || referenceSeq.size() % pathLength == pathLength + variants[variant].referenceEndPos - variants[variant].referenceStartPos);
		size_t leftClip = 0;
		size_t rightClip = 0;
		if (graph.nodeSeqs.at(variants[variant].referencePath[0].id()).size()+graph.nodeSeqs.at(variants[variant].referencePath.back().id()).size()+2 >= referenceSeq.size())
		{
			while (leftClip < variantSeq.size() && leftClip < referenceSeq.size() && leftClip < graph.nodeSeqs.at(variants[variant].referencePath[0].id()).size() && variantSeq[leftClip] == referenceSeq[leftClip]) leftClip += 1;
			while (rightClip < variantSeq.size() && rightClip < referenceSeq.size() && rightClip < graph.nodeSeqs.at(variants[variant].referencePath.back().id()).size() && variantSeq[variantSeq.size()-1-rightClip] == referenceSeq[referenceSeq.size()-1-rightClip]) rightClip += 1;
			assert(leftClip == graph.nodeSeqs.at(variants[variant].referencePath[0].id()).size());
			assert(rightClip == graph.nodeSeqs.at(variants[variant].referencePath.back().id()).size());
		}
		while (leftClip + rightClip < variantSeq.size() && leftClip + rightClip < referenceSeq.size() && variantSeq[leftClip] == referenceSeq[leftClip]) leftClip += 1;
		while (leftClip + rightClip < variantSeq.size() && leftClip + rightClip < referenceSeq.size() && variantSeq[variantSeq.size()-1-rightClip] == referenceSeq[referenceSeq.size()-1-rightClip]) rightClip += 1;
		if (leftClip + rightClip == variantSeq.size() || leftClip + rightClip == referenceSeq.size())
		{
			assert(leftClip >= 1);
			leftClip -= 1;
		}
		variants[variant].variantSeq = variantSeq.substr(leftClip, variantSeq.size()-leftClip-rightClip);
		variants[variant].referenceSeq = referenceSeq.substr(leftClip, referenceSeq.size()-leftClip-rightClip);
		if (variants[variant].variantSeq.size() == 0) variants[variant].variantSeq = ".";
		if (variants[variant].referenceSeq.size() == 0) variants[variant].referenceSeq = ".";
		if (variants[variant].referenceStartPos + leftClip < pathLength)
		{
			variants[variant].referenceStartPos += leftClip;
		}
		else
		{
			variants[variant].referenceStartPos += leftClip;
			variants[variant].referenceStartPos -= pathLength;
		}
		if (variants[variant].referenceEndPos >= rightClip)
		{
			variants[variant].referenceEndPos -= rightClip;
		}
		else
		{
			variants[variant].referenceEndPos = pathLength + variants[variant].referenceEndPos - rightClip;
		}
	}
	std::stable_sort(variants.begin(), variants.end(), [](const Variant& left, const Variant& right) { return left.coverage > right.coverage; });
	std::stable_sort(variants.begin(), variants.end(), [](const Variant& left, const Variant& right) { return left.referenceEndPos < right.referenceEndPos; });
	std::stable_sort(variants.begin(), variants.end(), [](const Variant& left, const Variant& right) { return left.referenceStartPos < right.referenceStartPos; });
	for (size_t variant = 0; variant < variants.size(); variant++)
	{
		variants[variant].name = "var" + std::to_string(variant) + "_" + std::to_string(variants[variant].referenceStartPos) + "_" + std::to_string(variants[variant].referenceEndPos);
	}
}

void writeVariantVCF(std::string filename, const Path& heavyPath, const GfaGraph& graph, const std::vector<Variant>& variants)
{
	size_t pathLength = heavyPath.getSequence(graph.nodeSeqs).size();
	std::ofstream file { filename };
	file << "##fileformat=VCFv4.2" << std::endl;
	file << "##contig=<ID=heavy_path,length=" << pathLength << ">" << std::endl;
	file << "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total read depth for each allele\">" << std::endl;
	file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;
	for (size_t i = 0; i < variants.size(); i++)
	{
		file << "heavy_path\t" << variants[i].referenceStartPos+1 << "\t" << variants[i].name << "\t" << variants[i].referenceSeq << "\t" << variants[i].variantSeq << "\t" << "." << "\tPASS\t" << "AD=" << variants[i].referenceCoverage << "," << variants[i].coverage << std::endl;
	}
}

void liftoverAnnotationsToMorphs(const std::string& basePath, const std::vector<MorphConsensus>& morphConsensuses, const std::string& annotationFasta, const std::string& annotationGff3)
{
	{
		std::ofstream typefile { basePath + "/liftoff_types.txt" };
		typefile << "rRNA" << std::endl;
		typefile << "misc_RNA" << std::endl;
		typefile << "repeat_region" << std::endl;
		typefile << "gene" << std::endl;
		typefile << "transcript" << std::endl;
		typefile << "exon" << std::endl;
		typefile << "pseudogene" << std::endl;
		typefile << "tandem_repeat" << std::endl;
	}
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		std::string tmpfilepath { basePath + "/tmp_seq.fa" };
		{
			std::ofstream tmpfile { tmpfilepath };
			tmpfile << ">" << morphConsensuses[i].name << std::endl;
			tmpfile << morphConsensuses[i].sequence << std::endl;
		}
		std::string command = "liftoff -f " + basePath + "/liftoff_types.txt -g " + annotationGff3 + " -o " + basePath + "/tmp-morph-annotations-part" + std::to_string(i) + ".gff3 -u "+ basePath + "/morph-unmapped_features" + std::to_string(i) + ".txt -dir " + basePath + "/liftoff_intermediate_files/ " + tmpfilepath + " " + annotationFasta + " 1> " + basePath + "/liftoff_morphs_stdout" + std::to_string(i) + ".txt 2> " + basePath + "/liftoff_morphs_stderr" + std::to_string(i) + ".txt";
		std::cerr << "running liftoff with command:" << std::endl;
		std::cerr << command << std::endl;
		int result = system(command.c_str());
		if (result != 0)
		{
			std::cerr << "liftoff did not run successfully" << std::endl;
			std::abort();
		}
		command = "rm -r " + basePath + "/liftoff_intermediate_files/ " + tmpfilepath + ".fai " + tmpfilepath + ".mmi";
		std::cerr << command << std::endl;
		result = system(command.c_str());
	}
	std::cerr << "combining liftoff results" << std::endl;
	std::string outputFile = basePath + "/morph-annotations.gff3";
	std::string command = "echo \"##gff-version 3\" > " + outputFile + "\ncat " + basePath + "/tmp-morph-annotations-part*.gff3 | grep -v '#' >> " + outputFile;
	std::cerr << command << std::endl;
	int result = system(command.c_str());
	if (result != 0)
	{
		std::cerr << "failed to combine liftoff results" << std::endl;
		std::abort();
	}
}

void liftoverAnnotationsToConsensus(const std::string& basePath, const std::string& consensusPath, const std::string& annotationFasta, const std::string& annotationGff3)
{
	{
		std::ofstream typefile { basePath + "/liftoff_types.txt" };
		typefile << "rRNA" << std::endl;
		typefile << "misc_RNA" << std::endl;
		typefile << "repeat_region" << std::endl;
		typefile << "gene" << std::endl;
		typefile << "transcript" << std::endl;
		typefile << "exon" << std::endl;
		typefile << "pseudogene" << std::endl;
		typefile << "tandem_repeat" << std::endl;
	}
	std::string command = "liftoff -f " + basePath + "/liftoff_types.txt -g " + annotationGff3 + " -o " + basePath + "/consensus-annotation.gff3 -u " + basePath + "/consensus-unmapped_features.txt -dir " + basePath + "/liftoff_intermediate_files/ " + consensusPath + " " + annotationFasta + " 1> " + basePath + "/liftoff_consensus_stdout.txt 2> " + basePath + "/liftoff_consensus_stderr.txt";
	std::cerr << "running liftoff with command:" << std::endl;
	std::cerr << command << std::endl;
	int result = system(command.c_str());
	if (result != 0)
	{
		std::cerr << "liftoff did not run successfully" << std::endl;
		std::abort();
	}
}

void AlignONTReads(std::string basePath, std::string graphAlignerPath, std::string ontReadPath, std::string graphPath, std::string outputAlnPath, size_t numThreads)
{
	std::string graphalignerCommand;
	graphalignerCommand = graphAlignerPath + " -g " + graphPath + " -f " + ontReadPath + " -a " + outputAlnPath + " -t " + std::to_string(numThreads) + " --seeds-mxm-length 30 --seeds-mem-count 10000 --bandwidth 15 --multimap-score-fraction 0.99 --precise-clipping 0.85 --min-alignment-score 5000 --discard-cigar --clip-ambiguous-ends 100 --overlap-incompatible-cutoff 0.15 --mem-index-no-wavelet-tree --max-trace-count 5 1> " + basePath + "/graphaligner_stdout.txt 2> " + basePath + "/graphaligner_stderr.txt";
	std::cerr << "GraphAligner command:" << std::endl;
	std::cerr << graphalignerCommand << std::endl;
	int result = system(graphalignerCommand.c_str());
	if (result != 0)
	{
		std::cerr << "GraphAligner did not run successfully" << std::endl;
		std::abort();
	}
}

std::string getPathSequence(const ReadPath& readPath, const GfaGraph& graph, const std::unordered_map<std::pair<Node, Node>, size_t>& edgeOverlaps)
{
	std::string result;
	for (size_t i = 0; i < readPath.path.size(); i++)
	{
		size_t overlap = 0;
		if (i > 0) overlap = edgeOverlaps.at(canon(readPath.path[i-1], readPath.path[i]));
		std::string add = graph.nodeSeqs.at(readPath.path[i].id());
		if (!readPath.path[i].forward())
		{
			add = revcomp(add);
		}
		if (overlap > 0) add = add.substr(overlap);
		result += add;
	}
	return result;
}

std::unordered_map<std::pair<Node, Node>, size_t> getEdgeOverlaps(const GfaGraph& graph)
{
	std::unordered_map<std::pair<Node, Node>, size_t> result;
	for (const auto& pair : graph.edges)
	{
		for (const auto& target : pair.second)
		{
			auto key = canon(pair.first, std::get<0>(target));
			size_t overlap = std::get<1>(target);
			assert(result.count(key) == 0 || result.at(key) == overlap);
			result[key] = overlap;
		}
	}
	return result;
}

std::unordered_map<size_t, bool> getNodeOrientations(const std::vector<Node>& referencePath)
{
	std::unordered_map<size_t, bool> result;
	for (const auto& str : referencePath)
	{
		assert(result.count(str.id()) == 0 || (str == referencePath[0] && str == referencePath.back()));
		result[str.id()] = str.forward();
	}
	return result;
}

bool orientPath(std::vector<Node>& path, const std::unordered_map<size_t, bool>& referenceOrientations)
{
	size_t fwMatches = 0;
	size_t bwMatches = 0;
	for (const auto& node : path)
	{
		if (referenceOrientations.count(node.id()) == 0) continue;
		bool orientHere = node.forward();
		bool refOrient = referenceOrientations.at(node.id());
		if (orientHere == refOrient)
		{
			fwMatches += 1;
		}
		else
		{
			bwMatches += 1;
		}
	}
	if (bwMatches > fwMatches)
	{
		path = reverse(path);
		return true;
	}
	return false;
}

std::vector<ReadPath> extractCorrectedONTPaths(std::string gafFile, const Path& heavyPath, const size_t minLength, const GfaGraph& graph)
{
	std::vector<ReadPath> result;
	std::ifstream file { gafFile };
	auto edgeOverlaps = getEdgeOverlaps(graph);
	auto nodeOrientations = getNodeOrientations(heavyPath.nodes);
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (line.size() == 0) continue;
		auto parts = split(line, '\t');
		ReadPath readPath;
		readPath.readName = parts[0];
		readPath.readLength = std::stoull(parts[1]);
		readPath.readStart = std::stoi(parts[2]);
		readPath.readEnd = std::stoi(parts[3]);
		readPath.pathStartClip = std::stoull(parts[7]);
		readPath.pathEndClip = std::stoull(parts[6]) - std::stoull(parts[8]);
		std::string pathstr = parts[5];
		size_t mapq = std::stoi(parts[11]);
		if (mapq < 20) continue;
		if (readPath.readEnd - readPath.readStart < minLength) continue;
		readPath.path = parsePath(pathstr, graph.nodeNameToId);
		bool reverse = orientPath(readPath.path, nodeOrientations);
		readPath.reverse = reverse;
		result.push_back(readPath);
	}
	return result;
}

std::unordered_set<size_t> getCoreNodes(const std::vector<OntLoop>& paths)
{
	std::unordered_set<size_t> existingNodes;
	std::unordered_set<size_t> notUniqueNodes;
	for (const auto& read : paths)
	{
		std::unordered_map<size_t, size_t> coverage;
		for (const auto& node : read.path)
		{
			coverage[node.id()] += 1;
		}
		for (auto pair : coverage)
		{
			if (pair.second == 1) existingNodes.insert(pair.first);
			if (pair.second != 1) notUniqueNodes.insert(pair.first);
		}
	}
	for (const auto& read : paths)
	{
		std::unordered_set<size_t> nodesHere;
		for (const auto& node : read.path) nodesHere.insert(node.id());
		for (auto node : existingNodes)
		{
			if (nodesHere.count(node) == 0) notUniqueNodes.insert(node);
		}
	}
	std::unordered_set<size_t> firstNodes;
	std::unordered_set<size_t> lastNodes;
	for (const auto& read : paths)
	{
		firstNodes.insert(read.path[0].id());
		lastNodes.insert(read.path.back().id());
	}
	std::unordered_set<size_t> potentialCoreNodes;
	for (auto node : existingNodes)
	{
		if (notUniqueNodes.count(node) == 1) continue;
		if (firstNodes.count(node) == 1 && lastNodes.count(node) == 1) continue;
		potentialCoreNodes.insert(node);
	}
	// static_assert(false, "todo: this did nothing before, check carefully");
	// std::unordered_map<size_t, size_t> uniqueCorePredecessor;
	// std::unordered_map<size_t, size_t> uniqueCoreSuccessor;
	// for (const auto& read : paths)
	// {
	// 	size_t lastCore = std::numeric_limits<size_t>::max();
	// 	for (const auto& node : read.path)
	// 	{
	// 		size_t nodename = node.id();
	// 		if (potentialCoreNodes.count(nodename) == 0) continue;
	// 		if (lastCore.size() == 0) continue;
	// 		if (uniqueCorePredecessor.count(nodename) == 1 && uniqueCorePredecessor.at(nodename) != lastCore) uniqueCorePredecessor[nodename] = std::numeric_limits<size_t>::max();
	// 		if (uniqueCoreSuccessor.count(lastCore) == 1 && uniqueCoreSuccessor.at(lastCore) != nodename) uniqueCoreSuccessor[lastCore] = std::numeric_limits<size_t>::max();
	// 		if (uniqueCorePredecessor.count(nodename) == 0) uniqueCorePredecessor[nodename] = lastCore;
	// 		if (uniqueCoreSuccessor.count(lastCore) == 0) uniqueCoreSuccessor[lastCore] = nodename;
	// 		lastCore = nodename;
	// 	}
	// }
	// for (const auto& pair : uniqueCoreSuccessor)
	// {
	// 	if (pair.second != std::numeric_limits<size_t>::max()) continue;
	// 	if (uniqueCorePredecessor.count(pair.first) == 0) continue;
	// 	if (uniqueCorePredecessor.at(pair.first) != std::numeric_limits<size_t>::max()) continue;
	// 	potentialCoreNodes.erase(pair.first);
	// }
	std::unordered_map<size_t, size_t> uniqueCoreIndex;
	for (const auto& read : paths)
	{
		size_t coreIndex = 0;
		for (const auto& node : read.path)
		{
			size_t nodename = node.id();
			if (potentialCoreNodes.count(nodename) == 0) continue;
			if (uniqueCoreIndex.count(nodename) == 1 && uniqueCoreIndex.at(nodename) != coreIndex) uniqueCoreIndex[nodename] = std::numeric_limits<size_t>::max();
			if (uniqueCoreIndex.count(nodename) == 0) uniqueCoreIndex[nodename] = coreIndex;
			coreIndex += 1;
		}
	}
	for (const auto& pair : uniqueCoreIndex)
	{
		if (pair.second != std::numeric_limits<size_t>::max()) continue;
		potentialCoreNodes.erase(pair.first);
	}
	return potentialCoreNodes;
}

size_t getDPRowBacktracePos(const std::string& compareQuerySeq, const std::string& compareRefSeq, const size_t columnNum)
{
	std::vector<std::vector<size_t>> DPmatrix;
	DPmatrix.resize(compareQuerySeq.size()+1);
	for (size_t i = 0; i < compareQuerySeq.size()+1; i++)
	{
		DPmatrix[i].resize(compareRefSeq.size()+1, 0);
	}
	for (size_t i = 0; i < compareRefSeq.size()+1; i++)
	{
		DPmatrix[0][i] = i;
	}
	for (size_t i = 1; i < compareQuerySeq.size()+1; i++)
	{
		DPmatrix[i][0] = i;
		for (size_t j = 1; j < compareRefSeq.size()+1; j++)
		{
			bool match = (compareQuerySeq[i-1] == compareRefSeq[j-1]);
			DPmatrix[i][j] = std::min(DPmatrix[i][j-1]+1, DPmatrix[i-1][j]+1);
			DPmatrix[i][j] = std::min(DPmatrix[i][j], DPmatrix[i-1][j-1] + (match ? 0 : 1));
		}
	}
	size_t i = compareQuerySeq.size();
	size_t j = compareRefSeq.size();
	while (j > columnNum)
	{
		assert(j > 0);
		assert(i > 0);
		bool match = (compareQuerySeq[i-1] == compareRefSeq[j-1]);
		if (DPmatrix[i-1][j-1] + (match ? 0 : 1) == DPmatrix[i][j])
		{
			i -= 1;
			j -= 1;
		}
		else if (DPmatrix[i-1][j] + 1 == DPmatrix[i][j])
		{
			i -= 1;
		}
		else if (DPmatrix[i][j-1] + 1 == DPmatrix[i][j])
		{
			j -= 1;
		}
		else
		{
			assert(false);
		}
	}
	return i;
}

size_t getExactBreakPos(const std::string& nodeseq, const std::string& consensusSeq, const size_t approxPosition)
{
	const size_t flankSize = 10;
	assert(nodeseq.size() >= flankSize);
	assert(consensusSeq.size() >= flankSize);
	std::string compareQuerySeq;
	std::string compareRefSeq;
	size_t midPos = 0;
	if (approxPosition < flankSize)
	{
		compareQuerySeq += nodeseq.substr(0, approxPosition);
		compareRefSeq += consensusSeq.substr(consensusSeq.size() - approxPosition);
		midPos = approxPosition;
	}
	else
	{
		compareQuerySeq += nodeseq.substr(approxPosition-flankSize, flankSize);
		compareRefSeq += consensusSeq.substr(consensusSeq.size() - flankSize);
		midPos = flankSize;
	}
	assert(compareRefSeq.size() == compareQuerySeq.size());
	if (approxPosition+flankSize > nodeseq.size())
	{
		compareQuerySeq += nodeseq.substr(approxPosition, flankSize);
		compareRefSeq += consensusSeq.substr(0, nodeseq.size() - approxPosition);
	}
	else
	{
		compareQuerySeq += nodeseq.substr(approxPosition, flankSize);
		compareRefSeq += consensusSeq.substr(0, flankSize);
	}
	assert(compareRefSeq.size() == compareQuerySeq.size());
	size_t zeroPos = approxPosition;
	if (approxPosition >= flankSize)
	{
		zeroPos = approxPosition - flankSize;
	}
	else
	{
		zeroPos = 0;
	}
	size_t exactMatchPos = zeroPos + getDPRowBacktracePos(compareQuerySeq, compareRefSeq, midPos);
	return exactMatchPos;
}

std::tuple<std::unordered_set<size_t>, std::unordered_map<Node, size_t>, std::unordered_map<Node, size_t>> getBorderNodes(const Path& heavyPath, const GfaGraph& graph)
{
	const size_t k = 101;
	const size_t borderSize = 200;
	std::string consensusSequence = heavyPath.getSequence(graph.nodeSeqs);
	std::unordered_map<std::string, size_t> startKmers;
	std::unordered_map<std::string, size_t> endKmers;
	for (size_t i = 0; i < borderSize; i++)
	{
		std::string subs = consensusSequence.substr(i, k);
		if (endKmers.count(subs) == 1 || startKmers.count(subs) == 1)
		{
			endKmers[subs] = std::numeric_limits<size_t>::max();
			startKmers[subs] = std::numeric_limits<size_t>::max();
		}
		else
		{
			startKmers[subs] = i;
		}
		subs = consensusSequence.substr(consensusSequence.size()-k-i, k);
		if (endKmers.count(subs) == 1 || startKmers.count(subs) == 1)
		{
			endKmers[subs] = std::numeric_limits<size_t>::max();
			startKmers[subs] = std::numeric_limits<size_t>::max();
		}
		else
		{
			endKmers[subs] = i;
		}
	}
	for (size_t i = borderSize; i < consensusSequence.size() - borderSize; i++)
	{
		std::string subs = consensusSequence.substr(i, k);
		if (endKmers.count(subs) == 1 || startKmers.count(subs) == 1)
		{
			endKmers[subs] = std::numeric_limits<size_t>::max();
			startKmers[subs] = std::numeric_limits<size_t>::max();
		}
	}
	std::unordered_set<size_t> borderNodes;
	std::unordered_map<Node, size_t> pathStartClip;
	std::unordered_map<Node, size_t> pathEndClip;
	for (size_t nodeid = 0; nodeid < graph.nodeSeqs.size(); nodeid++)
	{
		std::vector<size_t> fwBreaks;
		for (size_t i = 0; i < graph.nodeSeqs[nodeid].size()-k; i++)
		{
			std::string subs = graph.nodeSeqs[nodeid].substr(i, k);
			if (startKmers.count(subs) == 1)
			{
				size_t pos = startKmers.at(subs);
				if (pos != std::numeric_limits<size_t>::max() && pos <= i)
				{
					fwBreaks.push_back(i - pos);
				}
			}
			if (endKmers.count(subs) == 1)
			{
				size_t pos = endKmers.at(subs);
				if (pos != std::numeric_limits<size_t>::max() && i+pos+k < graph.nodeSeqs[nodeid].size())
				{
					fwBreaks.push_back(i+pos+k);
				}
			}
		}
		std::vector<size_t> bwBreaks;
		std::string reverseSeq = revcomp(graph.nodeSeqs[nodeid]);
		for (size_t i = 0; i < reverseSeq.size()-k; i++)
		{
			std::string subs = reverseSeq.substr(i, k);
			if (startKmers.count(subs) == 1)
			{
				size_t pos = startKmers.at(subs);
				if (pos != std::numeric_limits<size_t>::max() && pos <= i)
				{
					bwBreaks.push_back(i - pos);
				}
			}
			if (endKmers.count(subs) == 1)
			{
				size_t pos = endKmers.at(subs);
				if (pos != std::numeric_limits<size_t>::max() && i+pos+k < reverseSeq.size())
				{
					bwBreaks.push_back(i+pos+k);
				}
			}
		}
		if (fwBreaks.size() >= 3)
		{
			size_t breakPos = getExactBreakPos(graph.nodeSeqs[nodeid], consensusSequence, fwBreaks[fwBreaks.size()/2]);
			// size_t breakPos = fwBreaks[fwBreaks.size()/2];
			pathStartClip[Node { nodeid, true }] = breakPos;
			pathEndClip[Node { nodeid, true }] = graph.nodeSeqs[nodeid].size() - breakPos;
			breakPos = graph.nodeSeqs[nodeid].size() - 1 - breakPos;
			pathStartClip[Node { nodeid, false }] = breakPos;
			pathEndClip[Node { nodeid, false }] = graph.nodeSeqs[nodeid].size() - breakPos;
			borderNodes.insert(nodeid);
		}
		if (bwBreaks.size() >= 3)
		{
			size_t breakPos = getExactBreakPos(revcomp(graph.nodeSeqs[nodeid]), consensusSequence, bwBreaks[fwBreaks.size()/2]);
			// size_t breakPos = bwBreaks[bwBreaks.size()/2];
			pathStartClip[Node { nodeid, false }] = breakPos;
			pathEndClip[Node { nodeid, false }] = graph.nodeSeqs[nodeid].size() - breakPos;
			breakPos = graph.nodeSeqs[nodeid].size() - 1 - breakPos;
			pathStartClip[Node { nodeid, true }] = breakPos;
			pathEndClip[Node { nodeid, true }] = graph.nodeSeqs[nodeid].size() - breakPos;
			borderNodes.insert(nodeid);
		}
	}
	return std::make_tuple(borderNodes, pathStartClip, pathEndClip);
}

std::vector<OntLoop> extractLoopSequences(const std::vector<ReadPath>& correctedPaths, const Path& heavyPath, const size_t minLength, const GfaGraph& graph, const std::unordered_set<size_t>& borderNodes, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
{
	std::vector<OntLoop> result;
	for (const auto& read : correctedPaths)
	{
		size_t lastBreak = std::numeric_limits<size_t>::max();
		std::vector<size_t> approxReadStartPoses;
		std::vector<size_t> approxReadEndPoses;
		size_t pathLength = getPathLength(read.path, graph.nodeSeqs, graph.edges);
		assert(pathLength > read.pathEndClip);
		size_t currentPos = 0;
		for (size_t i = 0; i < read.path.size(); i++)
		{
			approxReadStartPoses.push_back(currentPos);
			approxReadEndPoses.push_back(currentPos + graph.nodeSeqs.at(read.path[i].id()).size());
			currentPos += graph.nodeSeqs.at(read.path[i].id()).size();
			if (i+1 < read.path.size()) currentPos -= getOverlap(read.path[i], read.path[i+1], graph.edges);
		}
		size_t minPos = pathLength;
		size_t maxPos = 0;
		for (size_t i = 0; i < read.path.size(); i++)
		{
			if (approxReadStartPoses[i] < read.pathStartClip) approxReadStartPoses[i] = read.pathStartClip;
			if (approxReadEndPoses[i] < read.pathStartClip) approxReadEndPoses[i] = read.pathStartClip;
			if (approxReadEndPoses[i] > pathLength - read.pathEndClip) approxReadEndPoses[i] = pathLength - read.pathEndClip;
			if (approxReadStartPoses[i] > pathLength - read.pathEndClip) approxReadStartPoses[i] = pathLength - read.pathEndClip;
			assert(approxReadStartPoses[i] <= approxReadEndPoses[i]);
			minPos = std::min(minPos, approxReadStartPoses[i]);
			maxPos = std::max(maxPos, approxReadEndPoses[i]);
		}
		double readScaleMultiplier = (double)(read.readEnd - read.readStart) / (double)(maxPos - minPos);
		for (size_t i = 0; i < read.path.size(); i++)
		{
			approxReadStartPoses[i] -= minPos;
			approxReadStartPoses[i] *= readScaleMultiplier;
			approxReadStartPoses[i] += read.readStart;
			approxReadEndPoses[i] -= minPos;
			approxReadEndPoses[i] *= readScaleMultiplier;
			approxReadEndPoses[i] += read.readStart;
		}
		for (size_t i = 0; i < read.path.size(); i++)
		{
			if (borderNodes.count(read.path[i].id()) == 0) continue;
			if (lastBreak == std::numeric_limits<size_t>::max())
			{
				lastBreak = i;
				continue;
			}
			std::vector<Node> looppath { read.path.begin()+lastBreak, read.path.begin()+i+1 };
			size_t len = getPathLength(looppath, graph.nodeSeqs, graph.edges);
			size_t clip = pathStartClip.at(read.path[lastBreak]) + pathEndClip.at(read.path[i]);
			if (len > clip && len - clip >= minLength)
			{
				result.emplace_back();
				result.back().originalReadLength = read.readLength;
				result.back().path = looppath;
				result.back().readName = nameWithoutTags(read.readName);
				result.back().approxStart = approxReadStartPoses[lastBreak] + pathStartClip.at(read.path[lastBreak]);
				result.back().approxEnd = approxReadEndPoses[i] - pathEndClip.at(read.path[i]);
				if (read.reverse)
				{
					result.back().approxStart = read.readEnd - (result.back().approxStart - read.readStart);
					result.back().approxEnd = read.readEnd - (result.back().approxEnd - read.readStart);
				}
			}
			lastBreak = i;
		}
	}
	return result;
}

size_t getEditDistanceWfa(const PathSequenceView& left, const PathSequenceView& right, const size_t maxEdits)
{
	if (left.size() > right.size()+maxEdits) return maxEdits+1;
	if (right.size() > left.size()+maxEdits) return maxEdits+1;
	std::vector<size_t> offsets;
	std::vector<size_t> nextOffsets;
	offsets.resize(1, 0);
	const size_t leftSize = left.size();
	const size_t rightSize = right.size();
	offsets[0] = getMatchLength(left, right, 0, 0);
	if (offsets[0] == leftSize) return rightSize - leftSize;
	if (offsets[0] == rightSize) return leftSize - rightSize;
	size_t zeroPos = 1;
	for (size_t score = 0; score < maxEdits; score++)
	{
		nextOffsets.resize(offsets.size()+2);
		for (size_t i = 0; i < nextOffsets.size(); i++)
		{
			nextOffsets[i] = 0;
			if (zeroPos > i && zeroPos-i > leftSize) continue;
			if (zeroPos+i >= rightSize) break;
			if (i >= 1 && i < nextOffsets.size()-1) nextOffsets[i] = offsets[i-1]+1;
			if (i >= 2) nextOffsets[i] = std::max(nextOffsets[i], offsets[i-2]);
			if (i < nextOffsets.size()-2) nextOffsets[i] = std::max(nextOffsets[i], offsets[i]+1);
			assert(nextOffsets[i] + i >= zeroPos);
			nextOffsets[i] = std::min(nextOffsets[i], rightSize - i + zeroPos);
			assert(nextOffsets[i]+i-zeroPos <= rightSize);
			nextOffsets[i] = std::min(nextOffsets[i], leftSize);
			assert(nextOffsets[i]+i-zeroPos <= rightSize);
			assert(nextOffsets[i] <= leftSize);
			nextOffsets[i] += getMatchLength(left, right, nextOffsets[i], nextOffsets[i]+i-zeroPos);
			if (nextOffsets[i] == leftSize && nextOffsets[i]+i-zeroPos == rightSize) return score;
		}
		std::swap(offsets, nextOffsets);
		zeroPos += 1;
		assert(score < leftSize+rightSize);
	}
	return maxEdits+1;
}

size_t getEditDistancePossiblyMemoized(const std::vector<Node>& left, const std::vector<Node>& right, const size_t leftStartClipBp, const size_t rightStartClipBp, const size_t leftEndClipBp, const size_t rightEndClipBp, const GfaGraph& graph, const size_t maxEdits, std::unordered_map<std::pair<std::vector<Node>, std::vector<Node>>, size_t>& memoizedEditDistances)
{
	if (left == right) return 0;
	auto key = canon(left, right);
	size_t add;
	if (memoizedEditDistances.count(key) == 0 || (leftStartClipBp != 0 || rightStartClipBp != 0 || leftEndClipBp != 0 || rightEndClipBp != 0))
	{
		auto leftSubseq = getSequenceView(left, graph.nodeSeqsTwobit, graph.revCompNodeSeqsTwobit, graph.edges, leftStartClipBp, leftEndClipBp);
		auto rightSubseq = getSequenceView(right, graph.nodeSeqsTwobit, graph.revCompNodeSeqsTwobit, graph.edges, rightStartClipBp, rightEndClipBp);
		add = getEditDistanceWfa(leftSubseq, rightSubseq, maxEdits);
		assert(add == getEditDistanceWfa(rightSubseq, leftSubseq, maxEdits));
		if (leftStartClipBp == 0 && rightStartClipBp == 0 && leftEndClipBp == 0 && rightEndClipBp == 0) memoizedEditDistances[key] = add;
	}
	else
	{
		add = memoizedEditDistances.at(key);
	}
	return add;
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

std::vector<std::pair<size_t, size_t>> getNodeMatches(const std::vector<Node>& left, size_t leftIndex, size_t leftStart, size_t leftEnd, const std::vector<Node>& right, size_t rightIndex, size_t rightStart, size_t rightEnd, const std::vector<phmap::flat_hash_map<Node, size_t>>& nodeCountIndex, const std::vector<phmap::flat_hash_map<Node, size_t>>& nodePosIndex)
{
	std::vector<std::pair<size_t, size_t>> unfilteredMatches;
	if (leftEnd == leftStart) return unfilteredMatches;
	if (rightEnd == rightStart) return unfilteredMatches;
	assert(leftEnd > leftStart);
	assert(rightEnd > rightStart);
	assert(leftEnd <= left.size());
	assert(rightEnd <= right.size());
	bool allIncreasing = true;
	for (size_t i = leftStart; i < leftEnd; i++)
	{
		auto found = nodeCountIndex[rightIndex].find(left[i]);
		if (found == nodeCountIndex[rightIndex].end()) continue;
		if (found->second != 1) continue;
		found = nodeCountIndex[leftIndex].find(left[i]);
		if (found == nodeCountIndex[leftIndex].end()) continue;
		if (found->second != 1) continue;
		size_t rightPos = nodePosIndex[rightIndex].at(left[i]);
		if (rightPos < rightStart || rightPos >= rightEnd) continue;
		unfilteredMatches.emplace_back(i, nodePosIndex[rightIndex].at(left[i]));
		if (unfilteredMatches.size() >= 2 && unfilteredMatches.back().second <= unfilteredMatches[unfilteredMatches.size()-2].second) allIncreasing = false;
	}
	if (allIncreasing) return unfilteredMatches;
	std::vector<bool> matchInConflict;
	matchInConflict.resize(unfilteredMatches.size(), false);
	for (size_t i = 1; i < unfilteredMatches.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (unfilteredMatches[i].first <= unfilteredMatches[j].first || unfilteredMatches[i].second <= unfilteredMatches[j].second)
			{
				matchInConflict[i] = true;
				matchInConflict[j] = true;
			}
		}
	}
	std::vector<std::pair<size_t, size_t>> result;
	for (size_t i = 0; i < unfilteredMatches.size(); i++)
	{
		if (matchInConflict[i]) continue;
		result.emplace_back(unfilteredMatches[i]);
	}
	assert(result.size()+2 <= unfilteredMatches.size());
	return result;
}

size_t getEditDistance(const std::vector<Node>& left, const size_t leftIndex, const std::vector<Node>& right, const size_t rightIndex, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t maxEdits, const std::unordered_set<size_t>& coreNodes, const std::vector<phmap::flat_hash_map<Node, size_t>>& nodeCountIndex, const std::vector<phmap::flat_hash_map<Node, size_t>>& nodePosIndex, std::unordered_map<std::pair<std::vector<Node>, std::vector<Node>>, size_t>& memoizedEditDistances)
{
	std::vector<size_t> leftCoreMatchPositions;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (coreNodes.count(left[i].id()) == 0) continue;
		leftCoreMatchPositions.push_back(i);
	}
	std::vector<size_t> rightCoreMatchPositions;
	for (size_t i = 0; i < right.size(); i++)
	{
		if (coreNodes.count(right[i].id()) == 0) continue;
		rightCoreMatchPositions.push_back(i);
	}
	assert(leftCoreMatchPositions.size() == coreNodes.size());
	assert(rightCoreMatchPositions.size() == coreNodes.size());
	for (size_t i = 0; i < leftCoreMatchPositions.size(); i++)
	{
		assert(left[leftCoreMatchPositions[i]] == right[rightCoreMatchPositions[i]]);
	}
	std::vector<std::pair<size_t, size_t>> nodeMatches;
	size_t lastLeftStart = 0;
	size_t lastRightStart = 0;
	for (size_t i = 0; i < leftCoreMatchPositions.size(); i++)
	{
		auto addMatches = getNodeMatches(left, leftIndex, lastLeftStart, leftCoreMatchPositions[i], right, rightIndex, lastRightStart, rightCoreMatchPositions[i], nodeCountIndex, nodePosIndex);
		nodeMatches.insert(nodeMatches.end(), addMatches.begin(), addMatches.end());
		nodeMatches.emplace_back(leftCoreMatchPositions[i], rightCoreMatchPositions[i]);
		lastLeftStart = leftCoreMatchPositions[i]+1;
		lastRightStart = rightCoreMatchPositions[i]+1;
	}
	auto addMatches = getNodeMatches(left, leftIndex, lastLeftStart, left.size(), right, rightIndex, lastRightStart, right.size(), nodeCountIndex, nodePosIndex);
	nodeMatches.insert(nodeMatches.end(), addMatches.begin(), addMatches.end());
	assert(nodeMatches.size() >= coreNodes.size());
	if (nodeMatches.size() == 0)
	{
		return getEditDistancePossiblyMemoized(left, right, pathStartClip.at(left[0]), pathStartClip.at(right[0]), pathEndClip.at(left.back()), pathEndClip.at(right.back()), graph, maxEdits, memoizedEditDistances);
	}
	assert(nodeMatches.size() >= 1);
	size_t result = 0;
	size_t add = 0;
	for (size_t i = 1; i < nodeMatches.size(); i++)
	{
		assert(nodeMatches[i].first > nodeMatches[i-1].first);
		assert(nodeMatches[i].second > nodeMatches[i-1].second);
		if (nodeMatches[i].first == nodeMatches[i-1].first+1 && nodeMatches[i].second == nodeMatches[i-1].second+1) continue;
		std::vector<Node> leftPath { left.begin() + nodeMatches[i-1].first, left.begin()+nodeMatches[i].first+1 };
		std::vector<Node> rightPath { right.begin() + nodeMatches[i-1].second, right.begin()+nodeMatches[i].second+1 };
		add = getEditDistancePossiblyMemoized(leftPath, rightPath, 0, 0, 0, 0, graph, maxEdits, memoizedEditDistances);
		result += add;
		if (result >= maxEdits) return maxEdits+1;
	}
	std::vector<Node> leftPath { left.begin(), left.begin()+nodeMatches[0].first+1 };
	std::vector<Node> rightPath { right.begin(), right.begin()+nodeMatches[0].second+1 };
	add = getEditDistancePossiblyMemoized(leftPath, rightPath, pathStartClip.at(leftPath[0]), pathStartClip.at(rightPath[0]), 0, 0, graph, maxEdits, memoizedEditDistances);
	result += add;
	if (result >= maxEdits) return maxEdits+1;
	leftPath = std::vector<Node> { left.begin()+nodeMatches.back().first, left.end() };
	rightPath = std::vector<Node> { right.begin()+nodeMatches.back().second, right.end() };
	add = getEditDistancePossiblyMemoized(leftPath, rightPath, 0, 0, pathEndClip.at(leftPath.back()), pathEndClip.at(rightPath.back()), graph, maxEdits, memoizedEditDistances);
	result += add;
	if (result >= maxEdits) return maxEdits+1;
	return result;
}

size_t countNeedsAligning(const std::vector<size_t>& loopLengths, size_t maxEdits)
{
	size_t result = 0;
	size_t j = 0;
	for (size_t i = 1; i < loopLengths.size(); i++)
	{
		while (loopLengths[j]+maxEdits < loopLengths[i]) j += 1;
		assert(j <= i);
		result += i-j;
	}
	return result;
}

size_t getMatchCount(const std::vector<size_t>& left, const std::vector<size_t>& right)
{
	size_t result = 0;
	size_t leftindex = 0;
	size_t rightindex = 0;
	while (leftindex < left.size() && rightindex < right.size())
	{
		if (left[leftindex] == right[rightindex])
		{
			result += 1;
			leftindex += 1;
			rightindex += 1;
			continue;
		}
		else if (left[leftindex] < right[rightindex])
		{
			leftindex += 1;
		}
		else
		{
			assert(right[rightindex] < left[leftindex]);
			rightindex += 1;
		}
	}
	return result;
}

bool bubbleIsPhased(const std::vector<std::vector<std::vector<size_t>>>& readsPerBubble, const size_t left, const size_t right)
{
	if (readsPerBubble[left].size() != readsPerBubble[right].size()) return false;
	if (readsPerBubble[left].size() < 2) return false;
	for (size_t i = 0; i < readsPerBubble[left].size(); i++)
	{
		if (readsPerBubble[left][i].size() < minPhaseCoverage) return false;
		if (readsPerBubble[right][i].size() < minPhaseCoverage) return false;
	}
	std::vector<bool> rightHasMatch;
	rightHasMatch.resize(readsPerBubble[right].size(), false);
	for (size_t i = 0; i < readsPerBubble[left].size(); i++)
	{
		size_t matchesRight = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < readsPerBubble[right].size(); j++)
		{
			size_t matchCount = getMatchCount(readsPerBubble[left][i], readsPerBubble[right][j]);
			if (matchCount == 0) continue;
			if (matchCount < minPhaseCoverage) return false;
			if (matchesRight != std::numeric_limits<size_t>::max()) return false;
			matchesRight = j;
		}
		if (matchesRight == std::numeric_limits<size_t>::max()) return false;
		if (rightHasMatch[matchesRight]) return false;
		rightHasMatch[matchesRight] = true;
	}
	return true;
}

void phaseAndAppend(const std::vector<OntLoop>& paths, std::vector<std::vector<OntLoop>>& result)
{
	if (paths.size() < minPhaseCoverage)
	{
		result.emplace_back(paths);
		return;
	}
	auto coreNodes = getCoreNodes(paths);
	std::vector<std::map<std::vector<Node>, size_t>> bubbleAlleles;
	bubbleAlleles.resize(coreNodes.size());
	std::vector<std::vector<std::vector<size_t>>> readsPerBubble;
	readsPerBubble.resize(coreNodes.size());
	std::vector<std::vector<size_t>> allelesPerRead;
	for (size_t i = 0; i < paths.size(); i++)
	{
		size_t lastCore = std::numeric_limits<size_t>::max();
		size_t coreNum = 0;
		allelesPerRead.emplace_back();
		for (size_t j = 0; j < paths[i].path.size(); j++)
		{
			if (coreNodes.count(paths[i].path[j].id()) == 0) continue;
			if (lastCore == std::numeric_limits<size_t>::max())
			{
				lastCore = j;
				continue;
			}
			std::vector<Node> subpath { paths[i].path.begin()+lastCore, paths[i].path.begin()+j+1 };
			if (bubbleAlleles[coreNum].count(subpath) == 0)
			{
				size_t count = readsPerBubble[coreNum].size();
				bubbleAlleles[coreNum][subpath] = count;
				readsPerBubble[coreNum].emplace_back();
			}
			readsPerBubble[coreNum][bubbleAlleles[coreNum].at(subpath)].push_back(i);
			allelesPerRead.back().push_back(bubbleAlleles[coreNum].at(subpath));
			coreNum += 1;
		}
	}
	for (size_t i = 0; i < readsPerBubble.size(); i++)
	{
		for (size_t j = 0; j < readsPerBubble[i].size(); j++)
		{
			std::sort(readsPerBubble[i][j].begin(), readsPerBubble[i][j].end());
		}
	}
	std::unordered_set<size_t> phaseInformativeBubblesSet;
	for (size_t i = 0; i < bubbleAlleles.size(); i++)
	{
		for (size_t j = i+1; j < bubbleAlleles.size(); j++)
		{
			if (!bubbleIsPhased(readsPerBubble, i, j)) continue;
			phaseInformativeBubblesSet.insert(i);
			phaseInformativeBubblesSet.insert(j);
		}
	}
	if (phaseInformativeBubblesSet.size() < 2)
	{
		result.emplace_back(paths);
		return;
	}
	std::vector<size_t> phaseInformativeBubbles { phaseInformativeBubblesSet.begin(), phaseInformativeBubblesSet.end() };
	std::sort(phaseInformativeBubbles.begin(), phaseInformativeBubbles.end());
	std::map<std::vector<size_t>, size_t> bubbleToCluster;
	std::vector<std::vector<OntLoop>> phasedClusters;
	for (size_t i = 0; i < paths.size(); i++)
	{
		std::vector<size_t> allelesHere;
		for (size_t index : phaseInformativeBubbles)
		{
			allelesHere.push_back(allelesPerRead[i][index]);
		}
		if (bubbleToCluster.count(allelesHere) == 0)
		{
			size_t count = bubbleToCluster.size();
			bubbleToCluster[allelesHere] = count;
			phasedClusters.emplace_back();
		}
		phasedClusters[bubbleToCluster.at(allelesHere)].emplace_back(paths[i]);
	}
	assert(phasedClusters.size() >= 2);
	for (size_t i = 0; i < phasedClusters.size(); i++)
	{
		phaseAndAppend(phasedClusters[i], result);
	}
}

std::vector<std::vector<OntLoop>> phaseClusters(const std::vector<std::vector<OntLoop>>& unphasedClusters)
{
	std::vector<std::vector<OntLoop>> result;
	for (size_t i = 0; i < unphasedClusters.size(); i++)
	{
		phaseAndAppend(unphasedClusters[i], result);
	}
	return result;
}

std::vector<std::vector<OntLoop>> roughClusterLoopSequences(const std::vector<OntLoop>& loops, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::unordered_set<size_t>& coreNodes, const size_t maxEdits)
{
	std::vector<size_t> parent;
	parent.resize(loops.size());
	for (size_t i = 0; i < loops.size(); i++)
	{
		parent[i] = i;
	}
	std::vector<size_t> loopLengths;
	loopLengths.reserve(loops.size());
	for (size_t i = 0; i < loops.size(); i++)
	{
		size_t len = getPathLength(loops[i].path, graph.nodeSeqs, graph.edges);
		len -= pathStartClip.at(loops[i].path[0]);
		len -= pathEndClip.at(loops[i].path.back());
		loopLengths.emplace_back(len);
	}
	std::vector<phmap::flat_hash_map<Node, size_t>> nodeCountIndex;
	std::vector<phmap::flat_hash_map<Node, size_t>> nodePosIndex;
	nodeCountIndex.resize(loops.size());
	nodePosIndex.resize(loops.size());
	for (size_t i = 0; i < loops.size(); i++)
	{
		for (size_t j = 0; j < loops[i].path.size(); j++)
		{
			nodeCountIndex[i][loops[i].path[j]] += 1;
			nodePosIndex[i][loops[i].path[j]] = j;
		}
	}
	std::unordered_map<std::pair<std::vector<Node>, std::vector<Node>>, size_t> memoizedEditDistances;
	size_t sumAligned = 0;
	size_t needsAligning = countNeedsAligning(loopLengths, maxEdits);
	for (size_t i = 0; i < loops.size(); i++)
	{
		while (parent[i] != parent[parent[i]]) parent[i] = parent[parent[i]];
		for (size_t j = i-1; j < i; j--)
		{
			if (sumAligned % 1000000 == 0) std::cerr << "aligning morph path pair " << sumAligned << " / " << needsAligning << std::endl;
			sumAligned += 1;
			assert(loopLengths[j] <= loopLengths[i]);
			if (loopLengths[j]+maxEdits < loopLengths[i]) break;
			while (parent[j] != parent[parent[j]]) parent[j] = parent[parent[j]];
			if (parent[i] == parent[j]) continue;
			size_t edits = getEditDistance(loops[i].path, i, loops[j].path, j, graph, pathStartClip, pathEndClip, maxEdits, coreNodes, nodeCountIndex, nodePosIndex, memoizedEditDistances);
			assert(edits == getEditDistance(loops[j].path, j, loops[i].path, i, graph, pathStartClip, pathEndClip, maxEdits, coreNodes, nodeCountIndex, nodePosIndex, memoizedEditDistances));
			if (edits > maxEdits) continue;
			parent[parent[j]] = parent[i];
		}
	}
	std::vector<size_t> parentToCluster;
	parentToCluster.resize(parent.size(), std::numeric_limits<size_t>::max());
	std::vector<std::vector<OntLoop>> result;
	for (size_t i = 0; i < parent.size(); i++)
	{
		while (parent[i] != parent[parent[i]]) parent[i] = parent[parent[i]];
		if (parent[i] != i) continue;
		parentToCluster[i] = result.size();
		result.emplace_back();
	}
	for (size_t i = 0; i < parent.size(); i++)
	{
		result[parentToCluster[parent[i]]].emplace_back(loops[i]);
	}
	return result;
}

std::vector<std::vector<OntLoop>> clusterByDbscan(const std::vector<OntLoop>& cluster, const size_t epsilon, const size_t minPoints, const std::vector<std::vector<uint16_t>>& editDistanceMatrix)
{
	// O(n^2) is fast enough for the cluster sizes here
	std::vector<bool> corePoint;
	corePoint.resize(cluster.size(), false);
	for (size_t i = 0; i < cluster.size(); i++)
	{
		size_t nearPoints = 0;
		for (size_t j = 0; j < cluster.size(); j++)
		{
			if (i == j) continue;
			size_t canoni = std::max(i, j);
			size_t canonj = std::min(i, j);
			if (editDistanceMatrix[canoni][canonj] > epsilon) continue;
			nearPoints += 1;
			if (nearPoints == minPoints)
			{
				corePoint[i] = true;
				break;
			}
		}
	}
	std::vector<size_t> parent;
	parent.resize(cluster.size(), 0);
	for (size_t i = 0; i < cluster.size(); i++)
	{
		parent[i] = i;
	}
	for (size_t i = 0; i < cluster.size(); i++)
	{
		if (!corePoint[i]) continue;
		for (size_t j = 0; j < cluster.size(); j++)
		{
			if (i == j) continue;
			if (!corePoint[j]) continue;
			size_t canoni = std::max(i, j);
			size_t canonj = std::min(i, j);
			if (editDistanceMatrix[canoni][canonj] > epsilon) continue;
			while (parent[parent[i]] != parent[i]) parent[i] = parent[parent[i]];
			while (parent[parent[j]] != parent[j]) parent[j] = parent[parent[j]];
			parent[parent[j]] = parent[i];
		}
	}
	for (size_t j = 0; j < parent.size(); j++)
	{
		while (parent[j] != parent[parent[j]]) parent[j] = parent[parent[j]];
	}
	std::vector<size_t> uniqueClusterForEdgePoint;
	uniqueClusterForEdgePoint.resize(cluster.size(), std::numeric_limits<size_t>::max()-1);
	for (size_t i = 0; i < cluster.size(); i++)
	{
		if (corePoint[i]) continue;
		assert(parent[i] == i);
		for (size_t j = 0; j < cluster.size(); j++)
		{
			if (i == j) continue;
			if (!corePoint[j]) continue;
			size_t canoni = std::max(i, j);
			size_t canonj = std::min(i, j);
			if (editDistanceMatrix[canoni][canonj] > epsilon) continue;
			if (uniqueClusterForEdgePoint[i] == std::numeric_limits<size_t>::max()-1)
			{
				assert(parent[j] == parent[parent[j]]);
				uniqueClusterForEdgePoint[i] = parent[j];
				parent[i] = parent[j];
			}
			else if (uniqueClusterForEdgePoint[i] != parent[j])
			{
				uniqueClusterForEdgePoint[i] = std::numeric_limits<size_t>::max();
			}
		}
	}
	std::vector<size_t> clusterMapping;
	std::vector<std::vector<OntLoop>> result;
	clusterMapping.resize(cluster.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < cluster.size(); i++)
	{
		if (!corePoint[i] && (uniqueClusterForEdgePoint[i] == std::numeric_limits<size_t>::max() || uniqueClusterForEdgePoint[i] == std::numeric_limits<size_t>::max()-1)) continue;
		if (clusterMapping[parent[i]] == std::numeric_limits<size_t>::max())
		{
			clusterMapping[parent[i]] = result.size();
			result.emplace_back();
		}
		result[clusterMapping[parent[i]]].emplace_back(cluster[i]);
	}
	return result;
}

std::vector<std::vector<OntLoop>> densityClusterLoops(const std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::unordered_set<size_t>& coreNodes, const size_t maxEdits, const size_t minPoints, const size_t minEdits)
{
	std::vector<size_t> editHistogram;
	editHistogram.resize(maxEdits, 0);
	std::vector<std::vector<std::vector<uint16_t>>> editDistanceMatrices;
	editDistanceMatrices.resize(clusters.size());
	for (size_t clusteri = 0; clusteri < clusters.size(); clusteri++)
	{
		editDistanceMatrices[clusteri].resize(clusters[clusteri].size());
		std::vector<phmap::flat_hash_map<Node, size_t>> nodeCountIndex;
		std::vector<phmap::flat_hash_map<Node, size_t>> nodePosIndex;
		nodeCountIndex.resize(clusters[clusteri].size());
		nodePosIndex.resize(clusters[clusteri].size());
		for (size_t i = 0; i < clusters[clusteri].size(); i++)
		{
			for (size_t j = 0; j < clusters[clusteri][i].path.size(); j++)
			{
				nodeCountIndex[i][clusters[clusteri][i].path[j]] += 1;
				nodePosIndex[i][clusters[clusteri][i].path[j]] = j;
			}
		}
		std::unordered_map<std::pair<std::vector<Node>, std::vector<Node>>, size_t> memoizedEditDistances;
		for (size_t i = 0; i < clusters[clusteri].size(); i++)
		{
			editDistanceMatrices[clusteri][i].resize(i);
			for (size_t j = 0; j < i; j++)
			{
				size_t edits = getEditDistance(clusters[clusteri][i].path, i, clusters[clusteri][j].path, j, graph, pathStartClip, pathEndClip, maxEdits, coreNodes, nodeCountIndex, nodePosIndex, memoizedEditDistances);
				assert(edits == getEditDistance(clusters[clusteri][j].path, j, clusters[clusteri][i].path, i, graph, pathStartClip, pathEndClip, maxEdits, coreNodes, nodeCountIndex, nodePosIndex, memoizedEditDistances));
				if (edits >= maxEdits) edits = maxEdits-1;
				editDistanceMatrices[clusteri][i][j] = edits;
				editHistogram[edits] += 1;
			}
		}
	}
	size_t histogramPeak = 0;
	for (size_t i = 0; i < editHistogram.size()-1; i++)
	{
		if (editHistogram[i] >= editHistogram[histogramPeak]) histogramPeak = i;
	}
	std::cerr << "edit distance peak at " << histogramPeak << std::endl;
	size_t newEditDistance = histogramPeak;
	if (newEditDistance >= maxEdits)
	{
		return clusters;
	}
	if (newEditDistance < minEdits) newEditDistance = minEdits;
	std::cerr << "recluster with max edit distance " << newEditDistance << ", min points " << minPoints << std::endl;
	std::vector<std::vector<OntLoop>> result;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		auto partialResult = clusterByDbscan(clusters[i], newEditDistance, minPoints, editDistanceMatrices[i]);
		std::cerr << "cluster " << i << " with " << clusters[i].size() << " reads reclustered to " << partialResult.size() << " clusters, sizes:";
		for (size_t j = 0; j < partialResult.size(); j++)
		{
			std::cerr << " " << partialResult[j].size();
		}
		std::cerr << std::endl;
		while (partialResult.size() > 0)
		{
			result.emplace_back();
			std::swap(result.back(), partialResult.back());
			partialResult.pop_back();
		}
	}
	return result;
}

std::unordered_set<size_t> getMajorityNodes(const std::vector<OntLoop>& paths)
{
	std::unordered_map<size_t, size_t> nodeCoverage;
	std::unordered_set<size_t> notCore;
	std::unordered_set<size_t> firstNodes;
	std::unordered_set<size_t> lastNodes;
	for (const auto& path : paths)
	{
		firstNodes.insert(path.path[0].id());
		lastNodes.insert(path.path.back().id());
		std::unordered_map<size_t, size_t> coverageHere;
		for (auto node : path.path)
		{
			coverageHere[node.id()] += 1;
		}
		for (auto pair : coverageHere)
		{
			if (pair.second != 1) notCore.insert(pair.first);
			nodeCoverage[pair.first] += 1;
		}
	}
	std::unordered_set<size_t> result;
	for (auto pair : nodeCoverage)
	{
		if (notCore.count(pair.first) == 1) continue;
		if (firstNodes.count(pair.first) == 1 && lastNodes.count(pair.first) == 1) continue;
		if (pair.second > paths.size()/2) result.insert(pair.first);
	}
	return result;
}

std::vector<size_t> getMajorityPath(const std::unordered_set<size_t>& coreNodes, const std::vector<OntLoop>& rawPaths)
{
	std::unordered_map<size_t, size_t> firstCount;
	std::unordered_map<size_t, size_t> coverage;
	std::unordered_map<size_t, std::unordered_map<size_t, size_t>> edges;
	for (const auto& path : rawPaths)
	{
		size_t lastCore = std::numeric_limits<size_t>::max();
		for (auto node : path.path)
		{
			if (coreNodes.count(node.id()) == 0) continue;
			coverage[node.id()] += 1;
			if (lastCore == std::numeric_limits<size_t>::max())
			{
				firstCount[node.id()] += 1;
			}
			else
			{
				edges[lastCore][node.id()] += 1;
			}
			lastCore = node.id();
		}
	}
	std::vector<size_t> result;
	size_t startNode = std::numeric_limits<size_t>::max();
	for (auto pair : firstCount)
	{
		if (startNode == std::numeric_limits<size_t>::max() || pair.second > firstCount.at(startNode))
		{
			startNode = pair.first;
		}
	}
	result.emplace_back(startNode);
	size_t pos = startNode;
	while (true)
	{
		if (edges.count(pos) == 0) break;
		size_t maxNode = std::numeric_limits<size_t>::max();
		for (auto pair : edges.at(pos))
		{
			if (maxNode == std::numeric_limits<size_t>::max() || pair.second > edges.at(pos).at(maxNode))
			{
				maxNode = pair.first;
			}
		}
		assert(maxNode != std::numeric_limits<size_t>::max());
		pos = maxNode;
		result.emplace_back(pos);
	}
	return result;
}

std::vector<Node> getConsensusPath(const std::vector<OntLoop>& rawPaths, const GfaGraph& graph)
{
	assert(rawPaths.size() >= 1);
	auto coreNodes = getCoreNodes(rawPaths);
	if (coreNodes.size() == 0)
	{
		coreNodes = getMajorityNodes(rawPaths);
	}
	assert(coreNodes.size() >= 1);
	std::unordered_map<std::vector<Node>, size_t> alleleCounts;
	Node fakeStart { graph.numNodes()+1, true };
	Node fakeEnd { graph.numNodes()+2, true };
	for (const auto& read : rawPaths)
	{
		const auto& path = read.path;
		size_t lastCore = std::numeric_limits<size_t>::max();
		for (size_t i = 0; i < path.size(); i++)
		{
			assert(path[i].id() != fakeStart.id());
			assert(path[i].id() != fakeEnd.id());
			if (coreNodes.count(path[i].id()) == 0) continue;
			std::vector<Node> subpath { path.begin() + (lastCore == std::numeric_limits<size_t>::max() ? 0 : lastCore), path.begin() + i + 1 };
			assert(subpath.size() >= 2 || (i == 0 && subpath.size() == 1));
			if (lastCore == std::numeric_limits<size_t>::max()) subpath.insert(subpath.begin(), fakeStart);
			alleleCounts[subpath] += 1;
			lastCore = i;
		}
		std::vector<Node> subpath { path.begin() + (lastCore == std::numeric_limits<size_t>::max() ? 0 : lastCore), path.end()};
		assert(subpath.size() >= 2 || (lastCore == path.size()-1 && subpath.size() == 1));
		if (lastCore == std::numeric_limits<size_t>::max()) subpath.insert(subpath.begin(), fakeStart);
		subpath.insert(subpath.end(), fakeEnd);
		alleleCounts[subpath] += 1;
	}
	auto corePath = getMajorityPath(coreNodes, rawPaths);
	std::unordered_map<size_t, size_t> coreNodePositionInPath;
	for (size_t i = 0; i < corePath.size(); i++)
	{
		coreNodePositionInPath[corePath[i]] = i+1;
	}
	coreNodePositionInPath[fakeStart.id()] = 0;
	coreNodePositionInPath[fakeEnd.id()] = corePath.size()+1;
	std::vector<std::pair<std::vector<Node>, size_t>> bestAlleles;
	bestAlleles.resize(corePath.size()+1, std::make_pair(std::vector<Node>{}, 0));
	for (auto pair : alleleCounts)
	{
		if (coreNodePositionInPath.count(pair.first[0].id()) == 0) continue;
		if (coreNodePositionInPath.count(pair.first.back().id()) == 0) continue;
		if (coreNodePositionInPath.at(pair.first.back().id()) != coreNodePositionInPath.at(pair.first[0].id()) + 1) continue;
		size_t index = coreNodePositionInPath.at(pair.first[0].id());
		if (pair.second > bestAlleles[index].second) bestAlleles[index] = pair;
	}
	std::vector<Node> consensusPath;
	for (size_t i = 0; i < bestAlleles.size(); i++)
	{
		assert(bestAlleles[i].second >= 1);
		assert(bestAlleles[i].first.size() >= 1);
		if (i == 0)
		{
			consensusPath.insert(consensusPath.end(), bestAlleles[i].first.begin(), bestAlleles[i].first.end());
		}
		else
		{
			assert(consensusPath.size() >= 1);
			assert(consensusPath.back() == bestAlleles[i].first[0]);
			consensusPath.insert(consensusPath.end(), bestAlleles[i].first.begin()+1, bestAlleles[i].first.end());
		}
	}
	assert(consensusPath.size() >= 3);
	assert(consensusPath[0] == fakeStart);
	assert(consensusPath.back() == fakeEnd);
	consensusPath.erase(consensusPath.begin());
	consensusPath.erase(consensusPath.begin()+consensusPath.size()-1);
	assert(consensusPath.size() >= 1);
	return consensusPath;
}

std::string getConsensusSequence(const std::vector<Node>& consensusPath, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
{
	std::string consensusSeq = getSequence(consensusPath, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
	consensusSeq = consensusSeq.substr(pathStartClip.at(consensusPath[0]), consensusSeq.size() - pathStartClip.at(consensusPath[0]) - pathEndClip.at(consensusPath.back()));
	return consensusSeq;
}

std::vector<std::tuple<size_t, size_t, std::string>> getEdits(const std::string& refSequence, const std::string& querySequence, size_t maxEdits)
{
	// reference is horizontal, first bp left, last bp right (one bp one column)
	// query is vertical, first bp top, last bp bottom (one bp one row)
	// edits moves diagonally down-right
	// j moves horizontally, plus right minus left
	// insertion is extra sequence on query (vertical backtrace)
	// deletion is missing sequence on query (horizontal backtrace)
	// wfa matrix pair coord first is edits (row), second is j (column)
	// wfa matrix is triangular, missing parts are implied: edit 0 j 0 is above edit 1 j 1 and up-right from edit 1 j 0
	// real matrix pair coord first is ref pos (column), second is query pos (row)
	if (refSequence == querySequence) return std::vector<std::tuple<size_t, size_t, std::string>>{};
	const std::pair<size_t, size_t> uninitialized { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
	const size_t mismatchCost = 2;
	const size_t gapOpenCost = 2; // only open cost, and 1bp indel costs gapOpenCost+gapExtendCost
	const size_t gapExtendCost = 1;
	// match cost 0
	// gap close cost 0
	Wfa::WfaMatrix matchMatrix;
	Wfa::WfaMatrix insertionMatrix;
	Wfa::WfaMatrix deletionMatrix;
	matchMatrix.resize(1);
	insertionMatrix.resize(1);
	deletionMatrix.resize(1);
	insertionMatrix[0].resize(1, uninitialized);
	deletionMatrix[0].resize(1, uninitialized);
	matchMatrix[0].resize(1);
	matchMatrix[0][0].first = 0;
	matchMatrix[0][0].second = 0;
	while (matchMatrix[0][0].second < refSequence.size() && matchMatrix[0][0].second < querySequence.size() && refSequence[matchMatrix[0][0].second] == querySequence[matchMatrix[0][0].second]) matchMatrix[0][0].second += 1;
	size_t edits = 0;
	std::pair<size_t, size_t> backtraceWfaPosition = uninitialized;
	for (edits = 1; edits < maxEdits; edits++)
	{
		matchMatrix.emplace_back();
		insertionMatrix.emplace_back();
		deletionMatrix.emplace_back();
		matchMatrix.back().resize(2*edits+1, uninitialized);
		insertionMatrix.back().resize(2*edits+1, uninitialized);
		deletionMatrix.back().resize(2*edits+1, uninitialized);
		for (size_t j = 0; j < 2*edits+1; j++)
		{
			//mismatch
			Wfa::updateMatrix(std::make_pair(edits, j), mismatchCost, 1, 1, matchMatrix, matchMatrix, refSequence.size(), querySequence.size());
			// gap open
			Wfa::updateMatrix(std::make_pair(edits, j), gapOpenCost, 0, 0, matchMatrix, insertionMatrix, refSequence.size(), querySequence.size());
			Wfa::updateMatrix(std::make_pair(edits, j), gapOpenCost, 0, 0, matchMatrix, deletionMatrix, refSequence.size(), querySequence.size());
			// gap extend
			Wfa::updateMatrix(std::make_pair(edits, j), gapExtendCost, 1, 0, insertionMatrix, insertionMatrix, refSequence.size(), querySequence.size());
			Wfa::updateMatrix(std::make_pair(edits, j), gapExtendCost, 0, 1, deletionMatrix, deletionMatrix, refSequence.size(), querySequence.size());
			insertionMatrix[edits][j].second = insertionMatrix[edits][j].first;
			deletionMatrix[edits][j].second = deletionMatrix[edits][j].first;
			// gap close
			Wfa::updateMatrix(std::make_pair(edits, j), 0, 0, 0, insertionMatrix, matchMatrix, refSequence.size(), querySequence.size());
			Wfa::updateMatrix(std::make_pair(edits, j), 0, 0, 0, deletionMatrix, matchMatrix, refSequence.size(), querySequence.size());
			//matches
			matchMatrix[edits][j].second = matchMatrix[edits][j].first;
			size_t offset = 0;
			if (matchMatrix[edits][j] == uninitialized) continue;
			std::pair<size_t, size_t> realPos = Wfa::getRealMatrixPosition(std::make_pair(edits, j), matchMatrix[edits][j].first);
			assert(realPos.first <= refSequence.size() && realPos.second <= querySequence.size());
			while (realPos.first+offset < refSequence.size() && realPos.second+offset < querySequence.size() && refSequence[realPos.first+offset] == querySequence[realPos.second+offset]) offset += 1;
			matchMatrix[edits][j].second += offset;
			if (realPos.first+offset == refSequence.size() && realPos.second+offset == querySequence.size())
			{
				backtraceWfaPosition.first = edits;
				backtraceWfaPosition.second = j;
				break;
			}
		}
		if (backtraceWfaPosition != uninitialized) break;
	}
	if (backtraceWfaPosition == uninitialized) return std::vector<std::tuple<size_t, size_t, std::string>> {};
	assert(backtraceWfaPosition != uninitialized);
	std::vector<std::tuple<size_t, size_t, std::string>> result;
	// 0 match 1 insertion 2 deletion
	size_t currentMatrix = 0;
	std::pair<size_t, size_t> lastMatchRealPosition = uninitialized;
	while (backtraceWfaPosition.first != 0)
	{
		switch(currentMatrix)
		{
			case 0:
			{
				assert(backtraceWfaPosition.first < matchMatrix.size());
				assert(backtraceWfaPosition.second < matchMatrix[backtraceWfaPosition.first].size());
				assert(matchMatrix[backtraceWfaPosition.first][backtraceWfaPosition.second] != uninitialized);
				std::pair<size_t, size_t> newMatchRealPosition = Wfa::getRealMatrixPosition(backtraceWfaPosition, matchMatrix[backtraceWfaPosition.first][backtraceWfaPosition.second].second);
				if (lastMatchRealPosition != uninitialized)
				{
					assert(newMatchRealPosition.first <= lastMatchRealPosition.first);
					assert(newMatchRealPosition.second <= lastMatchRealPosition.second);
					result.emplace_back(newMatchRealPosition.first, lastMatchRealPosition.first, querySequence.substr(newMatchRealPosition.second, lastMatchRealPosition.second - newMatchRealPosition.second));
				}
				lastMatchRealPosition = Wfa::getRealMatrixPosition(backtraceWfaPosition, matchMatrix[backtraceWfaPosition.first][backtraceWfaPosition.second].first);
				if (Wfa::canBacktrace(backtraceWfaPosition, mismatchCost, 1, 1, matchMatrix, matchMatrix, refSequence.size(), querySequence.size()))
				{
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, mismatchCost, 1, 1, matchMatrix, matchMatrix);
					break;
				}
				if (Wfa::canBacktrace(backtraceWfaPosition, 0, 0, 0, insertionMatrix, matchMatrix, refSequence.size(), querySequence.size()))
				{
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, 0, 0, 0, matchMatrix, insertionMatrix);
					currentMatrix = 1;
					break;
				}
				else
				{
					assert(Wfa::canBacktrace(backtraceWfaPosition, 0, 0, 0, deletionMatrix, matchMatrix, refSequence.size(), querySequence.size()));
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, 0, 0, 0, matchMatrix, deletionMatrix);
					currentMatrix = 2;
					break;
				}
			}
			case 1:
			{
				assert(backtraceWfaPosition.first < insertionMatrix.size());
				assert(backtraceWfaPosition.second < insertionMatrix[backtraceWfaPosition.first].size());
				assert(insertionMatrix[backtraceWfaPosition.first][backtraceWfaPosition.second] != uninitialized);
				if (Wfa::canBacktrace(backtraceWfaPosition, gapExtendCost, 1, 0, insertionMatrix, insertionMatrix, refSequence.size(), querySequence.size()))
				{
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, gapExtendCost, 1, 0, insertionMatrix, insertionMatrix);
					break;
				}
				else
				{
					assert(Wfa::canBacktrace(backtraceWfaPosition, gapOpenCost, 0, 0, matchMatrix, insertionMatrix, refSequence.size(), querySequence.size()));
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, gapOpenCost, 0, 0, insertionMatrix, matchMatrix);
					currentMatrix = 0;
					break;
				}
			}
			case 2:
			{
				assert(backtraceWfaPosition.first < deletionMatrix.size());
				assert(backtraceWfaPosition.second < deletionMatrix[backtraceWfaPosition.first].size());
				assert(deletionMatrix[backtraceWfaPosition.first][backtraceWfaPosition.second] != uninitialized);
				if (Wfa::canBacktrace(backtraceWfaPosition, gapExtendCost, 0, 1, deletionMatrix, deletionMatrix, refSequence.size(), querySequence.size()))
				{
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, gapExtendCost, 0, 1, deletionMatrix, deletionMatrix);
					break;
				}
				else
				{
					assert(Wfa::canBacktrace(backtraceWfaPosition, gapOpenCost, 0, 0, matchMatrix, deletionMatrix, refSequence.size(), querySequence.size()));
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, gapOpenCost, 0, 0, deletionMatrix, matchMatrix);
					currentMatrix = 0;
					break;
				}
			}
		}
	}
	return result;
}

std::string getPolishedSequence(const std::string& rawSequence, const std::vector<OntLoop>& rawPaths, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
{
	std::map<std::tuple<size_t, size_t, std::string>, size_t> editCounts;
	for (const auto& path : rawPaths)
	{
		auto sequence = getSequence(path.path, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
		sequence = sequence.substr(pathStartClip.at(path.path[0]), sequence.size() - pathStartClip.at(path.path[0]) - pathEndClip.at(path.path.back()));
		auto edits = getEdits(rawSequence, sequence, 2000);
		for (const auto& edit : edits)
		{
			assert(std::get<1>(edit) >= std::get<0>(edit));
			assert(std::get<1>(edit) < rawSequence.size());
			editCounts[edit] += 1;
		}
	}
	std::vector<std::tuple<size_t, size_t, std::string>> pickedEdits;
	for (const auto& pair : editCounts)
	{
		if (!(pair.second > rawPaths.size() * 0.6)) continue; // roundabout comparison to make sure rounding issues are handled more conservatively
		pickedEdits.emplace_back(pair.first);
	}
	std::sort(pickedEdits.begin(), pickedEdits.end(), [](auto& left, auto& right) { return std::get<0>(left) < std::get<0>(right); });
	std::string result;
	size_t lastMatch = 0;
	for (size_t i = 0; i < pickedEdits.size(); i++)
	{
		assert(i == 0 || std::get<0>(pickedEdits[i]) >= std::get<1>(pickedEdits[i-1]));
		assert(std::get<0>(pickedEdits[i]) >= lastMatch);
		result += rawSequence.substr(lastMatch, std::get<0>(pickedEdits[i]) - lastMatch);
		result += std::get<2>(pickedEdits[i]);
		lastMatch = std::get<1>(pickedEdits[i]);
	}
	result += rawSequence.substr(lastMatch);
	return result;
}

std::vector<MorphConsensus> getMorphConsensuses(const std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::string& namePrefix)
{
	std::vector<MorphConsensus> result;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		result.emplace_back();
		result.back().path = getConsensusPath(clusters[i], graph);
		result.back().sequence = getConsensusSequence(result.back().path, graph, pathStartClip, pathEndClip);
		// result.back().sequence = getPolishedSequence(result.back().sequence, clusters[i], graph, pathStartClip, pathEndClip);
		result.back().coverage = clusters[i].size();
		result.back().ontLoops = clusters[i];
		result.back().name = "morphconsensus" + std::to_string(i) + "_coverage" + std::to_string(result.back().coverage);
		if (namePrefix != "") result.back().name = namePrefix + "_" + result.back().name;
	}
	return result;
}

void writeMorphConsensuses(std::string outFile, const std::vector<MorphConsensus>& morphConsensuses)
{
	std::ofstream file { outFile };
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		file << ">" << morphConsensuses[i].name << std::endl;
		file << morphConsensuses[i].sequence << std::endl;
	}
}

void orderLoopsByLength(std::vector<OntLoop>& loops, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
{
	std::vector<std::pair<size_t, size_t>> indexAndLength;
	indexAndLength.reserve(loops.size());
	for (size_t i = 0; i < loops.size(); i++)
	{
		size_t len = getPathLength(loops[i].path, graph.nodeSeqs, graph.edges);
		size_t leftclip = pathStartClip.at(loops[i].path[0]);
		size_t rightclip = pathEndClip.at(loops[i].path.back());
		assert(len > leftclip+rightclip);
		indexAndLength.emplace_back(i, len-leftclip-rightclip);
	}
	std::sort(indexAndLength.begin(), indexAndLength.end(), [](auto left, auto right) { return left.second < right.second; });
	std::vector<OntLoop> result;
	result.resize(loops.size());
	for (size_t i = 0; i < loops.size(); i++)
	{
		std::swap(result[i], loops[indexAndLength[i].first]);
	}
	std::swap(result, loops);
}

void HandleCluster(const ClusterParams& params)
{
	std::cerr << "running MBG" << std::endl;
	runMBG(params.basePath, params.hifiReadPath, params.MBGPath, params.k, params.maxResolveLength, params.numThreads);
	std::cerr << "reading graph" << std::endl;
	GfaGraph graph;
	graph.loadFromFile(params.basePath + "/graph.gfa");
	std::cerr << "getting consensus" << std::endl;
	Path heavyPath = getHeavyPath(graph, params.maxResolveLength);
	std::cerr << "consensus length " << heavyPath.getSequence(graph.nodeSeqs).size() << "bp" << std::endl;
	if (params.orientReferencePath.size() > 0)
	{
		std::cerr << "orienting consensus" << std::endl;
		heavyPath = orientPath(graph, heavyPath, params.orientReferencePath, 101);
	}
	std::cerr << "writing consensus" << std::endl;
	writePathSequence(heavyPath, graph, params.basePath + "/consensus.fa", params.namePrefix);
	writePathGaf(heavyPath, graph, params.basePath + "/consensus_path.gaf");
	std::cerr << "reading read paths" << std::endl;
	std::vector<ReadPath> readPaths = loadReadPaths(params.basePath + "/paths.gaf", graph);
	std::cerr << "getting variants" << std::endl;
	std::vector<Variant> variants = getVariants(graph, heavyPath, readPaths, 3);
	nameVariants(variants, graph, heavyPath);
	std::cerr << variants.size() << " variants" << std::endl;
	std::cerr << "writing variants" << std::endl;
	writeVariants(heavyPath, graph, variants, params.basePath + "/variants.txt");
	std::cerr << "writing variant graph" << std::endl;
	writeVariantGraph(params.basePath + "/graph.gfa", graph, heavyPath, variants, params.basePath + "/variant-graph.gfa");
	std::cerr << "writing allele graph" << std::endl;
	writeAlleleGraph(graph, heavyPath, variants, params.basePath + "/allele-graph.gfa");
	std::cerr << "writing variant vcf" << std::endl;
	writeVariantVCF(params.basePath + "/variants.vcf", heavyPath, graph, variants);
	if (params.annotationFasta.size() > 0)
	{
		std::cerr << "lifting over annotations to consensus" << std::endl;
		liftoverAnnotationsToConsensus(params.basePath, params.basePath + "/consensus.fa", params.annotationFasta, params.annotationGff3);
	}
}

void writeMorphPaths(const std::string& outputFile, const std::vector<MorphConsensus>& morphConsensuses, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
{
	std::ofstream file { outputFile };
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		file << morphConsensuses[i].name << "\t" << morphConsensuses[i].sequence.size() << "\t" << 0 << "\t" << morphConsensuses[i].sequence.size() << "\t+\t";
		for (auto node : morphConsensuses[i].path)
		{
			file << (node.forward() ? ">" : "<") << graph.nodeNames.at(node.id());
		}
		size_t size = getSequence(morphConsensuses[i].path, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges).size();
		file << "\t" << size << "\t" << pathStartClip.at(morphConsensuses[i].path[0]) << "\t" << (size - pathEndClip.at(morphConsensuses[i].path.back())) << "\t" << morphConsensuses[i].sequence.size() << "\t" << morphConsensuses[i].sequence.size() << "\t60" << std::endl;
	}
}

void writeMorphGraphAndReadPaths(const std::string& graphFile, const std::string& pathsFile, const std::vector<MorphConsensus>& morphConsensuses)
{
	std::unordered_map<std::string, size_t> originalReadLength;
	std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t, size_t>>> readFwMatches;
	std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t, size_t>>> readBwMatches;
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		for (const auto& loop : morphConsensuses[i].ontLoops)
		{
			originalReadLength[loop.readName] = loop.originalReadLength;
			if (loop.approxEnd < loop.approxStart)
			{
				readBwMatches[loop.readName].emplace_back(loop.approxStart, loop.approxEnd, i);
			}
			else
			{
				readFwMatches[loop.readName].emplace_back(loop.approxStart, loop.approxEnd, i);
			}
		}
	}
	std::vector<ReadPath> readPaths;
	std::vector<size_t> pathLength;
	for (auto& pair : readFwMatches)
	{
		assert(pair.second.size() >= 1);
		std::sort(pair.second.begin(), pair.second.end(), [](auto left, auto right){ return std::get<0>(left) < std::get<0>(right); });
		for (size_t i = 0; i < pair.second.size(); i++)
		{
			if (i == 0 || std::get<0>(pair.second[i]) + 1000 < std::get<1>(pair.second[i-1]) || std::get<0>(pair.second[i]) > std::get<1>(pair.second[i-1]) + 1000)
			{
				readPaths.emplace_back();
				pathLength.emplace_back(0);
				readPaths.back().readName = pair.first;
				readPaths.back().pathStartClip = 0;
				readPaths.back().pathEndClip = 0;
				readPaths.back().reverse = false;
				readPaths.back().readStart = std::get<0>(pair.second[i]);
				readPaths.back().readEnd = std::get<1>(pair.second[i]);
				readPaths.back().path.emplace_back(std::get<2>(pair.second[i]), true);
				pathLength.back() = morphConsensuses[std::get<2>(pair.second[i])].sequence.size();
			}
			else
			{
				readPaths.back().readEnd = std::get<1>(pair.second[i]);
				readPaths.back().path.emplace_back(std::get<2>(pair.second[i]), true);
				pathLength.back() += morphConsensuses[std::get<2>(pair.second[i])].sequence.size();
			}
		}
	}
	for (auto& pair : readBwMatches)
	{
		assert(pair.second.size() >= 1);
		std::sort(pair.second.begin(), pair.second.end(), [](auto left, auto right){ return std::get<0>(left) < std::get<0>(right); });
		for (size_t i = 0; i < pair.second.size(); i++)
		{
			if (i == 0 || std::get<1>(pair.second[i]) + 1000 < std::get<0>(pair.second[i-1]) || std::get<1>(pair.second[i]) > std::get<0>(pair.second[i-1]) + 1000)
			{
				readPaths.emplace_back();
				pathLength.emplace_back(0);
				readPaths.back().readLength = originalReadLength.at(pair.first);
				readPaths.back().readName = pair.first;
				readPaths.back().pathStartClip = 0;
				readPaths.back().pathEndClip = 0;
				readPaths.back().reverse = false;
				readPaths.back().readStart = std::get<1>(pair.second[i]);
				readPaths.back().readEnd = std::get<0>(pair.second[i]);
				readPaths.back().path.emplace_back(std::get<2>(pair.second[i]), false);
				pathLength.back() = morphConsensuses[std::get<2>(pair.second[i])].sequence.size();
			}
			else
			{
				readPaths.back().readEnd = std::get<0>(pair.second[i]);
				readPaths.back().path.emplace_back(std::get<2>(pair.second[i]), false);
				pathLength.back() += morphConsensuses[std::get<2>(pair.second[i])].sequence.size();
			}
		}
	}
	std::map<std::pair<Node, Node>, size_t> edgeCoverage;
	for (const auto& path : readPaths)
	{
		for (size_t i = 1; i < path.path.size(); i++)
		{
			auto key = canon(path.path[i-1], path.path[i]);
			edgeCoverage[key] += 1;
		}
	}
	{
		std::ofstream file { graphFile };
		for (size_t i = 0; i < morphConsensuses.size(); i++)
		{
			file << "S\t" << morphConsensuses[i].name << "\t*\tLN:i:" << morphConsensuses[i].sequence.size() << "\tll:f:" << morphConsensuses[i].coverage << "\tFC:i:" << (morphConsensuses[i].coverage * morphConsensuses[i].sequence.size()) << std::endl;
		}
		for (auto pair : edgeCoverage)
		{
			file << "L\t" << morphConsensuses[pair.first.first.id()].name << "\t" << (pair.first.first.forward() ? "+" : "-") << "\t" << morphConsensuses[pair.first.second.id()].name << "\t" << (pair.first.second.forward() ? "+" : "-") << "\t0M\tec:i:" << pair.second << std::endl;
		}
	}
	{
		std::ofstream file { pathsFile };
		for (size_t i = 0; i < readPaths.size(); i++)
		{
			file << readPaths[i].readName << "\t" << readPaths[i].readLength << " \t" << readPaths[i].readStart << "\t" << readPaths[i].readEnd << "\t+\t";
			for (auto node : readPaths[i].path)
			{
				file << (node.forward() ? ">" : "<") << morphConsensuses[node.id()].name;
			}
			file << "\t" << pathLength[i] << "\t" << 0 << "\t" << pathLength[i] << "\t" << pathLength[i] << "\t" << pathLength[i] << "\t" << 60 << std::endl;
		}
	}
}

template <typename T>
void filterVector(std::vector<T>& vec, const std::vector<bool>& kept)
{
	std::vector<size_t> newIndex;
	size_t totalKept = 0;
	for (size_t i = 0; i < kept.size(); i++)
	{
		newIndex.emplace_back(totalKept);
		if (kept[i]) totalKept += 1;
	}
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (!kept[i]) continue;
		vec[newIndex[i]] = vec[i];
	}
	vec.resize(totalKept);
}

size_t hammingDistance(const std::vector<uint8_t>& left, const std::vector<uint8_t>& right)
{
	assert(left.size() == right.size());
	size_t result = 0;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i] != right[i] || left[i] == 'N')
		{
			result += 1;
		}
	}
	return result;
}

bool sitesArePhaseInformative(const std::vector<std::vector<uint8_t>>& readSNPMSA, const size_t left, const size_t right, const size_t third, const size_t solidThreshold)
{
	phmap::flat_hash_map<std::tuple<uint8_t, uint8_t, uint8_t>, size_t> pairCounts;
	for (size_t i = 0; i < readSNPMSA.size(); i++)
	{
		pairCounts[std::make_tuple(readSNPMSA[i][left], readSNPMSA[i][right], readSNPMSA[i][third])] += 1;
	}
	phmap::flat_hash_set<std::tuple<uint8_t, uint8_t, uint8_t>> coveredPairs;
	for (auto pair : pairCounts)
	{
		if (pair.second < solidThreshold) return false;
		coveredPairs.emplace(pair.first);
	}
	if (coveredPairs.size() < 2) return false;
	for (auto pair : coveredPairs)
	{
		for (auto pair2 : coveredPairs)
		{
			if (pair == pair2) continue;
			size_t hammingDistance = 0;
			if (std::get<0>(pair) != std::get<0>(pair2)) hammingDistance += 1;
			if (std::get<1>(pair) != std::get<1>(pair2)) hammingDistance += 1;
			if (std::get<2>(pair) != std::get<2>(pair2)) hammingDistance += 1;
			assert(hammingDistance >= 1);
			if (hammingDistance == 1) return false;
		}
	}
	return true;
}

bool sitesArePhaseInformative(const std::vector<std::vector<uint8_t>>& readSNPMSA, const size_t left, const size_t right, const size_t solidThreshold)
{
	phmap::flat_hash_map<std::pair<uint8_t, uint8_t>, size_t> pairCounts;
	for (size_t i = 0; i < readSNPMSA.size(); i++)
	{
		pairCounts[std::make_pair(readSNPMSA[i][left], readSNPMSA[i][right])] += 1;
	}
	phmap::flat_hash_set<std::pair<uint8_t, uint8_t>> coveredPairs;
	for (auto pair : pairCounts)
	{
		if (pair.second < solidThreshold) return false;
		coveredPairs.emplace(pair.first);
	}
	if (coveredPairs.size() < 2) return false;
	for (auto pair : coveredPairs)
	{
		for (auto pair2 : coveredPairs)
		{
			if (pair == pair2) continue;
			if (pair.first == pair2.first) return false;
			if (pair.second == pair2.second) return false;
		}
	}
	return true;
}

std::pair<uint8_t, uint8_t> getBiallelicMatchAlleles(const std::vector<std::vector<uint8_t>>& MSA, const size_t left, const size_t right)
{
	phmap::flat_hash_map<std::pair<uint8_t, uint8_t>, size_t> pairCounts;
	for (size_t i = 0; i < MSA.size(); i++)
	{
		pairCounts[std::make_pair(MSA[i][left], MSA[i][right])] += 1;
	}
	std::vector<std::pair<std::pair<uint8_t, uint8_t>, size_t>> counts { pairCounts.begin(), pairCounts.end() };
	assert(counts.size() >= 2);
	std::sort(counts.begin(), counts.end(), [](const auto& left, const auto& right) {
		if (left.second > right.second) return true;
		if (left.second < right.second) return false;
		if (left.first.first < right.first.first) return true;
		return false;
	});
	assert(counts[0].second >= counts[1].second);
	assert(counts[0].first.first != counts[1].first.first);
	assert(counts[0].first.second != counts[1].first.second);
	if (counts.size() >= 3)
	{
		assert(counts[1].second >= counts[2].second * 2);
		counts.erase(counts.begin()+2, counts.end());
	}
	assert(counts.size() == 2);
	std::sort(counts.begin(), counts.end());
	assert(counts[0].first.first < counts[1].first.first);
	return std::make_pair(counts[0].first.second, counts[1].first.second);
}

bool sitesApproximatelyBiallelicallyMatch(const std::vector<std::vector<uint8_t>>& MSA, const size_t left, const size_t right)
{
	phmap::flat_hash_map<std::pair<uint8_t, uint8_t>, size_t> pairCounts;
	for (size_t i = 0; i < MSA.size(); i++)
	{
		pairCounts[std::make_pair(MSA[i][left], MSA[i][right])] += 1;
	}
	std::vector<std::pair<std::pair<uint8_t, uint8_t>, size_t>> counts { pairCounts.begin(), pairCounts.end() };
	std::sort(counts.begin(), counts.end(), [](const auto& left, const auto& right) { return left.second > right.second; });
	assert(counts[0].second >= counts[1].second);
	if (counts.size() < 2) return false;
	if (counts[0].second + counts[1].second < MSA.size() * 0.9) return false;
	if (counts.size() >= 3 && counts[1].second < counts[2].second * 2) return false;
	if (counts[0].second < MSA.size() * 0.1) return false;
	if (counts[1].second < MSA.size() * 0.1) return false;
	if (counts[0].second < 5) return false;
	if (counts[1].second < 5) return false;
	if (counts[0].first.first == counts[1].first.first) return false;
	if (counts[0].first.second == counts[1].first.second) return false;
	return true;
}

std::vector<std::vector<size_t>> tryBiallelicPhasing(const std::vector<std::vector<uint8_t>>& readSNPMSA)
{
	std::cerr << "coverage " << readSNPMSA.size() << " try biallelic phasing" << std::endl;
	std::vector<size_t> coveredSiteIndices;
	for (size_t i = 0; i < readSNPMSA[0].size(); i++)
	{
		phmap::flat_hash_map<uint8_t, size_t> alleleCounts;
		for (size_t j = 0; j < readSNPMSA.size(); j++)
		{
			alleleCounts[readSNPMSA[j][i]] += 1;
		}
		std::vector<std::pair<uint8_t, size_t>> alleles { alleleCounts.begin(), alleleCounts.end() };
		std::sort(alleles.begin(), alleles.end(), [](const auto& left, const auto& right) { return left.second > right.second; });
		if (alleles.size() < 2) continue;
		if (alleles.size() >= 3 && alleles[1].second < alleles[2].second * 2) continue;
		if (alleles[0].second + alleles[1].second < readSNPMSA.size() * 0.9) continue;
		if (alleles[0].second < 5) continue;
		if (alleles[1].second < 5) continue;
		if (alleles[0].first == '-') continue;
		if (alleles[0].first == 'N') continue;
		if (alleles[1].first == '-') continue;
		if (alleles[1].first == 'N') continue;
		coveredSiteIndices.emplace_back(i);
	}
	std::vector<bool> siteChecked;
	std::vector<std::vector<size_t>> alleleIndices;
	alleleIndices.resize(readSNPMSA.size());
	siteChecked.resize(coveredSiteIndices.size(), false);
	for (size_t i = 0; i < coveredSiteIndices.size(); i++)
	{
		if (siteChecked[i]) continue;
		std::vector<size_t> matchers;
		for (size_t j = 0; j < coveredSiteIndices.size(); j++)
		{
			if (i == j) continue;
			if (coveredSiteIndices[i]+50 > coveredSiteIndices[j] && coveredSiteIndices[j]+50 > coveredSiteIndices[i]) continue;
			if (sitesApproximatelyBiallelicallyMatch(readSNPMSA, coveredSiteIndices[i], coveredSiteIndices[j]))
			{
				matchers.emplace_back(j);
			}
		}
		if (matchers.size() >= 1)
		{
			std::cerr << "coverage " << readSNPMSA.size() << " biallelic phasing site " << coveredSiteIndices[i] << " count " << matchers.size() << ":";
			for (auto j : matchers)
			{
				std::cerr << " " << coveredSiteIndices[j];
			}
			std::cerr << std::endl;
		}
		if (matchers.size() < 2) continue; // at least three sites used for phasing
		std::sort(matchers.begin(), matchers.end());
		std::vector<uint8_t> alleleOne;
		std::vector<uint8_t> alleleTwo;
		{
			phmap::flat_hash_map<uint8_t, size_t> alleleCounts;
			for (size_t j = 0; j < readSNPMSA.size(); j++)
			{
				alleleCounts[readSNPMSA[j][coveredSiteIndices[i]]] += 1;
			}
			std::vector<std::pair<uint8_t, size_t>> alleles { alleleCounts.begin(), alleleCounts.end() };
			std::sort(alleles.begin(), alleles.end(), [](const auto& left, const auto& right) { return left.second > right.second; });
			assert(alleles.size() >= 2);
			if (alleles.size() >= 3)
			{
				assert(alleles[1].second >= alleles[2].second*2);
				alleles.erase(alleles.begin()+2, alleles.end());
			}
			std::sort(alleles.begin(), alleles.end());
			alleleOne.emplace_back(alleles[0].first);
			alleleTwo.emplace_back(alleles[1].first);
		}
		for (size_t j : matchers)
		{
			uint8_t alleleA = 0;
			uint8_t alleleB = 0;
			std::tie(alleleA, alleleB) = getBiallelicMatchAlleles(readSNPMSA, coveredSiteIndices[i], coveredSiteIndices[j]);
			assert(alleleA != 0);
			assert(alleleB != 0);
			assert(alleleA != alleleB);
			alleleOne.emplace_back(alleleA);
			alleleTwo.emplace_back(alleleB);
		}
		std::vector<size_t> assignments;
		size_t countSkipped = 0;
		size_t countOne = 0;
		size_t countTwo = 0;
		for (size_t j = 0; j < readSNPMSA.size(); j++)
		{
			std::vector<uint8_t> alleleHere;
			alleleHere.emplace_back(readSNPMSA[j][coveredSiteIndices[i]]);
			for (size_t k : matchers)
			{
				alleleHere.emplace_back(readSNPMSA[j][coveredSiteIndices[k]]);
			}
			size_t distanceOne = hammingDistance(alleleOne, alleleHere);
			size_t distanceTwo = hammingDistance(alleleTwo, alleleHere);
			if (distanceOne < distanceTwo)
			{
				countOne += 1;
				assignments.emplace_back(0);
			}
			else if (distanceTwo < distanceOne)
			{
				countTwo += 1;
				assignments.emplace_back(1);
			}
			else
			{
				assert(distanceOne == distanceTwo);
				assignments.emplace_back(2);
				countSkipped += 1;
			}
		}
		if (matchers.size() < 4 && countSkipped > 0) continue; // five matchers keeps everything, three/four only keeps very consistent phasings
		if (countSkipped > readSNPMSA.size() * 0.05) continue;
		if (countSkipped > 5) continue;
		if (countOne < 5) continue;
		if (countTwo < 5) continue;
		for (auto j : matchers)
		{
			siteChecked[j] = true;
		}
		std::cerr << "coverage " << readSNPMSA.size() << " biallelic phasing site " << coveredSiteIndices[i] << " does phase, counts " << countOne << " " << countTwo << ", skipped " << countSkipped << std::endl;
		for (size_t j = 0; j < assignments.size(); j++)
		{
			alleleIndices[j].emplace_back(assignments[j]);
		}
	}
	std::vector<size_t> parent;
	for (size_t j = 0; j < readSNPMSA.size(); j++)
	{
		parent.emplace_back(j);
	}
	phmap::flat_hash_map<size_t, size_t> clusterToIndex;
	size_t nextNum = 0;
	for (size_t j = 1; j < readSNPMSA.size(); j++)
	{
		for (size_t k = 0; k < j; k++)
		{
			if (alleleIndices[j] != alleleIndices[k]) continue;
			merge(parent, j, k);
		}
	}
	for (size_t i = 0; i < readSNPMSA.size(); i++)
	{
		if (clusterToIndex.count(find(parent, i)) == 1) continue;
		clusterToIndex[find(parent, i)] = nextNum;
		nextNum += 1;
	}
	assert(nextNum > 0);
	std::vector<std::vector<size_t>> result;
	result.resize(nextNum);
	for (size_t i = 0; i < readSNPMSA.size(); i++)
	{
		result[clusterToIndex.at(find(parent, i))].emplace_back(i);
	}
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) { return left.size() < right.size(); });
	std::cerr << "coverage " << readSNPMSA.size() << " biallelic phasing split to " << nextNum << " clusters" << std::endl;
	std::cerr << "sizes:";
	for (size_t i = 0; i < nextNum; i++)
	{
		std::cerr << " " << result[i].size();
	}
	std::cerr << std::endl;
	return result;
}

std::vector<std::vector<size_t>> SNPSplitAndAppend(const std::vector<OntLoop>& paths, const std::vector<std::string>& rawReads, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
{
	if (paths.size() < 10)
	{
//		std::cerr << "cluster with coverage " << rawReads.size() << " split to 1 cluster" << std::endl;
		std::vector<std::vector<size_t>> result;
		result.emplace_back();
		for (size_t i = 0; i < paths.size(); i++)
		{
			result.back().emplace_back(i);
		}
		return result;
	}
	assert(rawReads.size() == paths.size());
	assert(paths.size() >= 1);
	std::vector<std::vector<uint8_t>> readSNPMSA;
	readSNPMSA.resize(rawReads.size());
	auto path = getConsensusPath(paths, graph);
	std::string consensusSeq = getConsensusSequence(path, graph, pathStartClip, pathEndClip);
/*	{
		std::ofstream file { "ref_coverage_" + std::to_string(paths.size()) + ".fa" };
		file << ">ref_coverage_" << paths.size() << std::endl;
		file << consensusSeq << std::endl;
	}
	{
		std::ofstream file { "reads_coverage_" + std::to_string(paths.size()) + ".fa" };
		for (size_t i = 0; i < rawReads.size(); i++)
		{
			file << ">read_coverage_" << paths.size() << "_index_" << i << std::endl;
			file << rawReads[i] << std::endl;
		}
	}*/
	bool anyAlignmentError = false;
	size_t sequenceLengthSum = 0;
	size_t totalMismatches = 0;
	for (size_t i = 0; i < rawReads.size(); i++)
	{
		EdlibAlignResult result = edlibAlign(consensusSeq.data(), consensusSeq.size(), rawReads[i].data(), rawReads[i].size(), edlibNewAlignConfig(std::max(consensusSeq.size(), rawReads[i].size()), EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		if (result.status != EDLIB_STATUS_OK)
		{
			anyAlignmentError = true;
			edlibFreeAlignResult(result);
			break;
		}
		assert(result.numLocations >= 1);
		assert(result.endLocations[0] <= rawReads[i].size());
		assert(result.endLocations[0] > result.startLocations[0]);
		sequenceLengthSum += result.endLocations[0] - result.startLocations[0] + 1;
		totalMismatches += result.editDistance;
		size_t consensusPos = 0;
		size_t seqPos = result.startLocations[0];
		readSNPMSA[i].resize(consensusSeq.size(), 'N');
		assert(result.alignmentLength >= std::max((size_t)(result.endLocations[0]-result.startLocations[0]), consensusSeq.size()));
		for (size_t k = 0; k < result.alignmentLength; k++)
		{
			switch(result.alignment[k])
			{
			case 0:
				assert(seqPos < rawReads[i].size());
				assert(consensusPos < consensusSeq.size());
				if (!(rawReads[i][seqPos] == consensusSeq[consensusPos]))
				{
					std::cerr << rawReads[i][seqPos] << " " << consensusSeq[consensusPos] << " " << seqPos << " " << consensusPos << " " << rawReads.size() << " " << i << std::endl;
				}
				assert(rawReads[i][seqPos] == consensusSeq[consensusPos]);
				readSNPMSA[i][consensusPos] = rawReads[i][seqPos];
				seqPos += 1;
				consensusPos += 1;
				continue;
			case 1:
				assert(consensusPos < consensusSeq.size());
				readSNPMSA[i][consensusPos] = '-';
				consensusPos += 1;
				continue;
			case 2:
				assert(seqPos < rawReads[i].size());
				seqPos += 1;
				continue;
			case 3:
				assert(seqPos < rawReads[i].size());
				assert(consensusPos < consensusSeq.size());
				assert(rawReads[i][seqPos] != consensusSeq[consensusPos]);
				readSNPMSA[i][consensusPos] = rawReads[i][seqPos];
				seqPos += 1;
				consensusPos += 1;
				continue;
			default:
				assert(false);
			}
		}
		assert(seqPos == result.endLocations[0]+1);
		assert(consensusPos == consensusSeq.size());
		edlibFreeAlignResult(result);
	}
	assert(!anyAlignmentError);
	double approxErrorRate = (double)totalMismatches / (double)sequenceLengthSum;
//	std::cerr << "cluster with coverage " << rawReads.size() << " approx error rate " << approxErrorRate << std::endl;
	std::vector<bool> coveredSite;
	coveredSite.resize(consensusSeq.size(), false);
	std::vector<size_t> coveredSiteIndices;
	for (size_t i = 0; i < consensusSeq.size(); i++)
	{
		phmap::flat_hash_map<size_t, size_t> alleleCounts;
		for (size_t j = 0; j < readSNPMSA.size(); j++)
		{
			alleleCounts[readSNPMSA[j][i]] += 1;
		}
		size_t countCovered = 0;
		size_t countN = 0;
		for (auto pair : alleleCounts)
		{
//			std::cerr << "cluster with coverage " << rawReads.size() << " site " << i << " allele " << (char)pair.first << " count " << pair.second << std::endl;
			if (pair.first == 'N')
			{
				countN = pair.second;
				continue;
			}
			if (pair.first == '-') continue;
			if (pair.second < 5) continue;
			if (pair.second < rawReads.size() * approxErrorRate + 1) continue;
			if (pair.first != consensusSeq[i])
			{
				if ((i > 0 && consensusSeq[i-1] == consensusSeq[i]) || (i+1 < consensusSeq.size() && consensusSeq[i+1] == consensusSeq[i]))
				{
					if ((i > 0 && consensusSeq[i-1] == pair.first) || (i+1 < consensusSeq.size() && consensusSeq[i+1] == pair.first))
					{
//						std::cerr << "hpc" << std::endl;
						continue;
					}
				}
			}
			countCovered += 1;
		}
		if (countCovered < 2) continue;
		if (countN > 5) continue;
		coveredSite[i] = true;
		coveredSiteIndices.emplace_back(i);
//		std::cerr << "covered " << i << std::endl;
	}
//	std::cerr << "coverage " << rawReads.size() << " covered sites before phasing filtering " << coveredSiteIndices.size() << std::endl;
	std::vector<bool> phaseInformativeSite;
	phaseInformativeSite.resize(coveredSite.size(), false);
	for (size_t i = 0; i < coveredSiteIndices.size(); i++)
	{
		assert(coveredSite[coveredSiteIndices[i]]);
		for (size_t j = 0; j < i; j++)
		{
			assert(coveredSite[coveredSiteIndices[j]]);
			if (phaseInformativeSite[coveredSiteIndices[i]] && phaseInformativeSite[coveredSiteIndices[j]]) continue;
			if (sitesArePhaseInformative(readSNPMSA, coveredSiteIndices[i], coveredSiteIndices[j], std::max<size_t>(5, rawReads.size() * approxErrorRate + 1)))
			{
//				std::cerr << "phase informative two " << coveredSiteIndices[i] << " " << coveredSiteIndices[j] << std::endl;
				phaseInformativeSite[coveredSiteIndices[i]] = true;
				phaseInformativeSite[coveredSiteIndices[j]] = true;
			}
		}
	}
	for (size_t i = 2; i < coveredSiteIndices.size(); i++)
	{
		for (size_t j = 1; j < i ; j++)
		{
			for (size_t k = 0; k < j; k++)
			{
				if (phaseInformativeSite[coveredSiteIndices[i]] && phaseInformativeSite[coveredSiteIndices[j]] && phaseInformativeSite[coveredSiteIndices[k]]) continue;
				if (sitesArePhaseInformative(readSNPMSA, coveredSiteIndices[i], coveredSiteIndices[j], coveredSiteIndices[k], std::max<size_t>(5, rawReads.size() * approxErrorRate + 1)))
				{
//					std::cerr << "phase informative three " << coveredSiteIndices[i] << " " << coveredSiteIndices[j] << " " << coveredSiteIndices[k] << std::endl;
					phaseInformativeSite[coveredSiteIndices[i]] = true;
					phaseInformativeSite[coveredSiteIndices[j]] = true;
					phaseInformativeSite[coveredSiteIndices[k]] = true;
				}
			}
		}
	}
/*	std::cerr << "cluster with coverage " << rawReads.size() << " covered sites:";
	for (size_t i = 0; i < coveredSite.size(); i++)
	{
		if (!coveredSite[i]) continue;
		std::cerr << " " << i;
	}
	std::cerr << std::endl;
	std::cerr << "cluster with coverage " << rawReads.size() << " total sites " << readSNPMSA[0].size() << std::endl;*/
	auto filteredMSA = readSNPMSA;
	for (size_t i = 0; i < readSNPMSA.size(); i++)
	{
		filterVector(filteredMSA[i], phaseInformativeSite);
	}
//	std::cerr << "cluster with coverage " << rawReads.size() << " count covered sites " << filteredMSA[0].size() << std::endl;
	if (filteredMSA[0].size() < 2)
	{
//		std::cerr << "cluster with coverage " << rawReads.size() << " split to 1 cluster" << std::endl;
		return tryBiallelicPhasing(readSNPMSA);
	}
	std::vector<size_t> parent;
	for (size_t i = 0; i < rawReads.size(); i++)
	{
		parent.emplace_back(i);
	}
/*	{
		std::ofstream file { "MSA_coverage" + std::to_string(rawReads.size()) + ".txt" };
		for (size_t i = 0; i < rawReads.size(); i++)
		{
			for (size_t j = 0; j < filteredMSA[i].size(); j++)
			{
				file << (char)filteredMSA[i][j];
			}
			file << std::endl;
		}
	}*/
	assert(filteredMSA[0].size() >= 2);
	size_t maxErrorsToMerge = filteredMSA[0].size() * 0.5;
	assert(maxErrorsToMerge >= 1);
//	std::cerr << "cluster with coverage " << rawReads.size() << " max errors " << maxErrorsToMerge << std::endl;
	for (size_t i = 1; i < rawReads.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (find(parent, i) == find(parent, j)) continue;
//			if (i == 1 && j == 0) std::cerr << "hamming distance " << hammingDistance(filteredMSA[i], filteredMSA[j]) << " vs " << maxErrorsToMerge << std::endl;
			if (hammingDistance(filteredMSA[i], filteredMSA[j]) > maxErrorsToMerge) continue;
			merge(parent, i, j);
		}
	}
	phmap::flat_hash_map<size_t, size_t> clusterToIndex;
	size_t nextNum = 0;
	for (size_t i = 0; i < rawReads.size(); i++)
	{
		if (clusterToIndex.count(find(parent, i)) == 1) continue;
		clusterToIndex[find(parent, i)] = nextNum;
		nextNum += 1;
	}
	assert(nextNum > 0);
	if (nextNum == 1)
	{
		return tryBiallelicPhasing(readSNPMSA);
	}
	if (nextNum == 1)
	{
//		std::cerr << "cluster with coverage " << rawReads.size() << " split to 1 cluster" << std::endl;
		std::vector<std::vector<size_t>> result;
		result.emplace_back();
		for (size_t i = 0; i < paths.size(); i++)
		{
			result.back().emplace_back(i);
		}
		return result;
	}
	std::vector<std::vector<size_t>> result;
	result.resize(nextNum);
	for (size_t i = 0; i < paths.size(); i++)
	{
		result[clusterToIndex.at(find(parent, i))].emplace_back(i);
	}
/*	std::cerr << "cluster with coverage " << rawReads.size() << " split to " << nextNum << " clusters" << std::endl;
	std::cerr << "sizes:";
	for (size_t i = 0; i < nextNum; i++)
	{
		std::cerr << " " << result[i].size();
	}
	std::cerr << std::endl;*/
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) { return left.size() < right.size(); });
	return result;
}

std::vector<std::vector<std::string>> getRawLoopSequences(const std::vector<std::vector<OntLoop>>& clusters, const std::string& rawOntPath)
{
	phmap::flat_hash_map<std::string, std::vector<std::tuple<size_t, size_t, size_t, size_t>>> loopPositionsInReads;
	std::vector<std::vector<std::string>> result;
	result.resize(clusters.size());
	for (size_t i = 0; i < clusters.size(); i++)
	{
		result[i].resize(clusters[i].size(), "");
		for (size_t j = 0; j < clusters[i].size(); j++)
		{
			loopPositionsInReads[clusters[i][j].readName].emplace_back(i, j, clusters[i][j].approxStart, clusters[i][j].approxEnd);
		}
	}
	FastQ::streamFastqFromFile(rawOntPath, false, [&result, &loopPositionsInReads](FastQ& fastq)
	{
		std::string readname = nameWithoutTags(fastq.seq_id);
		if (loopPositionsInReads.count(readname) == 0) return;
		for (auto t : loopPositionsInReads.at(readname))
		{
			size_t startPos = std::get<2>(t);
			size_t endPos = std::get<3>(t);
			bool reverse = false;
			if (startPos > endPos)
			{
				std::swap(startPos, endPos);
				reverse = true;
			}
			if (startPos > 500)
			{
				startPos -= 500;
			}
			else
			{
				startPos = 0;
			}
			if (endPos+500 < fastq.sequence.size())
			{
				endPos += 500;
			}
			else
			{
				endPos = fastq.sequence.size();
			}
			std::string str = fastq.sequence.substr(startPos, endPos-startPos);
			if (reverse)
			{
				str = revcomp(str);
			}
			assert(result[std::get<0>(t)][std::get<1>(t)] == "");
			assert(str.size() >= 1);
			result[std::get<0>(t)][std::get<1>(t)] = str;
		}
	});
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].size(); j++)
		{
			assert(result[i][j] != "");
		}
	}
	return result;
}

std::vector<std::vector<OntLoop>> SNPsplitClusters(const std::vector<std::vector<OntLoop>>& unphasedClusters, const std::string rawOntPath, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t numThreads)
{
	std::vector<std::pair<std::vector<OntLoop>, std::vector<std::string>>> checkStack;
	checkStack.resize(unphasedClusters.size());
	for (size_t i = 0; i < unphasedClusters.size(); i++)
	{
		checkStack[i].first = unphasedClusters[i];
	}
	{
		auto rawSequences = getRawLoopSequences(unphasedClusters, rawOntPath);
		assert(rawSequences.size() == unphasedClusters.size());
		assert(rawSequences.size() == checkStack.size());
		for (size_t i = 0; i < unphasedClusters.size(); i++)
		{
			assert(checkStack[i].first.size() == rawSequences[i].size());
			std::swap(checkStack[i].second, rawSequences[i]);
		}
	}
	for (size_t i = 0; i < unphasedClusters.size(); i++)
	{
		assert(checkStack[i].first.size() == checkStack[i].second.size());
		for (size_t j = 0; j < unphasedClusters[i].size(); j++)
		{
			assert(checkStack[i].second[j].size() >= 1);
		}
	}
	std::vector<std::vector<OntLoop>> result;
	std::vector<std::vector<std::string>> rawSequencesPerResult;
	std::vector<std::thread> threads;
	std::mutex popMutex;
	std::atomic<size_t> threadsRunning;
	threadsRunning = 0;
	std::sort(checkStack.begin(), checkStack.end(), [](const auto& left, const auto& right) { return left.first.size() < right.first.size(); });
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&result, &checkStack, &popMutex, &threadsRunning, &graph, &pathStartClip, &pathEndClip, &rawSequencesPerResult](){
			while (true)
			{
				std::vector<OntLoop> loops;
				std::vector<std::string> rawSequences;
				{
					std::lock_guard<std::mutex> lock { popMutex };
					if (checkStack.size() == 0)
					{
						if (threadsRunning == 0)
						{
							return;
						}
					}
					else
					{
						std::swap(loops, checkStack.back().first);
						std::swap(rawSequences, checkStack.back().second);
						checkStack.pop_back();
						threadsRunning += 1;
					}
				}
				if (loops.size() == 0)
				{
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
					continue;
				}
				assert(loops.size() == rawSequences.size());
				auto splitted = SNPSplitAndAppend(loops, rawSequences, graph, pathStartClip, pathEndClip);
				assert(splitted.size() >= 1);
				if (splitted.size() == 1)
				{
					std::lock_guard<std::mutex> lock { popMutex };
					result.emplace_back();
					std::swap(result.back(), loops);
					rawSequencesPerResult.emplace_back();
					std::swap(rawSequencesPerResult.back(), rawSequences);
					threadsRunning -= 1;
					continue;
				}
				{
					std::lock_guard<std::mutex> lock { popMutex };
					for (size_t i = 0; i < splitted.size(); i++)
					{
						assert(splitted[i].size() >= 1);
						checkStack.emplace_back();
						for (size_t n : splitted[i])
						{
							checkStack.back().first.emplace_back();
							checkStack.back().second.emplace_back();
							std::swap(checkStack.back().first.back(), loops[n]);
							std::swap(checkStack.back().second.back(), rawSequences[n]);
						}
						for (size_t j = 0; j < checkStack.back().first.size(); j++)
						{
							assert(checkStack.back().second[j].size() != 0);
						}
					}
					threadsRunning -= 1;
					continue;
				}
			}
		});
	}
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		auto path = getConsensusPath(result[i], graph);
		assert(rawSequencesPerResult[i].size() == result[i].size());
		std::string consensusSeq = getConsensusSequence(path, graph, pathStartClip, pathEndClip);
		std::ofstream file { "ref_" + std::to_string(i) + "_coverage_" + std::to_string(result[i].size()) + ".fa" };
		file << ">ref_" << i << "_coverage_" << result[i].size() << std::endl;
		file << consensusSeq << std::endl;
		std::ofstream file2 { "reads_" + std::to_string(i) + "_coverage_" + std::to_string(result[i].size()) + ".fa" };
		for (size_t j = 0; j < rawSequencesPerResult[i].size(); j++)
		{
			file2 << ">read_" << i << "_" << j << std::endl;
			file2 << rawSequencesPerResult[i][j] << std::endl;
		}
	}
	return result;
}

std::vector<std::vector<std::pair<std::string, std::string>>> getRawOntLoops(const std::vector<MorphConsensus>& morphConsensuses, const std::string rawOntPath)
{
	std::vector<std::vector<OntLoop>> loops;
	loops.resize(morphConsensuses.size());
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		loops[i] = morphConsensuses[i].ontLoops;
	}
	std::vector<std::vector<std::string>> sequences = getRawLoopSequences(loops, rawOntPath);
	std::vector<std::vector<std::pair<std::string, std::string>>> result;
	result.resize(sequences.size());
	for (size_t i = 0; i < sequences.size(); i++)
	{
		result[i].resize(sequences[i].size());
		for (size_t j = 0; j < sequences[i].size(); j++)
		{
			std::swap(result[i][j].first, sequences[i][j]);
			result[i][j].second = "morph_" + std::to_string(i) + "_rawloop_" + std::to_string(j) + "_" + loops[i][j].readName + "_" + std::to_string(loops[i][j].approxStart) + "_" + std::to_string(loops[i][j].approxEnd);
		}
	}
	return result;
}

void writeRawOntLoopSequences(const std::string outputFile, const std::vector<std::vector<std::pair<std::string, std::string>>>& loopSequences)
{
	std::ofstream file { outputFile };
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		for (size_t j = 0; j < loopSequences[i].size(); j++)
		{
			file << ">" << loopSequences[i][j].second << std::endl;
			file << loopSequences[i][j].first << std::endl;
		}
	}
}

void alignRawOntLoopsToMorphConsensuses(const std::vector<std::vector<std::pair<std::string, std::string>>>& loopSequences, std::string outputFile, std::string tmppath, const size_t numThreads, const std::vector<MorphConsensus>& morphConsensuses, const std::string& winnowmapPath, const std::string& samtoolsPath)
{
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		{
			std::ofstream tmpReadsFile { tmppath + "/tmpreads.fa" };
			std::ofstream tmpConsensusFile { tmppath + "/tmpconsensus.fa" };
			for (size_t j = 0; j < loopSequences[i].size(); j++)
			{
				tmpReadsFile << ">" << loopSequences[i][j].second << std::endl;
				tmpReadsFile << loopSequences[i][j].first << std::endl;
			}
			tmpConsensusFile << ">" << morphConsensuses[i].name << std::endl;
			tmpConsensusFile << morphConsensuses[i].sequence << std::endl;
		}
		std::string winnowmapCommand = winnowmapPath + " -x map-ont -a -t " + std::to_string(numThreads) + " " + tmppath + "/tmpconsensus.fa " + tmppath + "/tmpreads.fa | " + samtoolsPath + " view -b | " + samtoolsPath + " sort > " + tmppath + "/tmpalns_" + std::to_string(i) + ".bam";
		std::cerr << "winnowmap command:" << std::endl;
		std::cerr << winnowmapCommand << std::endl;
		int result = system(winnowmapCommand.c_str());
		if (result != 0)
		{
			std::cerr << "alignment did not run successfully" << std::endl;
			std::abort();
		}
	}
	std::string mergeCommand = samtoolsPath + " merge -o " + outputFile;
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		mergeCommand += " " + tmppath + "/tmpalns_" + std::to_string(i) + ".bam";
	}
	std::cerr << "samtools command:" << std::endl;
	std::cerr << mergeCommand << std::endl;
	int result = system(mergeCommand.c_str());
	if (result != 0)
	{
		std::cerr << "samtools merge did not run successfully" << std::endl;
		std::abort();
	}
	std::string indexCommand = samtoolsPath + " index -b " + outputFile;
	std::cerr << "samtools command:" << std::endl;
	std::cerr << indexCommand << std::endl;
	result = system(indexCommand.c_str());
	if (result != 0)
	{
		std::cerr << "samtools index did not run successfully" << std::endl;
		std::abort();
	}
}

void DoClusterONTAnalysis(const ClusterParams& params)
{
	std::filesystem::create_directories(params.basePath + "/tmp");
	std::cerr << "reading allele graph" << std::endl;
	GfaGraph graph;
	graph.loadFromFile(params.basePath + "/allele-graph.gfa");
	std::cerr << "reading consensus" << std::endl;
	Path heavyPath = readHeavyPath(graph, params.basePath + "/consensus_path.gaf");
	std::cerr << "extract corrected ultralong paths" << std::endl;
	size_t heavyPathLength = heavyPath.getSequence(graph.nodeSeqs).size();
	size_t minLength = heavyPathLength * 0.5;
	std::cerr << "consensus path length " << heavyPathLength << ", using " << minLength << " as minimum morph length" << std::endl;
	auto ontPaths = extractCorrectedONTPaths(params.basePath + "/ont-alns.gaf", heavyPath, minLength, graph);
	std::cerr << ontPaths.size() << " corrected paths" << std::endl;
	std::cerr << "extract loops from ONTs" << std::endl;
	std::unordered_map<Node, size_t> pathStartClip;
	std::unordered_map<Node, size_t> pathEndClip;
	std::unordered_set<size_t> borderNodes;
	std::tie(borderNodes, pathStartClip, pathEndClip) = getBorderNodes(heavyPath, graph);
	assert(borderNodes.size() > 0);
	auto loopSequences = extractLoopSequences(ontPaths, heavyPath, minLength, graph, borderNodes, pathStartClip, pathEndClip);
	auto coreNodes = getCoreNodes(loopSequences);
	std::cerr << loopSequences.size() << " loops in ONTs" << std::endl;
	{
		std::ofstream file { params.basePath + "/loops.fa" };
		for (size_t i = 0; i < loopSequences.size(); i++)
		{
			file << ">loop_" << i << "_read_" << loopSequences[i].readName << "_start_" << loopSequences[i].approxStart << "_end_" << loopSequences[i].approxEnd << std::endl;
			std::string loopSeq = getSequence(loopSequences[i].path, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
			size_t leftClipBp = pathStartClip.at(loopSequences[i].path[0]);
			size_t rightClipBp = pathEndClip.at(loopSequences[i].path.back());
			loopSeq = loopSeq.substr(leftClipBp, loopSeq.size() - leftClipBp - rightClipBp);
			file << loopSeq << std::endl;
		}
	}
	std::cerr << "cluster loops roughly" << std::endl;
	orderLoopsByLength(loopSequences, graph, pathStartClip, pathEndClip);
	std::cerr << "max clustering edit distance " << params.maxClusterDifference << std::endl;
	auto unphasedClusters = roughClusterLoopSequences(loopSequences, graph, pathStartClip, pathEndClip, coreNodes, params.maxClusterDifference);
	std::cerr << unphasedClusters.size() << " rough clusters" << std::endl;
	std::cerr << "cluster loops by density" << std::endl;
	auto clusters = densityClusterLoops(unphasedClusters, graph, pathStartClip, pathEndClip, coreNodes, params.maxClusterDifference, 5, params.minReclusterDistance);
	std::cerr << clusters.size() << " density clusters" << std::endl;
	std::cerr << "phase clusters" << std::endl;
	clusters = phaseClusters(clusters);
	std::cerr << clusters.size() << " phased clusters" << std::endl;
	std::cerr << "SNP-split clusters" << std::endl;
	clusters = SNPsplitClusters(clusters, params.ontReadPath, graph, pathStartClip, pathEndClip, params.numThreads);
	std::cerr << clusters.size() << " phased clusters" << std::endl;
	std::sort(clusters.begin(), clusters.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
	std::cerr << "getting morph consensuses" << std::endl;
	auto morphConsensuses = getMorphConsensuses(clusters, graph, pathStartClip, pathEndClip, params.namePrefix);
	std::cerr << "write morph consensuses" << std::endl;
	writeMorphConsensuses(params.basePath + "/morphs.fa", morphConsensuses);
	std::cerr << "write morph paths" << std::endl;
	writeMorphPaths(params.basePath + "/morphs.gaf", morphConsensuses, graph, pathStartClip, pathEndClip);
	{
		auto ontLoopSequences = getRawOntLoops(morphConsensuses, params.ontReadPath);
		std::cerr << "write raw ONT loop sequences" << std::endl;
		writeRawOntLoopSequences(params.basePath + "/raw_loops.fa", ontLoopSequences);
		if (params.winnowmapPath != "" && params.samtoolsPath != "")
		{
			std::cerr << "realign raw ONT loop sequences to morph consensuses" << std::endl;
			alignRawOntLoopsToMorphConsensuses(ontLoopSequences, params.basePath + "/raw_loop_to_morphs_alignments.bam", params.basePath + "/tmp", params.numThreads, morphConsensuses, params.winnowmapPath, params.samtoolsPath);
		}
		else
		{
			if (params.winnowmapPath == "")
			{
				std::cerr << "winnowmap not found, ";
			}
			if (params.samtoolsPath != "")
			{
				std::cerr << "bamtools not found, ";
			}
			std::cerr << "skipping alignment of raw ONT loop sequences to morph consensuses" << std::endl;
		}
	}
	std::cerr << "write morph graph and read paths" << std::endl;
	writeMorphGraphAndReadPaths(params.basePath + "/morphgraph.gfa", params.basePath + "/readpaths-morphgraph.gaf", morphConsensuses);
	if (params.annotationFasta.size() > 0)
	{
		std::cerr << "lifting over annotations to morphs" << std::endl;
		liftoverAnnotationsToMorphs(params.basePath, morphConsensuses, params.annotationFasta, params.annotationGff3);
	}
}
