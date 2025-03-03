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

std::chrono::time_point<std::chrono::steady_clock> getTime()
{
	return std::chrono::steady_clock::now();
}

std::string formatTime(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end)
{
	size_t milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
	return std::to_string(milliseconds / 1000) + "," + std::to_string(milliseconds % 1000) + " s";
}

class CyclicGraphException : public std::exception
{
};

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
	template <> struct hash<__uint128_t>
	{
		size_t operator()(const __uint128_t& x) const
		{
			return std::hash<size_t>{}(x) ^ std::hash<size_t>{}(x >> 64);
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
	std::string rawSequence;
	std::string rawLoopName;
	std::string selfCorrectedSequence;
};

class MorphConsensus
{
public:
	std::vector<Node> path;
	std::string sequence;
	std::string preCorrectionSequence;
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
	bool nonrDNAAnchorVariant;
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

__attribute__((always_inline)) inline size_t getBitvectorMatchLength(uint64_t leftFirst, uint64_t leftSecond, uint64_t rightFirst, uint64_t rightSecond)
{
	uint64_t mismatch = (leftFirst ^ rightFirst) | (leftSecond ^ rightSecond);
	uint64_t prefixZerosPlusOneOne = (mismatch-1) ^ mismatch;
	uint64_t prefixMatchPlusOne = popcount(prefixZerosPlusOneOne);
	assert(prefixMatchPlusOne >= 1);
	return prefixMatchPlusOne-1;
}

__attribute__((always_inline)) inline size_t getMatchLength(const std::string& left, const std::string& right, const size_t leftStart, const size_t rightStart)
{
	if (leftStart == left.size() || rightStart == right.size()) return 0;
	size_t result = 0;
	while (leftStart+result < left.size() && rightStart+result < right.size())
	{
		if (left[leftStart+result] != right[rightStart+result]) break;
		result += 1;
	}
	return result;
}

__attribute__((always_inline)) inline size_t getMatchLength(const TwobitString& left, const TwobitString& right, const size_t leftStart, const size_t rightStart)
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
	if (predecessor.count(start) == 0)
	{
		throw CyclicGraphException {};
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

std::vector<Node> getShortestPath(const phmap::flat_hash_map<Node, Node>& shortestPathPredecessor, const Node startNode, const size_t endNode)
{
	std::vector<Node> result;
	result.emplace_back(startNode);
	while (true)
	{
		assert(shortestPathPredecessor.count(result.back()) == 1);
		result.emplace_back(shortestPathPredecessor.at(result.back()));
		if (result.back().id() == endNode) return result;
		assert(result.size() < shortestPathPredecessor.size() + 5);
	}
	assert(false);
	return result;
}

void filterOutLocalCycles(std::unordered_set<size_t>& coveredNodes, std::unordered_map<Node, std::unordered_set<Node>>& coveredEdges, const size_t topCoverageNode, const GfaGraph& graph)
{
	phmap::flat_hash_map<Node, Node> shortestPathPredecessor;
	std::vector<std::tuple<size_t, Node, Node>> checkStack;
	shortestPathPredecessor[Node { topCoverageNode, true }] = Node { topCoverageNode, true };
	shortestPathPredecessor[Node { topCoverageNode, false }] = Node { topCoverageNode, false };
	for (auto node : coveredEdges.at(Node { topCoverageNode, true }))
	{
		assert(coveredNodes.count(node.id()) == 1);
		size_t overlap = getOverlap(Node { topCoverageNode, true }, node, graph.edges);
		checkStack.emplace_back(graph.nodeSeqs[node.id()].size() - overlap, Node { topCoverageNode, true }, node);
	}
	for (auto node : coveredEdges.at(Node { topCoverageNode, false }))
	{
		assert(coveredNodes.count(node.id()) == 1);
		size_t overlap = getOverlap(Node { topCoverageNode, false }, node, graph.edges);
		checkStack.emplace_back(graph.nodeSeqs[node.id()].size() - overlap, Node { topCoverageNode, false }, node);
	}
	std::sort(checkStack.begin(), checkStack.end(), [](auto left, auto right) { return std::get<0>(left) > std::get<0>(right); });
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (shortestPathPredecessor.count(std::get<2>(top)) == 1) continue;
		shortestPathPredecessor[std::get<2>(top)] = std::get<1>(top);
		assert(coveredNodes.count(std::get<2>(top).id()) == 1);
		assert(coveredEdges.count(std::get<2>(top)) == 1);
		for (auto edge : coveredEdges.at(std::get<2>(top)))
		{
			size_t overlap = getOverlap(std::get<2>(top), edge, graph.edges);
			checkStack.emplace_back(graph.nodeSeqs[edge.id()].size() - overlap, std::get<2>(top), edge);
		}
	}
	std::unordered_set<size_t> repeatNodes;
	for (size_t node : coveredNodes)
	{
		assert(shortestPathPredecessor.count(Node { node, true }) == 1);
		assert(shortestPathPredecessor.count(Node { node, false }) == 1);
		std::vector<Node> shortestPathToHere = getShortestPath(shortestPathPredecessor, Node { node, true }, topCoverageNode);
		std::vector<Node> shortestPathFromHere = getShortestPath(shortestPathPredecessor, Node { node, false }, topCoverageNode);
		phmap::flat_hash_set<size_t> nodesInPath;
		for (Node pathnode : shortestPathToHere)
		{
			if (pathnode.id() == node) continue;
			if (pathnode.id() == topCoverageNode) continue;
			assert(nodesInPath.count(pathnode.id()) == 0);
			nodesInPath.emplace(pathnode.id());
		}
		for (Node pathnode : shortestPathFromHere)
		{
			if (pathnode.id() == node) continue;
			if (pathnode.id() == topCoverageNode) continue;
			if (nodesInPath.count(pathnode.id()) == 1)
			{
				repeatNodes.emplace(node);
				break;
			}
		}
	}
	std::cerr << "filter out " << repeatNodes.size() << " local repeat nodes" << std::endl;
	filterOut(coveredNodes, repeatNodes);
	filterOut(coveredEdges, repeatNodes);
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
	Path result;
	try
	{
		result = getRecWidestPath(coveredNodes, coveredEdges, graph, Node { topCoverageNode, true });
	}
	catch (CyclicGraphException& e)
	{
		std::cerr << "cyclic heavy path, remove local cycles and retry" << std::endl;
		filterOutLocalCycles(coveredNodes, coveredEdges, topCoverageNode, graph);
		while (true)
		{
			size_t sizeBeforeFilter = coveredNodes.size();
			filterOutDefinitelyNonConsensusNodes(coveredNodes, coveredEdges, topCoverageNode, graph, localResolveLength);
			filterOutNonCircularParts(coveredNodes, coveredEdges, topCoverageNode);
			if (coveredNodes.size() == sizeBeforeFilter) break;
		}
		result = getRecWidestPath(coveredNodes, coveredEdges, graph, Node { topCoverageNode, true });
	}
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
		here.readName = readname;
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
		result.back().nonrDNAAnchorVariant = false;
	}
	return result;
}

void writeVariants(const Path& heavyPath, const GfaGraph& graph, const std::vector<Variant>& variants, const std::string& filename)
{
	std::ofstream file { filename };
	for (const auto& variant : variants)
	{
		if (variant.nonrDNAAnchorVariant) continue;
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
		for (size_t j = 1; j < variants[i].path.size() - (variants[i].nonrDNAAnchorVariant ? 0 : 1); j++)
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
			if (j != variants[i].path.size()-1 || variants[i].nonrDNAAnchorVariant)
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
		if (variants[i].path.size() == 2 && !variants[i].nonrDNAAnchorVariant)
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
		if (variants[i].path.size() == 2 && variants[i].nonrDNAAnchorVariant)
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
			assert(edgeCoverage.count(std::make_pair(fromname, toname)) == 1);
			out << "S\t" << toname << "\t" << fullGraph.nodeSeqs[variants[i].path[1].id()] << "\tll:f:" << nodeCoverage.at(key) << "\tFC:f:" << (nodeCoverage.at(key) * fullGraph.nodeSeqs[variants[i].path[1].id()].size()) << std::endl;
			out << "L\t" << fromname << "\t" << (prev.forward() ? "+" : "-") << "\t" << toname << "\t" << (curr.forward() ? "+" : "-") << "\t" << overlap << "M" << "\tec:i:" << edgeCoverage.at(std::make_pair(fromname, toname)) << std::endl;
			continue;
		}
		assert(referenceNode[variants[i].path[0].id()]);
		assert(referenceNode[variants[i].path.back().id()] || variants[i].nonrDNAAnchorVariant);
		for (size_t j = 1; j < variants[i].path.size() - (variants[i].nonrDNAAnchorVariant ? 0 : 1); j++)
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
			assert(edgeCoverage.count(std::make_pair(fromname, toname)) == 1);
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
			assert(edgeCoverage.count(std::make_pair(fromname, toname)) == 1);
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
			if (variants[i].nonrDNAAnchorVariant)
			{
				std::pair<size_t, size_t> tokey = parent[i][variants[i].path.size()-1];
				while (parent[tokey.first][tokey.second] != tokey)
				{
					tokey = parent[tokey.first][tokey.second];
				}
				toname = fullGraph.nodeNames[curr.id()] + "_" + std::to_string(tokey.first) + "_" + std::to_string(tokey.second);
			}
			assert(edgeCoverage.count(std::make_pair(fromname, toname)) == 1);
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
		if (variants[variant].nonrDNAAnchorVariant) continue;
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
		if (variants[i].nonrDNAAnchorVariant) continue;
		file << "heavy_path\t" << variants[i].referenceStartPos+1 << "\t" << variants[i].name << "\t" << variants[i].referenceSeq << "\t" << variants[i].variantSeq << "\t" << "." << "\tPASS\t" << "AD=" << variants[i].referenceCoverage << "," << variants[i].coverage << std::endl;
	}
}

void liftoverAnnotationsToMorphs(const std::string& basePath, const std::vector<MorphConsensus>& morphConsensuses, const std::string& annotationFasta, const std::string& annotationGff3, const std::string& tmpPath, const std::string& liftoffPath)
{
	{
		std::ofstream typefile { tmpPath + "/liftoff_types.txt" };
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
		std::string tmpfilepath { tmpPath + "/tmp_seq.fa" };
		{
			std::ofstream tmpfile { tmpfilepath };
			tmpfile << ">" << morphConsensuses[i].name << std::endl;
			tmpfile << morphConsensuses[i].sequence << std::endl;
		}
		std::string command = liftoffPath + " -f " + tmpPath + "/liftoff_types.txt -g " + annotationGff3 + " -o " + tmpPath + "/tmp-morph-annotations-part" + std::to_string(i) + ".gff3 -u "+ tmpPath + "/morph-unmapped_features" + std::to_string(i) + ".txt -dir " + tmpPath + "/liftoff_intermediate_files/ " + tmpfilepath + " " + annotationFasta + " 1> " + tmpPath + "/liftoff_morphs_stdout" + std::to_string(i) + ".txt 2> " + tmpPath + "/liftoff_morphs_stderr" + std::to_string(i) + ".txt";
		std::cerr << "running liftoff with command:" << std::endl;
		std::cerr << command << std::endl;
		int result = system(command.c_str());
		if (result != 0)
		{
			std::cerr << "liftoff did not run successfully" << std::endl;
			std::abort();
		}
		command = "rm -r " + tmpPath + "/liftoff_intermediate_files/ " + tmpfilepath + ".fai " + tmpfilepath + ".mmi";
		std::cerr << command << std::endl;
		result = system(command.c_str());
	}
	std::cerr << "combining liftoff results" << std::endl;
	std::string outputFile = basePath + "/morph-annotations.gff3";
	std::string command = "echo \"##gff-version 3\" > " + outputFile + "\ncat " + tmpPath + "/tmp-morph-annotations-part*.gff3 | grep -v '#' >> " + outputFile;
	std::cerr << command << std::endl;
	int result = system(command.c_str());
	if (result != 0)
	{
		std::cerr << "failed to combine liftoff results" << std::endl;
		std::abort();
	}
}

void liftoverAnnotationsToConsensus(const std::string& basePath, const std::string& consensusPath, const std::string& annotationFasta, const std::string& annotationGff3, const std::string& tmpPath, const std::string& liftoffPath)
{
	{
		std::ofstream typefile { tmpPath + "/liftoff_types.txt" };
		typefile << "rRNA" << std::endl;
		typefile << "misc_RNA" << std::endl;
		typefile << "repeat_region" << std::endl;
		typefile << "gene" << std::endl;
		typefile << "transcript" << std::endl;
		typefile << "exon" << std::endl;
		typefile << "pseudogene" << std::endl;
		typefile << "tandem_repeat" << std::endl;
	}
	std::string command = liftoffPath + " -f " + tmpPath + "/liftoff_types.txt -g " + annotationGff3 + " -o " + basePath + "/consensus-annotation.gff3 -u " + tmpPath + "/consensus-unmapped_features.txt -dir " + tmpPath + "/liftoff_intermediate_files/ " + consensusPath + " " + annotationFasta + " 1> " + tmpPath + "/liftoff_consensus_stdout.txt 2> " + tmpPath + "/liftoff_consensus_stderr.txt";
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
		if (i == 0)
		{
			j -= 1;
			continue;
		}
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

std::tuple<std::unordered_set<size_t>, std::unordered_set<size_t>, std::unordered_map<Node, size_t>, std::unordered_map<Node, size_t>> getBorderNodes(const Path& heavyPath, const GfaGraph& graph)
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
	std::unordered_set<size_t> anchorNodes;
	std::unordered_map<Node, size_t> pathStartClip;
	std::unordered_map<Node, size_t> pathEndClip;
	for (size_t nodeid = 0; nodeid < graph.nodeSeqs.size(); nodeid++)
	{
		if (graph.edges.count(Node { nodeid, true }) == 0)
		{
			pathStartClip[Node { nodeid, false }] = 0;
			pathEndClip[Node { nodeid, true }] = graph.nodeSeqs[nodeid].size();
			anchorNodes.emplace(nodeid);
			borderNodes.emplace(nodeid);
		}
		if (graph.edges.count(Node { nodeid, false }) == 0)
		{
			pathStartClip[Node { nodeid, true }] = 0;
			pathEndClip[Node { nodeid, false }] = graph.nodeSeqs[nodeid].size();
			anchorNodes.emplace(nodeid);
			borderNodes.emplace(nodeid);
		}
	}
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
	return std::make_tuple(borderNodes, anchorNodes, pathStartClip, pathEndClip);
}

std::vector<OntLoop> extractLoopSequences(const std::vector<ReadPath>& correctedPaths, const Path& heavyPath, const size_t minLength, const GfaGraph& graph, const std::unordered_set<size_t>& borderNodes, const std::unordered_set<size_t>& anchorNodes, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
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
			bool anchored = anchorNodes.count(read.path[lastBreak].id()) == 1 || anchorNodes.count(read.path[i].id()) == 1;
			if (len > clip && (len - clip >= minLength || anchored))
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

size_t getEditDistanceWfa(const std::string& left, const std::string& right)
{
	std::vector<size_t> offsets;
	std::vector<size_t> nextOffsets;
	offsets.resize(1, 0);
	const size_t leftSize = left.size();
	const size_t rightSize = right.size();
	offsets[0] = getMatchLength(left, right, 0, 0);
	if (offsets[0] == leftSize) return rightSize - leftSize;
	if (offsets[0] == rightSize) return leftSize - rightSize;
	size_t zeroPos = 1;
	for (size_t score = 0; score < std::max(left.size(), right.size()); score++)
	{
		nextOffsets.resize(offsets.size()+2);
		for (size_t i = 0; i < nextOffsets.size(); i++)
		{
			nextOffsets[i] = 0;
			if (i < zeroPos && zeroPos-i >= left.size()) continue;
			if (i > zeroPos && i-zeroPos >= right.size()) continue;
			if (i >= 1 && i < nextOffsets.size()-1) nextOffsets[i] = offsets[i-1]+1;
			if (i >= 2) nextOffsets[i] = std::max(nextOffsets[i], offsets[i-2]);
			if (i < nextOffsets.size()-2) nextOffsets[i] = std::max(nextOffsets[i], offsets[i]+1);
			assert(nextOffsets[i] + i >= zeroPos);
			nextOffsets[i] = std::min(nextOffsets[i], rightSize - i + zeroPos);
			nextOffsets[i] = std::min(nextOffsets[i], leftSize);
			assert(nextOffsets[i] + i >= zeroPos);
			assert(nextOffsets[i]+i-zeroPos <= rightSize);
			assert(nextOffsets[i] <= leftSize);
			nextOffsets[i] += getMatchLength(left, right, nextOffsets[i], nextOffsets[i]+i-zeroPos);
			assert(nextOffsets[i]+i-zeroPos <= rightSize);
			assert(nextOffsets[i] <= leftSize);
			if (nextOffsets[i] == leftSize && nextOffsets[i]+i-zeroPos == rightSize) return score;
			assert(nextOffsets[i] + i >= zeroPos);
		}
		std::swap(offsets, nextOffsets);
		zeroPos += 1;
		assert(score < leftSize+rightSize);
	}
	assert(false);
	return std::max(left.size(), right.size())+1;
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

template <typename T>
T find(phmap::flat_hash_map<T, T>& parent, T key)
{
	while (parent.at(key) != parent.at(parent.at(key))) parent[key] = parent[parent[key]];
	return parent[key];
}

template <typename T>
void merge(phmap::flat_hash_map<T, T>& parent, T left, T right)
{
	left = find(parent, left);
	right = find(parent, right);
	assert(parent.at(left) == left);
	assert(parent.at(right) == right);
	parent[right] = left;
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

template <typename T>
T find(std::map<T, T>& parent, const T& key)
{
	while (parent.at(key) != parent.at(parent.at(key)))
	{
		parent[key] = parent.at(parent.at(key));
	}
	return parent.at(key);
}

template <typename T>
void merge(std::map<T, T>& parent, const T& left, const T& right)
{
	auto leftp = find(parent, left);
	auto rightp = find(parent, right);
	parent[rightp] = leftp;
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
		size_t countInClusters = 0;
		for (size_t j = 0; j < partialResult.size(); j++)
		{
			std::cerr << " " << partialResult[j].size();
			countInClusters += partialResult[j].size();
		}
		if (partialResult.size() < 2 || countInClusters < clusters[i].size()*0.9)
		{
			std::cerr << ", don't split" << std::endl;
			result.emplace_back(clusters[i]);
			continue;
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

std::vector<std::tuple<size_t, size_t, std::string, size_t>> getEdits(const std::string_view& refSequence, const std::string_view& querySequence, size_t maxEdits)
{
	if (refSequence == querySequence) return std::vector<std::tuple<size_t, size_t, std::string, size_t>>{};
	// special case for very small differences
	{
		size_t countEqualAtStart = 0;
		size_t countEqualAtEnd = 0;
		while (countEqualAtStart < refSequence.size() && countEqualAtStart < querySequence.size() && refSequence[countEqualAtStart] == querySequence[countEqualAtStart]) countEqualAtStart += 1;
		while (countEqualAtEnd < refSequence.size() && countEqualAtEnd < querySequence.size() && refSequence[refSequence.size()-1-countEqualAtEnd] == querySequence[querySequence.size()-1-countEqualAtEnd]) countEqualAtEnd += 1;
		if (refSequence.size() == querySequence.size() && countEqualAtEnd+countEqualAtStart+1 == refSequence.size())
		{
			// just a single SNP
			std::vector<std::tuple<size_t, size_t, std::string, size_t>> result;
			result.emplace_back(countEqualAtStart, countEqualAtStart+1, querySequence.substr(countEqualAtStart, 1), countEqualAtStart);
			return result;
		}
		if (querySequence.size() > refSequence.size() && countEqualAtEnd+countEqualAtStart == refSequence.size())
		{
			// just a single insertion
			std::vector<std::tuple<size_t, size_t, std::string, size_t>> result;
			result.emplace_back(countEqualAtStart, countEqualAtStart, querySequence.substr(countEqualAtStart, querySequence.size() - countEqualAtStart - countEqualAtEnd), countEqualAtStart);
			return result;
		}
		if (querySequence.size() < refSequence.size() && countEqualAtEnd+countEqualAtStart == querySequence.size())
		{
			// just a single deletion
			std::vector<std::tuple<size_t, size_t, std::string, size_t>> result;
			result.emplace_back(countEqualAtStart, refSequence.size()-countEqualAtEnd, "", countEqualAtStart);
			return result;
		}
	}
	// reference is horizontal, first bp left, last bp right (one bp one column)
	// query is vertical, first bp top, last bp bottom (one bp one row)
	// edits moves diagonally down-right
	// j moves horizontally, plus right minus left
	// insertion is extra sequence on query (vertical backtrace)
	// deletion is missing sequence on query (horizontal backtrace)
	// wfa matrix pair coord first is edits (row), second is j (column)
	// wfa matrix is triangular, missing parts are implied: edit 0 j 0 is above edit 1 j 1 and up-right from edit 1 j 0
	// real matrix pair coord first is ref pos (column), second is query pos (row)
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
	if (backtraceWfaPosition == uninitialized) return std::vector<std::tuple<size_t, size_t, std::string, size_t>> {};
	assert(backtraceWfaPosition != uninitialized);
	std::vector<std::tuple<size_t, size_t, std::string, size_t>> result;
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
					result.emplace_back(newMatchRealPosition.first, lastMatchRealPosition.first, querySequence.substr(newMatchRealPosition.second, lastMatchRealPosition.second - newMatchRealPosition.second), newMatchRealPosition.second);
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
	std::pair<size_t, size_t> newMatchRealPosition = Wfa::getRealMatrixPosition(backtraceWfaPosition, matchMatrix[backtraceWfaPosition.first][backtraceWfaPosition.second].second);
	if (lastMatchRealPosition != uninitialized)
	{
		assert(newMatchRealPosition.first <= lastMatchRealPosition.first);
		assert(newMatchRealPosition.second <= lastMatchRealPosition.second);
		result.emplace_back(newMatchRealPosition.first, lastMatchRealPosition.first, querySequence.substr(newMatchRealPosition.second, lastMatchRealPosition.second - newMatchRealPosition.second), newMatchRealPosition.second);
	}
	assert(result.size() >= 1);
	std::reverse(result.begin(), result.end());
	return result;
}

phmap::flat_hash_map<uint64_t, size_t> getRefKmers(const std::string_view& ref, const size_t k)
{
	assert(k <= 31);
	const uint64_t mask = (1ull << (2ull * k)) - 1;
	phmap::flat_hash_map<uint64_t, size_t> result;
	assert(ref.size() >= k);
	uint64_t kmer = 0;
	for (size_t i = 0; i < k; i++)
	{
		kmer <<= 2;
		switch(ref[i])
		{
			case 'A':
			case 'a':
				kmer += 0;
				break;
			case 'C':
			case 'c':
				kmer += 1;
				break;
			case 'G':
			case 'g':
				kmer += 2;
				break;
			case 'T':
			case 't':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
	}
	result[kmer] = 0;
	phmap::flat_hash_set<uint64_t> eraseThese;
	for (size_t i = k; i < ref.size(); i++)
	{
		kmer <<= 2;
		kmer &= mask;
		switch(ref[i])
		{
			case 'A':
			case 'a':
				kmer += 0;
				break;
			case 'C':
			case 'c':
				kmer += 1;
				break;
			case 'G':
			case 'g':
				kmer += 2;
				break;
			case 'T':
			case 't':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
		if (result.count(kmer) == 1) eraseThese.emplace(kmer);
		result[kmer] = i-k+1;
	}
	for (uint64_t kmer : eraseThese)
	{
		assert(result.count(kmer) == 1);
		result.erase(kmer);
	}
	return result;
}

void removeRepeatKmers(std::vector<std::pair<uint64_t, uint64_t>>& matches)
{
	phmap::flat_hash_set<size_t> foundOnce;
	phmap::flat_hash_set<size_t> foundTwice;
	for (const auto& pair : matches)
	{
		if (foundOnce.count(pair.first) == 1) foundTwice.emplace(pair.first);
		foundOnce.emplace(pair.first);
	}
	if (foundTwice.size() == 0) return;
	for (size_t i = matches.size()-1; i < matches.size(); i--)
	{
		if (foundTwice.count(matches[i].first) == 0) continue;
		std::swap(matches[i], matches.back());
		matches.pop_back();
	}
	std::sort(matches.begin(), matches.end());
}

void removeRepeatKmersEitherRepeat(std::vector<std::pair<uint64_t, uint64_t>>& matches)
{
	phmap::flat_hash_set<size_t> foundOnceFirst;
	phmap::flat_hash_set<size_t> foundTwiceFirst;
	phmap::flat_hash_set<size_t> foundOnceSecond;
	phmap::flat_hash_set<size_t> foundTwiceSecond;
	for (const auto& pair : matches)
	{
		if (foundOnceFirst.count(pair.first) == 1) foundTwiceFirst.emplace(pair.first);
		foundOnceFirst.emplace(pair.first);
		if (foundOnceSecond.count(pair.second) == 1) foundTwiceSecond.emplace(pair.second);
		foundOnceSecond.emplace(pair.second);
	}
	if (foundTwiceFirst.size() == 0 && foundTwiceSecond.size() == 0) return;
	for (size_t i = matches.size()-1; i < matches.size(); i--)
	{
		if (foundTwiceFirst.count(matches[i].first) == 0 && foundTwiceSecond.count(matches[i].second) == 0) continue;
		std::swap(matches[i], matches.back());
		matches.pop_back();
	}
	std::sort(matches.begin(), matches.end());
}

std::vector<std::pair<size_t, size_t>> getIncreasingChain(const std::vector<std::pair<size_t, size_t>>& kmerMatches)
{
	if (kmerMatches.size() == 0) return kmerMatches;
	std::vector<std::vector<size_t>> paretoFrontier;
	paretoFrontier.emplace_back();
	paretoFrontier.emplace_back();
	paretoFrontier[1].emplace_back(0);
	for (size_t i = 1; i < kmerMatches.size(); i++)
	{
		bool found = false;
		for (size_t length = paretoFrontier.size()-1; length < paretoFrontier.size(); length++)
		{
			for (size_t j : paretoFrontier[length])
			{
				assert(kmerMatches[i].first > kmerMatches[j].first);
				assert(kmerMatches[i].second != kmerMatches[j].second);
				if (kmerMatches[i].second < kmerMatches[j].second) continue;
				found = true;
				break;
			}
			if (found)
			{
				assert(length+1 <= paretoFrontier.size());
				if (length+1 == paretoFrontier.size()) paretoFrontier.emplace_back();
				bool thisIsInParetoFrontier = true;
				for (size_t j : paretoFrontier[length+1])
				{
					assert(kmerMatches[i].first > kmerMatches[j].first);
					assert(kmerMatches[i].second != kmerMatches[j].second);
					if (kmerMatches[i].second < kmerMatches[j].second) continue;
					thisIsInParetoFrontier = false;
					break;
				}
				if (thisIsInParetoFrontier) paretoFrontier[length+1].emplace_back(i);
				break;
			}
		}
		if (!found)
		{
			bool thisIsInParetoFrontier = true;
			for (size_t j : paretoFrontier[1])
			{
				assert(kmerMatches[i].first > kmerMatches[j].first);
				assert(kmerMatches[i].second != kmerMatches[j].second);
				if (kmerMatches[i].second < kmerMatches[j].second) continue;
				thisIsInParetoFrontier = false;
				break;
			}
			if (thisIsInParetoFrontier) paretoFrontier[1].emplace_back(i);
		}
	}
	assert(paretoFrontier.back().size() >= 1);
	size_t maxIndex = paretoFrontier.back()[0];
	std::vector<std::pair<size_t, size_t>> result;
	result.emplace_back(kmerMatches[maxIndex]);
	for (size_t length = paretoFrontier.size()-2; length > 0; length--)
	{
		bool found = false;
		for (size_t j : paretoFrontier[length])
		{
			if (kmerMatches[j].first > kmerMatches[maxIndex].first) continue;
			if (kmerMatches[j].second > kmerMatches[maxIndex].second) continue;
			found = true;
			result.emplace_back(kmerMatches[j]);
			maxIndex = j;
			break;
		}
		assert(found);
	}
	std::reverse(result.begin(), result.end());
	return result;
}

void removeNonDiagonalKmers(std::vector<std::pair<size_t, size_t>>& kmerMatches)
{
	std::sort(kmerMatches.begin(), kmerMatches.end(), [](auto left, auto right){ return (int)left.first - (int)left.second < (int)right.first - (int)right.second; });
	size_t bestClusterSize = 0;
	int bestClusterStart = 0;
	int bestClusterEnd = 0;
	size_t currentClusterStart = 0;
	for (size_t i = 1; i <= kmerMatches.size(); i++)
	{
		if (i == kmerMatches.size() || (int)kmerMatches[i].first - (int)kmerMatches[i].second > (int)kmerMatches[i-1].first - (int)kmerMatches[i-1].second + 100)
		{
			size_t currentClusterSize = i - currentClusterStart;
			if (currentClusterSize > bestClusterSize)
			{
				bestClusterSize = currentClusterSize;
				bestClusterStart = (int)kmerMatches[currentClusterStart].first - (int)kmerMatches[currentClusterStart].second;
				bestClusterEnd = (int)kmerMatches[i-1].first - (int)kmerMatches[i-1].second;
			}
			currentClusterSize = 0;
			currentClusterStart = i;
		}
	}
	for (size_t i = kmerMatches.size()-1; i < kmerMatches.size(); i--)
	{
		int diagonal = kmerMatches[i].first - kmerMatches[i].second;
		if (diagonal >= bestClusterStart && diagonal <= bestClusterEnd) continue;
		std::swap(kmerMatches[i], kmerMatches.back());
		kmerMatches.pop_back();
	}
	std::sort(kmerMatches.begin(), kmerMatches.end());
}

std::vector<std::pair<size_t, size_t>> getKmerAnchors(const std::string_view& ref, const phmap::flat_hash_map<uint64_t, size_t>& refKmers, const std::string_view& query, const size_t k)
{
	assert(k <= 31);
	const uint64_t mask = (1ull << (2ull * k)) - 1;
	assert(ref.size() >= k);
	uint64_t kmer = 0;
	std::vector<std::pair<size_t, size_t>> kmerMatches;
	for (size_t i = 0; i < k; i++)
	{
		kmer <<= 2;
		switch(query[i])
		{
			case 'A':
			case 'a':
				kmer += 0;
				break;
			case 'C':
			case 'c':
				kmer += 1;
				break;
			case 'G':
			case 'g':
				kmer += 2;
				break;
			case 'T':
			case 't':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
	}
	if (refKmers.count(kmer) == 1) kmerMatches.emplace_back(refKmers.at(kmer), 0);
	for (size_t i = k; i < query.size(); i++)
	{
		kmer <<= 2;
		kmer &= mask;
		switch(query[i])
		{
			case 'A':
			case 'a':
				kmer += 0;
				break;
			case 'C':
			case 'c':
				kmer += 1;
				break;
			case 'G':
			case 'g':
				kmer += 2;
				break;
			case 'T':
			case 't':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
		if (refKmers.count(kmer) == 1) kmerMatches.emplace_back(refKmers.at(kmer), i-k+1);
	}
	for (auto pair : kmerMatches)
	{
		assert(ref.substr(pair.first, k) == query.substr(pair.second, k));
	}
	std::sort(kmerMatches.begin(), kmerMatches.end());
	removeNonDiagonalKmers(kmerMatches);
	removeRepeatKmers(kmerMatches);
	std::vector<std::pair<size_t, size_t>> chain = getIncreasingChain(kmerMatches);
	return chain;
}

std::vector<std::tuple<size_t, size_t, std::string, size_t>> getEditsRec(const std::string_view& refSeq, const std::string_view& querySeq, const size_t k)
{
	assert(refSeq.size() >= k);
	assert(querySeq.size() >= k);
	assert(refSeq.substr(0, k) == querySeq.substr(0, k));
	assert(refSeq.substr(refSeq.size()-k) == querySeq.substr(querySeq.size()-k));
	if (k <= 11 || k > 31)
	{
		return getEdits(refSeq, querySeq, refSeq.size()*2);
	}
	auto refKmers = getRefKmers(refSeq, k);
	auto anchors = getKmerAnchors(refSeq, refKmers, querySeq, k);
	if (anchors.size() >= 1 && anchors.back().first+k == refSeq.size() && anchors.back().second+k == querySeq.size())
	{
		anchors.pop_back();
	}
	if (anchors.size() >= 1 && anchors[0].first == 0 && anchors[0].second == 0)
	{
		anchors.erase(anchors.begin());
	}
	if (anchors.size() == 0)
	{
		return getEdits(refSeq, querySeq, refSeq.size()*2);
	}
	std::vector<std::tuple<size_t, size_t, std::string, size_t>> result;
	size_t lastRefPos = 0;
	size_t lastQueryPos = 0;
	for (size_t i = 0; i < anchors.size(); i++)
	{
		if (anchors[i].first - lastRefPos == anchors[i].second - lastQueryPos && anchors[i].first - lastRefPos < k)
		{
			lastRefPos = anchors[i].first;
			lastQueryPos = anchors[i].second;
			continue;
		}
		auto partialResult = getEditsRec(refSeq.substr(lastRefPos, anchors[i].first - lastRefPos + k), querySeq.substr(lastQueryPos, anchors[i].second - lastQueryPos + k), k-10);
		for (size_t j = 0; j < partialResult.size(); j++)
		{
			result.emplace_back();
			std::swap(result.back(), partialResult[j]);
			std::get<0>(result.back()) += lastRefPos;
			std::get<1>(result.back()) += lastRefPos;
			std::get<3>(result.back()) += lastQueryPos;
		}
		lastRefPos = anchors[i].first;
		lastQueryPos = anchors[i].second;
	}
	if (refSeq.size()-lastRefPos == querySeq.size()-lastQueryPos && refSeq.size()-lastRefPos < k)
	{
	}
	else
	{
		auto partialResult = getEditsRec(refSeq.substr(lastRefPos), querySeq.substr(lastQueryPos), k-10);
		for (size_t j = 0; j < partialResult.size(); j++)
		{
			result.emplace_back();
			std::swap(result.back(), partialResult[j]);
			std::get<0>(result.back()) += lastRefPos;
			std::get<1>(result.back()) += lastRefPos;
			std::get<3>(result.back()) += lastQueryPos;
		}
	}
	return result;
}

template <typename F>
void iterateEdits(const std::string& rawConsensus, const phmap::flat_hash_map<uint64_t, size_t>& refKmers, const std::string& sequence, const size_t k, const size_t threadIndex, const size_t readIndex, F callback)
{
	auto anchors = getKmerAnchors(std::string_view { rawConsensus }, refKmers, std::string_view { sequence }, k);
	if (anchors.size() == 0) return;
	for (size_t i = 1; i < anchors.size(); i++)
	{
		assert(anchors[i-1].first < anchors[i].first);
		assert(anchors[i-1].second < anchors[i].second);
		assert(anchors[i].first+k <= rawConsensus.size());
		assert(anchors[i].second+k <= sequence.size());
		assert(rawConsensus.substr(anchors[i].first, k) == sequence.substr(anchors[i].second, k));
		assert(rawConsensus.substr(anchors[i-1].first, k) == sequence.substr(anchors[i-1].second, k));
		if (anchors[i].first - anchors[i-1].first == anchors[i].second - anchors[i-1].second && anchors[i].first - anchors[i-1].first < k) continue;
		std::string_view refSubstr { rawConsensus.data()+anchors[i-1].first, anchors[i].first - anchors[i-1].first + k };
		std::string_view querySubstr { sequence.data()+anchors[i-1].second, anchors[i].second - anchors[i-1].second + k };
		auto edits = getEditsRec(refSubstr, querySubstr, 25);
		for (auto edit : edits)
		{
			assert(std::get<1>(edit) >= std::get<0>(edit));
			assert(std::get<1>(edit) < refSubstr.size());
			std::get<0>(edit) += anchors[i-1].first;
			std::get<1>(edit) += anchors[i-1].first;
			assert(std::get<1>(edit) < rawConsensus.size());
			std::get<3>(edit) += anchors[i-1].second;
			callback(threadIndex, readIndex, anchors[0].first, anchors.back().first+k, anchors[0].second, anchors.back().second+k, edit);
		}
	}
}

template <typename F>
void iterateEdits(const std::string& rawConsensus, const std::vector<OntLoop>& rawLoops, const size_t numThreads, F callback)
{
	const size_t k = 31;
	auto refKmers = getRefKmers(std::string_view { rawConsensus }, k);
	std::atomic<size_t> nextIndex;
	nextIndex = 0;
	std::vector<std::thread> threads;
	for (size_t threadIndex = 0; threadIndex < numThreads; threadIndex++)
	{
		threads.emplace_back([&rawConsensus, &rawLoops, &nextIndex, &refKmers, threadIndex, callback]()
		{
			while (true)
			{
				size_t readIndex = nextIndex++;
				if (readIndex >= rawLoops.size()) break;
				iterateEdits(rawConsensus, refKmers, rawLoops[readIndex].rawSequence, k, threadIndex, readIndex, callback);
			}
		});
	}
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
}

template <typename F>
void iterateEdits(const std::string& rawConsensus, const std::vector<std::string>& rawLoops, const size_t numThreads, F callback)
{
	const size_t k = 31;
	auto refKmers = getRefKmers(std::string_view { rawConsensus }, k);
	std::atomic<size_t> nextIndex;
	nextIndex = 0;
	std::vector<std::thread> threads;
	for (size_t threadIndex = 0; threadIndex < numThreads; threadIndex++)
	{
		threads.emplace_back([&rawConsensus, &rawLoops, &nextIndex, &refKmers, threadIndex, callback]()
		{
			while (true)
			{
				size_t readIndex = nextIndex++;
				if (readIndex >= rawLoops.size()) break;
				iterateEdits(rawConsensus, refKmers, rawLoops[readIndex], k, threadIndex, readIndex, callback);
			}
		});
	}
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
}

std::string getPolishedSequence(const std::string& rawConsensus, const std::vector<std::string>& loopSequences, const size_t numThreads)
{
	std::vector<std::map<std::tuple<size_t, size_t, std::string>, size_t>> editCountsPerThread;
	editCountsPerThread.resize(numThreads);
	iterateEdits(rawConsensus, loopSequences, numThreads, [&editCountsPerThread](const size_t threadIndex, const size_t readIndex, const size_t firstMatchRef, const size_t lastMatchRef, const size_t firstMatchRead, const size_t lastMatchRead, const std::tuple<size_t, size_t, std::string, size_t>& edit)
	{
		assert(threadIndex < editCountsPerThread.size());
		auto filterEdit = std::make_tuple(std::get<0>(edit), std::get<1>(edit), std::get<2>(edit));
		editCountsPerThread[threadIndex][filterEdit] += 1;
	});
	std::map<std::tuple<size_t, size_t, std::string>, size_t> editCounts;
	for (size_t i = 0; i < editCountsPerThread.size(); i++)
	{
		for (const auto& pair : editCountsPerThread[i])
		{
			editCounts[pair.first] += pair.second;
		}
	}
	std::vector<std::tuple<size_t, size_t, std::string>> pickedEdits;
	for (const auto& pair : editCounts)
	{
		if (pair.second*1.5 <= loopSequences.size()) continue; // only pick edits supported by strictly >2/3 reads
		pickedEdits.emplace_back(pair.first);
	}
	std::sort(pickedEdits.begin(), pickedEdits.end());
	std::cerr << "picked " << pickedEdits.size() << " edits" << std::endl;
	std::string result;
	size_t lastMatch = 0;
	for (size_t i = 0; i < pickedEdits.size(); i++)
	{
		if (!(i == 0 || std::get<0>(pickedEdits[i]) >= std::get<1>(pickedEdits[i-1])))
		{
			std::cerr << pickedEdits.size() << " " << i << " " << loopSequences.size() << std::endl;
			for (size_t j = 0; j < pickedEdits.size(); j++)
			{
				std::cerr << j << " " << std::get<0>(pickedEdits[j]) << " " << std::get<1>(pickedEdits[j]) << " " << std::get<2>(pickedEdits[j]) << " " << editCounts.at(pickedEdits[j]) << std::endl;
			}
		}
		assert(i == 0 || std::get<0>(pickedEdits[i]) >= std::get<1>(pickedEdits[i-1]));
		assert(std::get<0>(pickedEdits[i]) >= lastMatch);
		result += rawConsensus.substr(lastMatch, std::get<0>(pickedEdits[i]) - lastMatch);
		result += std::get<2>(pickedEdits[i]);
		lastMatch = std::get<1>(pickedEdits[i]);
	}
	result += rawConsensus.substr(lastMatch);
	return result;
}

std::vector<MorphConsensus> getMorphConsensuses(const std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::string& namePrefix)
{
	std::vector<MorphConsensus> result;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		result.emplace_back();
		result.back().path = getConsensusPath(clusters[i], graph);
		result.back().preCorrectionSequence = getConsensusSequence(result.back().path, graph, pathStartClip, pathEndClip);
		result.back().sequence = result.back().preCorrectionSequence;
		result.back().coverage = clusters[i].size();
		result.back().ontLoops = clusters[i];
		result.back().name = "morphconsensus" + std::to_string(i) + "_coverage" + std::to_string(result.back().coverage);
		if (namePrefix != "") result.back().name = namePrefix + "_" + result.back().name;
	}
	return result;
}

std::string polishByDBG(const std::string& refSequence, const std::vector<std::string>& ontSequences, const std::string MBGPath, const std::string tmpPrefix, const size_t numThreads)
{
	std::string tmpReadsFile = tmpPrefix + "reads.fa";
	std::string tmpPathsFile = tmpPrefix + "paths.gaf";
	std::string tmpGraphFile = tmpPrefix + "graph.gfa";
	{
		std::ofstream file { tmpReadsFile };
		// two copies of consensus guarantees a contiguous path through the graph exists with min coverage 2 no matter any errors in reads, also slightly biases ambiguous sites towards no correction ie hifi consensus
		file << ">consensus" << std::endl;
		file << refSequence << std::endl;
		file << ">consensus_dupe" << std::endl;
		file << refSequence << std::endl;
		for (size_t i = 0; i < ontSequences.size(); i++)
		{
			file << ">sequence" << i << std::endl;
			file << ontSequences[i] << std::endl;
		}
	}
	std::string mbgCommand = MBGPath + " -i " + tmpReadsFile + " -o " + tmpGraphFile + " --output-sequence-paths " + tmpPathsFile + " -k 31 -a 2 -u 2 --error-masking=msat --no-multiplex-cleaning -r 18446744073709551615 -R 18446744073709551615 --only-local-resolve --do-unsafe-guesswork-resolutions 1> " + tmpPrefix + "stdout.txt 2> " + tmpPrefix + "stderr.txt";
	std::cerr << "MBG command:" << std::endl;
	std::cerr << mbgCommand << std::endl;
	int mbgresult = system(mbgCommand.c_str());
	if (mbgresult != 0)
	{
		std::cerr << "MBG did not run successfully" << std::endl;
		std::abort();
	}
	GfaGraph graph;
	graph.loadFromFile(tmpGraphFile);
	std::vector<ReadPath> readPaths = loadReadPaths(tmpPathsFile, graph);
	std::vector<Node> refPath;
	phmap::flat_hash_map<Node, size_t> nodeIndexInReference;
	phmap::flat_hash_set<Node> eraseThese;
	size_t startClip = 0;
	size_t endClip = 0;
	size_t startPlus = 0;
	size_t endPlus = 0;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		if (readPaths[i].readName != "consensus") continue;
		assert(readPaths[i].readStart < readPaths[i].readEnd);
		assert(readPaths[i].readEnd <= refSequence.size());
		assert(refPath.size() == 0);
		startPlus = readPaths[i].readStart;
		endPlus = refSequence.size() - readPaths[i].readEnd;
		startClip = readPaths[i].pathStartClip;
		endClip = readPaths[i].pathEndClip;
		refPath = readPaths[i].path;
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (nodeIndexInReference.count(readPaths[i].path[j]) == 1)
			{
				nodeIndexInReference[readPaths[i].path[j]] = std::numeric_limits<size_t>::max();
				eraseThese.emplace(readPaths[i].path[j]);
			}
			else
			{
				nodeIndexInReference[readPaths[i].path[j]] = j;
			}
		}
	}
	for (Node node : eraseThese)
	{
		assert(nodeIndexInReference.count(node) == 1);
		nodeIndexInReference.erase(node);
	}
	assert(refPath.size() >= 1);
	std::map<std::vector<Node>, size_t> alternatePathCounts;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		if (readPaths[i].readName == "consensus") continue;
		if (readPaths[i].readName == "consensus_dupe") continue;
		size_t lastRefMatch = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (nodeIndexInReference.count(readPaths[i].path[j]) == 0) continue;
			if (lastRefMatch == std::numeric_limits<size_t>::max())
			{
				lastRefMatch = j;
				continue;
			}
			std::vector<Node> altPath { readPaths[i].path.begin()+lastRefMatch, readPaths[i].path.begin()+j+1 };
			alternatePathCounts[altPath] += 1;
			lastRefMatch = j;
		}
		if (lastRefMatch == std::numeric_limits<size_t>::max())
		{
			// read is bw
			for (size_t j = readPaths[i].path.size()-1; j < readPaths[i].path.size(); j--)
			{
				if (nodeIndexInReference.count(reverse(readPaths[i].path[j])) == 0) continue;
				if (lastRefMatch == std::numeric_limits<size_t>::max())
				{
					lastRefMatch = j;
					continue;
				}
				std::vector<Node> altPath { readPaths[i].path.begin()+j, readPaths[i].path.begin()+lastRefMatch+1 };
				altPath = reverse(altPath);
				alternatePathCounts[altPath] += 1;
				lastRefMatch = j;
			}
		}
	}
	std::vector<std::vector<std::tuple<size_t, std::vector<Node>, size_t>>> inNeighbors;
	inNeighbors.resize(refPath.size());
	for (size_t i = 1; i < refPath.size(); i++)
	{
		std::vector<Node> thisPath { refPath[i-1], refPath[i] };
		size_t coverage = 1;
		if (alternatePathCounts.count(thisPath) == 1)
		{
			coverage += alternatePathCounts.at(thisPath);
		}
		inNeighbors[i].emplace_back(i-1, thisPath, coverage);
	}
	for (const auto& pair : alternatePathCounts)
	{
		assert(pair.first.size() >= 2);
		assert(nodeIndexInReference.count(pair.first[0]) == 1);
		assert(nodeIndexInReference.count(pair.first.back()) == 1);
		if (nodeIndexInReference.at(pair.first[0]) >= nodeIndexInReference.at(pair.first.back())) continue;
		inNeighbors[nodeIndexInReference.at(pair.first.back())].emplace_back(nodeIndexInReference.at(pair.first[0]), pair.first, pair.second);
	}
	assert(inNeighbors[0].size() == 0);
	std::vector<std::tuple<size_t, std::vector<Node>, size_t>> pickedInNeighbor;
	pickedInNeighbor.emplace_back(0, std::vector<Node>{}, std::numeric_limits<size_t>::max());
	for (size_t i = 1; i < refPath.size(); i++)
	{
		assert(inNeighbors[i].size() >= 1);
		pickedInNeighbor.emplace_back(0, std::vector<Node>{}, 0);
		for (const auto& neighbor : inNeighbors[i])
		{
			size_t coverageHere = std::min(std::get<2>(neighbor), std::get<2>(pickedInNeighbor[std::get<0>(neighbor)]));
			if (coverageHere <= std::get<2>(pickedInNeighbor.back())) continue;
			pickedInNeighbor.back() = neighbor;
			assert(coverageHere <= std::get<2>(pickedInNeighbor.back()));
			std::get<2>(pickedInNeighbor.back()) = coverageHere;
		}
	}
	std::vector<Node> pickedPath;
	size_t index = refPath.size()-1;
	while (index != 0)
	{
		assert(std::get<0>(pickedInNeighbor[index]) < index);
		std::vector<Node> addedHere = std::get<1>(pickedInNeighbor[index]);
		addedHere = reverse(addedHere);
		assert(addedHere.size() >= 2);
		assert(addedHere.back() == reverse(refPath[std::get<0>(pickedInNeighbor[index])]));
		pickedPath.insert(pickedPath.end(), addedHere.begin(), addedHere.end()-1);
		index = std::get<0>(pickedInNeighbor[index]);
	}
	pickedPath.emplace_back(reverse(refPath[0]));
	pickedPath = reverse(pickedPath);
	std::string result = getSequence(pickedPath, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
	result = refSequence.substr(0, startPlus) + result.substr(startClip, result.size() - startClip - endClip) + refSequence.substr(refSequence.size()-endPlus);
	std::cerr << "DBG try to polish from size " << refSequence.size() << " to " << result.size() << std::endl;
	if (result.size() > refSequence.size()+100) return refSequence;
	if (result.size()+100 < refSequence.size()) return refSequence;
	return result;
}

void polishMorphConsensuses(std::vector<MorphConsensus>& morphConsensuses, const std::vector<std::vector<OntLoop>>& ontLoopSequences, const std::string MBGPath, const std::string tmpPath, const size_t numThreads)
{
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		if (ontLoopSequences[i].size() < 3)
		{
			std::cerr << "skip polishing morph " << i << " with coverage " << ontLoopSequences[i].size() << std::endl;
			continue;
		}
		std::vector<std::string> seqs;
		for (size_t k = 0; k < ontLoopSequences[i].size(); k++)
		{
			seqs.emplace_back(ontLoopSequences[i][k].rawSequence);
		}
		morphConsensuses[i].sequence = polishByDBG(morphConsensuses[i].sequence, seqs, MBGPath, tmpPath + "/polish" + std::to_string(i), numThreads);
		for (size_t j = 0; j < 5; j++)
		{
			std::cerr << "polish morph " << i << " consensus round " << j << "/5" << std::endl;
			std::string newCorrection = getPolishedSequence(morphConsensuses[i].sequence, seqs, numThreads);
			if (newCorrection == morphConsensuses[i].sequence) break;
			morphConsensuses[i].sequence = newCorrection;
		}
	}
}

void writeMorphConsensuses(std::string outFile, std::string preconsensusOutputFile, const std::vector<MorphConsensus>& morphConsensuses)
{
	std::ofstream file { outFile };
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		file << ">" << morphConsensuses[i].name << std::endl;
		file << morphConsensuses[i].sequence << std::endl;
	}
	std::ofstream filePreCorrection { preconsensusOutputFile };
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		filePreCorrection << ">" << morphConsensuses[i].name << std::endl;
		filePreCorrection << morphConsensuses[i].preCorrectionSequence << std::endl;
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

void addNonrDNAAnchorsAsVariants(std::vector<Variant>& variants, const GfaGraph& graph, const Path& heavyPath, const std::vector<ReadPath>& readPaths)
{
	const size_t minAnchorLength = 5000;
	phmap::flat_hash_set<size_t> nodesInHeavyPath;
	for (auto node : heavyPath.nodes)
	{
		nodesInHeavyPath.emplace(node.id());
	}
	size_t separatorNode = heavyPath.nodes[0].id();
	phmap::flat_hash_set<Node> reachableFromSeparator;
	std::vector<Node> checkStack;
	checkStack.emplace_back(separatorNode, true);
	checkStack.emplace_back(separatorNode, false);
	while (checkStack.size() > 0)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (reachableFromSeparator.count(top) == 1) continue;
		reachableFromSeparator.insert(top);
		if (graph.edges.count(top) == 1)
		{
			for (auto edge : graph.edges.at(top))
			{
				checkStack.emplace_back(std::get<0>(edge));
			}
		}
	}
	std::map<std::vector<Node>, size_t> danglerCount;
	for (const ReadPath& path : readPaths)
	{
		size_t prefixDanglerLength = 0;
		for (size_t i = 0; i < path.path.size(); i++)
		{
			if (reachableFromSeparator.count(path.path[i]) == 0)
			{
				prefixDanglerLength += 1;
			}
			else
			{
				break;
			}
		}
		size_t firstHeavyPathIndex = std::numeric_limits<size_t>::max();
		size_t lastHeavyPathIndex = std::numeric_limits<size_t>::max();
		for (size_t i = 0; i < path.path.size(); i++)
		{
			if (nodesInHeavyPath.count(path.path[i].id()) == 1)
			{
				if (firstHeavyPathIndex == std::numeric_limits<size_t>::max()) firstHeavyPathIndex = i;
				lastHeavyPathIndex = i;
			}
		}
		if (prefixDanglerLength > 0 && prefixDanglerLength < path.path.size() && firstHeavyPathIndex != std::numeric_limits<size_t>::max())
		{
			assert(firstHeavyPathIndex >= prefixDanglerLength);
			std::vector<Node> danglerPath { path.path.begin(), path.path.begin() + prefixDanglerLength };
			size_t danglerPathLength = getPathLength(danglerPath, graph.nodeSeqs, graph.edges);
			if (danglerPathLength >= minAnchorLength)
			{
				danglerPath = reverse(danglerPath);
				while (danglerPath.size() >= 2 && getPathLength(std::vector<Node> { danglerPath.begin(), danglerPath.end()-1 }, graph.nodeSeqs, graph.edges) >= minAnchorLength) danglerPath.pop_back();
				std::vector<Node> pathFromHeavyPathUntilDangler = reverse(std::vector<Node> { path.path.begin() + prefixDanglerLength, path.path.begin() + firstHeavyPathIndex + 1});
				danglerPath.insert(danglerPath.begin(), pathFromHeavyPathUntilDangler.begin(), pathFromHeavyPathUntilDangler.end());
				assert(danglerPath.size() >= 2);
				assert(nodesInHeavyPath.count(danglerPath[0].id()) == 1);
				assert(nodesInHeavyPath.count(danglerPath[1].id()) == 0);
				danglerCount[danglerPath] += 1;
			}
		}
		size_t suffixDanglerLength = 0;
		for (size_t i = path.path.size()-1; i < path.path.size(); i--)
		{
			if (reachableFromSeparator.count(reverse(path.path[i])) == 0)
			{
				suffixDanglerLength += 1;
			}
			else
			{
				break;
			}
		}
		if (suffixDanglerLength > 0 && suffixDanglerLength < path.path.size() && lastHeavyPathIndex != std::numeric_limits<size_t>::max())
		{
			assert(lastHeavyPathIndex <= path.path.size() - suffixDanglerLength);
			std::vector<Node> danglerPath { path.path.end() - suffixDanglerLength, path.path.end() };
			size_t danglerPathLength = getPathLength(danglerPath, graph.nodeSeqs, graph.edges);
			if (danglerPathLength >= minAnchorLength)
			{
				while (danglerPath.size() >= 2 && getPathLength(std::vector<Node> { danglerPath.begin(), danglerPath.end()-1 }, graph.nodeSeqs, graph.edges) >= minAnchorLength) danglerPath.pop_back();
				std::vector<Node> pathFromHeavyPathUntilDangler = std::vector<Node> { path.path.begin() + lastHeavyPathIndex, path.path.end() - suffixDanglerLength };
				danglerPath.insert(danglerPath.begin(), pathFromHeavyPathUntilDangler.begin(), pathFromHeavyPathUntilDangler.end());
				assert(danglerPath.size() >= 2);
				assert(nodesInHeavyPath.count(danglerPath[0].id()) == 1);
				assert(nodesInHeavyPath.count(danglerPath[1].id()) == 0);
				danglerCount[danglerPath] += 1;
			}
		}
	}
	for (const auto& pair : danglerCount)
	{
		if (pair.second < 3) continue;
		assert(pair.first.size() >= 2);
		assert(nodesInHeavyPath.count(pair.first[0].id()) == 1);
		for (size_t j = 1; j < pair.first.size(); j++)
		{
			assert(nodesInHeavyPath.count(pair.first[j].id()) == 0);
		}
		assert(reachableFromSeparator.count(reverse(pair.first.back())) == 0);
		variants.emplace_back();
		variants.back().path = pair.first;
		variants.back().coverage = pair.second;
		variants.back().referencePath = std::vector<Node> { pair.first.begin(), pair.first.begin()+1 };
		variants.back().referenceCoverage = 0;
		variants.back().nonrDNAAnchorVariant = true;
	}
}

void processGraphAndWrite(const Path& heavyPath, const GfaGraph& graph, const std::string outputGraphName)
{
	const size_t minCoverage = 5;
	const size_t minAnchorLength = 10000;
	phmap::flat_hash_set<size_t> coveredNodes;
	std::vector<Node> checkStack;
	phmap::flat_hash_set<Node> reachableNodes;
	for (size_t i = 0; i < heavyPath.nodes.size(); i++)
	{
		coveredNodes.insert(heavyPath.nodes[i].id());
		checkStack.emplace_back(Node { heavyPath.nodes[i].id(), true });
		checkStack.emplace_back(Node { heavyPath.nodes[i].id(), false });
	}
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (graph.nodeCoverages[i] < minCoverage) continue;
		coveredNodes.insert(i);
	}
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (reachableNodes.count(top) == 1) continue;
		reachableNodes.insert(top);
		if (graph.edges.count(top) == 1)
		{
			for (auto edge : graph.edges.at(top))
			{
				auto target = std::get<0>(edge);
				if (coveredNodes.count(target.id()) == 0) continue;
				checkStack.emplace_back(target);
			}
		}
	}
	phmap::flat_hash_set<size_t> keptNodesInCycles;
	for (size_t node : coveredNodes)
	{
		if (reachableNodes.count(Node { node, true }) == 0) continue;
		if (reachableNodes.count(Node { node, false }) == 0) continue;
		keptNodesInCycles.insert(node);
	}
	phmap::flat_hash_map<Node, size_t> shortestDistanceToCyclicComponent;
	std::vector<std::pair<Node, size_t>> checkStack2;
	for (size_t node : keptNodesInCycles)
	{
		checkStack2.emplace_back(Node { node, true }, 0);
		checkStack2.emplace_back(Node { node, false }, 0);
	}
	while (checkStack2.size() >= 1)
	{
		auto top = checkStack2.back();
		checkStack2.pop_back();
		if (shortestDistanceToCyclicComponent.count(top.first) == 1)
		{
			assert(shortestDistanceToCyclicComponent.at(top.first) <= top.second);
			continue;
		}
		shortestDistanceToCyclicComponent[top.first] = top.second;
		bool addedAny = false;
		if (graph.edges.count(top.first) == 1)
		{
			for (auto edge : graph.edges.at(top.first))
			{
				if (keptNodesInCycles.count(std::get<0>(edge).id()) == 1) continue;
				checkStack2.emplace_back(std::get<0>(edge), top.second + graph.nodeSeqs[std::get<0>(edge).id()].size() - std::get<1>(edge));
				assert(checkStack2.back().second > top.second);
				addedAny = true;
			}
		}
		if (addedAny) std::sort(checkStack2.begin(), checkStack2.end(), [](auto left, auto right) { return left.second > right.second; });
	}
	phmap::flat_hash_set<size_t> keptNodesInAnchors;
	std::vector<Node> checkStack3;
	phmap::flat_hash_set<Node> nodesWithEnoughDistance;
	for (auto pair : shortestDistanceToCyclicComponent)
	{
		if (keptNodesInCycles.count(pair.first.id()) == 1) continue;
		if (shortestDistanceToCyclicComponent.count(reverse(pair.first)) == 1) continue;
		if (pair.second < minAnchorLength) continue;
		checkStack3.emplace_back(reverse(pair.first));
		nodesWithEnoughDistance.emplace(reverse(pair.first));
	}
	while (checkStack3.size() >= 1)
	{
		auto top = checkStack3.back();
		checkStack3.pop_back();
		if (keptNodesInAnchors.count(top.id()) == 1) continue;
		keptNodesInAnchors.insert(top.id());
		assert(graph.edges.count(top) == 1);
		for (auto edge : graph.edges.at(top))
		{
			if (keptNodesInCycles.count(std::get<0>(edge).id()) == 1) continue;
			if (shortestDistanceToCyclicComponent.count(reverse(std::get<0>(edge))) == 0) continue;
			checkStack3.emplace_back(std::get<0>(edge));
		}
	}
	for (Node node : nodesWithEnoughDistance)
	{
		assert(keptNodesInAnchors.count(node.id()) == 1);
		assert(graph.edges.count(node) == 1);
		bool allNeighborsHaveEnoughDistance = true;
		for (auto edge : graph.edges.at(node))
		{
			if (nodesWithEnoughDistance.count(std::get<0>(edge)) == 0)
			{
				allNeighborsHaveEnoughDistance = false;
				break;
			}
		}
		if (allNeighborsHaveEnoughDistance)
		{
			keptNodesInAnchors.erase(node.id());
		}
	}
	phmap::flat_hash_set<size_t> finalKeptNodes;
	finalKeptNodes.insert(keptNodesInCycles.begin(), keptNodesInCycles.end());
	finalKeptNodes.insert(keptNodesInAnchors.begin(), keptNodesInAnchors.end());
	std::ofstream file { outputGraphName };
	for (size_t node : finalKeptNodes)
	{
		file << "S\t" << graph.nodeNames[node] << "\t" << graph.nodeSeqs[node] << "\tll:f:" << graph.nodeCoverages[node] << "\tFC:i:" << (graph.nodeCoverages[node] * graph.nodeSeqs[node].size()) << std::endl;
	}
	for (const auto& pair : graph.edges)
	{
		if (finalKeptNodes.count(pair.first.id()) == 0) continue;
		for (const auto& target : pair.second)
		{
			if (finalKeptNodes.count(std::get<0>(target).id()) == 0) continue;
			file << "L\t" << graph.nodeNames[pair.first.id()] << "\t" << (pair.first.forward() ? "+" : "-") << "\t" << graph.nodeNames[std::get<0>(target).id()] << "\t" << (std::get<0>(target).forward() ? "+" : "-") << "\t" << std::get<1>(target) << "M\tec:i:" << std::get<2>(target) << std::endl;
		}
	}
}

void HandleCluster(const ClusterParams& params)
{
	std::filesystem::create_directories(params.basePath + "/tmp");
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
	std::cerr << "process graph" << std::endl;
	processGraphAndWrite(heavyPath, graph, params.basePath + "/processed-graph.gfa");
	// std::cerr << "getting variants" << std::endl;
	// std::vector<Variant> variants = getVariants(graph, heavyPath, readPaths, 3);
	// addNonrDNAAnchorsAsVariants(variants, graph, heavyPath, readPaths);
	// nameVariants(variants, graph, heavyPath);
	// std::cerr << variants.size() << " variants" << std::endl;
	// std::cerr << "writing variants" << std::endl;
	// writeVariants(heavyPath, graph, variants, params.basePath + "/variants.txt");
	// std::cerr << "writing variant graph" << std::endl;
	// writeVariantGraph(params.basePath + "/graph.gfa", graph, heavyPath, variants, params.basePath + "/variant-graph.gfa");
	// std::cerr << "writing allele graph" << std::endl;
	// writeAlleleGraph(graph, heavyPath, variants, params.basePath + "/allele-graph.gfa");
	// std::cerr << "writing variant vcf" << std::endl;
	// writeVariantVCF(params.basePath + "/variants.vcf", heavyPath, graph, variants);
	if (params.annotationFasta.size() > 0)
	{
		std::cerr << "lifting over annotations to consensus" << std::endl;
		liftoverAnnotationsToConsensus(params.basePath, params.basePath + "/consensus.fa", params.annotationFasta, params.annotationGff3, params.basePath + "/tmp", params.liftoffPath);
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
			assert(loop.originalReadLength != 0);
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
				readPaths.back().readLength = originalReadLength.at(pair.first);
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

phmap::flat_hash_map<uint64_t, size_t> getConsensusKmerPositions(const std::string& consensusSeq, const size_t refiningK)
{
	const uint64_t mask = (1ull << (2ull*refiningK))-1;
	phmap::flat_hash_map<uint64_t, size_t> result;
	assert(consensusSeq.size() > refiningK);
	uint64_t kmer = 0;
	for (size_t i = 0; i < refiningK; i++)
	{
		kmer <<= 2;
		switch(consensusSeq[i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
	}
	result[kmer] = 0;
	for (size_t i = refiningK; i < consensusSeq.size(); i++)
	{
		kmer <<= 2;
		kmer &= mask;
		switch(consensusSeq[i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
		if (result.count(kmer) == 1)
		{
			result[kmer] = std::numeric_limits<size_t>::max();
		}
		else
		{
			result[kmer] = i-refiningK+1;
		}
	}
	phmap::flat_hash_set<uint64_t> removeThese;
	for (auto pair : result)
	{
		if (pair.second != std::numeric_limits<size_t>::max()) continue;
		removeThese.emplace(pair.first);
	}
	for (uint64_t kmer : removeThese)
	{
		assert(result.count(kmer) == 1);
		result.erase(kmer);
	}
	return result;
}

std::pair<size_t, size_t> refineStartEndPoses(const std::string& readSequence, const std::string& consensusSequence, const phmap::flat_hash_map<uint64_t, size_t>& consensusKmerPositions, const size_t refiningK, const size_t initialStartPos, const size_t initialEndPos)
{
	const uint64_t mask = (1ull << (2ull*refiningK)) - 1ull;
	size_t checkAreaStart = initialStartPos;
	if (checkAreaStart > 5000)
	{
		checkAreaStart -= 5000;
	}
	else
	{
		checkAreaStart = 0;
	}
	size_t checkAreaEnd = std::min(initialEndPos + 5000, readSequence.size());
	uint64_t kmer = 0;
	for (size_t i = 0; i < refiningK; i++)
	{
		kmer <<= 2;
		switch(readSequence[checkAreaStart+i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
	}
	std::vector<std::pair<size_t, size_t>> matches;
	if (consensusKmerPositions.count(kmer) == 1)
	{
		matches.emplace_back(checkAreaStart, consensusKmerPositions.at(kmer));
	}
	for (size_t i = refiningK; checkAreaStart+i < checkAreaEnd; i++)
	{
		kmer <<= 2;
		kmer &= mask;
		switch(readSequence[checkAreaStart+i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
		if (consensusKmerPositions.count(kmer) == 1)
		{
			matches.emplace_back(checkAreaStart+i-refiningK+1, consensusKmerPositions.at(kmer));
		}
	}
	removeNonDiagonalKmers(matches);
	removeRepeatKmersEitherRepeat(matches);
	matches = getIncreasingChain(matches);
	if (matches.size() < 100)
	{
		return std::make_pair(initialStartPos, initialEndPos);
	}
	if (matches[0].second > 200)
	{
		return std::make_pair(initialStartPos, initialEndPos);
	}
	if (matches.back().second+200 < consensusSequence.size())
	{
		return std::make_pair(initialStartPos, initialEndPos);
	}
	size_t approxStart = matches[0].first;
	if (matches[0].second < approxStart)
	{
		approxStart -= matches[0].second;
	}
	else
	{
		approxStart = 0;
	}
	size_t approxEnd = std::min(matches.back().first + (consensusSequence.size() - matches.back().second), readSequence.size());
	approxStart = getExactBreakPos(readSequence, consensusSequence, approxStart);
	approxEnd = getExactBreakPos(readSequence, consensusSequence, approxEnd);
	return std::make_pair(approxStart, approxEnd);
}

std::vector<std::vector<std::string>> getRawLoopSequences(const std::vector<std::vector<OntLoop>>& clusters, const std::string& rawOntPath, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t numThreads)
{
	const size_t refiningK = 15;
	phmap::flat_hash_map<std::string, std::vector<std::tuple<size_t, size_t, size_t, size_t>>> loopPositionsInReads;
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> clusterConsensusKmerPositions;
	std::vector<std::string> clusterConsensuses;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		auto path = getConsensusPath(clusters[i], graph);
		clusterConsensuses.emplace_back(getConsensusSequence(path, graph, pathStartClip, pathEndClip));
		clusterConsensusKmerPositions.emplace_back(getConsensusKmerPositions(clusterConsensuses.back(), refiningK));
	}
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
	std::vector<std::thread> threads;
	std::atomic<bool> readDone;
	readDone = false;
	std::mutex stackMutex;
	std::vector<FastQ> readStack;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&readDone, &stackMutex, &readStack, &result, &clusterConsensusKmerPositions, &clusterConsensuses, &loopPositionsInReads]()
		{
			while (true)
			{
				FastQ fastq;
				{
					std::lock_guard<std::mutex> lock { stackMutex };
					if (readStack.size() == 0)
					{
						if (readDone)
						{
							if (readStack.size() == 0)
							{
								break;
							}
						}
						continue;
					}
					std::swap(readStack.back(), fastq);
					readStack.pop_back();
				}
				std::string readname = nameWithoutTags(fastq.seq_id);
				if (loopPositionsInReads.count(readname) == 0) continue;
				std::string revcompRead = revcomp(fastq.sequence);
				for (auto t : loopPositionsInReads.at(readname))
				{
					size_t startPos = std::get<2>(t);
					size_t endPos = std::get<3>(t);
					bool reverse = false;
					if (startPos > endPos)
					{
						reverse = true;
					}
					if (reverse)
					{
						startPos = fastq.sequence.size() - startPos;
						endPos = fastq.sequence.size() - endPos;
						assert(endPos > startPos);
						std::tie(startPos, endPos) = refineStartEndPoses(revcompRead, clusterConsensuses[std::get<0>(t)], clusterConsensusKmerPositions[std::get<0>(t)], refiningK, startPos, endPos);
						startPos = fastq.sequence.size() - startPos;
						endPos = fastq.sequence.size() - endPos;
						std::swap(startPos, endPos);
						assert(endPos > startPos);
					}
					else
					{
						std::tie(startPos, endPos) = refineStartEndPoses(fastq.sequence, clusterConsensuses[std::get<0>(t)], clusterConsensusKmerPositions[std::get<0>(t)], refiningK, startPos, endPos);
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
			}
		});
	}
	FastQ::streamFastqFromFile(rawOntPath, false, [&readStack, &stackMutex](FastQ& fastq)
	{
		std::lock_guard<std::mutex> lock { stackMutex };
		readStack.emplace_back();
		std::swap(readStack.back(), fastq);
	});
	readDone = true;
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].size(); j++)
		{
			assert(result[i][j] != "");
		}
	}
	return result;
}

void writeSelfCorrectedOntLoopSequences(const std::string outputFile, const std::vector<std::vector<OntLoop>>& loopSequences)
{
	std::ofstream file { outputFile };
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		for (size_t j = 0; j < loopSequences[i].size(); j++)
		{
			file << ">" << loopSequences[i][j].rawLoopName << std::endl;
			file << loopSequences[i][j].selfCorrectedSequence << std::endl;
		}
	}
}

void writeRawOntLoopSequences(const std::string outputFile, const std::vector<std::vector<OntLoop>>& loopSequences)
{
	std::ofstream file { outputFile };
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		for (size_t j = 0; j < loopSequences[i].size(); j++)
		{
			file << ">" << loopSequences[i][j].rawLoopName << std::endl;
			file << loopSequences[i][j].rawSequence << std::endl;
		}
	}
}

void alignSelfCorrectedOntLoopsToMorphConsensuses(const std::vector<std::vector<OntLoop>>& loopSequences, std::string outputFile, std::string tmppath, const size_t numThreads, const std::vector<MorphConsensus>& morphConsensuses, const std::string& winnowmapPath, const std::string& samtoolsPath)
{
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		{
			std::ofstream tmpReadsFile { tmppath + "/tmpreads.fa" };
			std::ofstream tmpConsensusFile { tmppath + "/tmpconsensus.fa" };
			for (size_t j = 0; j < loopSequences[i].size(); j++)
			{
				tmpReadsFile << ">" << loopSequences[i][j].rawLoopName << std::endl;
				tmpReadsFile << loopSequences[i][j].selfCorrectedSequence << std::endl;
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

void alignRawOntLoopsToMorphConsensuses(const std::vector<std::vector<OntLoop>>& loopSequences, std::string outputFile, std::string tmppath, const size_t numThreads, const std::vector<MorphConsensus>& morphConsensuses, const std::string& winnowmapPath, const std::string& samtoolsPath)
{
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		{
			std::ofstream tmpReadsFile { tmppath + "/tmpreads.fa" };
			std::ofstream tmpConsensusFile { tmppath + "/tmpconsensus.fa" };
			for (size_t j = 0; j < loopSequences[i].size(); j++)
			{
				tmpReadsFile << ">" << loopSequences[i][j].rawLoopName << std::endl;
				tmpReadsFile << loopSequences[i][j].rawSequence << std::endl;
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

bool isSNP(const std::tuple<size_t, size_t, std::string>& variant)
{
	if (std::get<1>(variant)-std::get<0>(variant) == 1 && std::get<2>(variant).size() == 1) return true;
	return false;
}

bool isBigIndel(const std::tuple<size_t, size_t, std::string>& variant)
{
	if (std::get<1>(variant)-std::get<0>(variant) >= 10) return true;
	if (std::get<2>(variant).size() >= 10) return true;
	return false;
}

template <typename F>
void iterateKmers(const std::string seq, const size_t k, F callback)
{
	if (seq.size() < k) return;
	const __uint128_t mask = ((__uint128_t)1 << (__uint128_t)(2ull * k)) - (__uint128_t)1;
	__uint128_t kmer = 0;
	for (size_t i = 0; i < k; i++)
	{
		kmer <<= 2;
		switch(seq[i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
	}
	callback(kmer, 0);
	for (size_t i = k; i < seq.size(); i++)
	{
		kmer <<= 2;
		kmer &= mask;
		switch(seq[i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
		callback(kmer, i-k+1);
	}
}

std::vector<std::vector<size_t>> getAllPresentKmerChain(const std::string consensusSeq, const std::vector<std::string>& loopSequences, const size_t k)
{
	phmap::flat_hash_map<uint64_t, size_t> kmerPosition;
	phmap::flat_hash_set<uint64_t> duplicateKmers;
	phmap::flat_hash_set<uint64_t> kmersEverywhere;
	iterateKmers(consensusSeq, k, [&kmerPosition, &duplicateKmers](const uint64_t kmer, const size_t position)
	{
		if (kmerPosition.count(kmer) == 1)
		{
			duplicateKmers.insert(kmer);
			return;
		}
		kmerPosition[kmer] = position;
	});
	for (auto pair : kmerPosition)
	{
		if (duplicateKmers.count(pair.first) == 1) continue;
		kmersEverywhere.insert(pair.first);
	}
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		phmap::flat_hash_set<uint64_t> kmersHere;
		iterateKmers(loopSequences[i], k, [&duplicateKmers, &kmersHere](const uint64_t kmer, const size_t position)
		{
			if (kmersHere.count(kmer) == 1)
			{
				duplicateKmers.insert(kmer);
				return;
			}
			kmersHere.insert(kmer);
		});
		for (auto kmer : kmersEverywhere)
		{
			if (kmersHere.count(kmer) == 0) duplicateKmers.insert(kmer);
		}
	}
	std::vector<std::pair<size_t, uint64_t>> kmerPositionsInConsensus;
	for (auto pair : kmerPosition)
	{
		if (duplicateKmers.count(pair.first) == 1) continue;
		kmerPositionsInConsensus.emplace_back(pair.second, pair.first);
	}
	std::sort(kmerPositionsInConsensus.begin(), kmerPositionsInConsensus.end());
	phmap::flat_hash_map<uint64_t, size_t> kmerToIndex;
	for (size_t i = 0; i < kmerPositionsInConsensus.size(); i++)
	{
		kmerToIndex[kmerPositionsInConsensus[i].second] = i;
	}
	std::vector<std::vector<size_t>> kmerPositionWithinLoop;
	kmerPositionWithinLoop.resize(kmerPositionsInConsensus.size());
	for (size_t i = 0; i < kmerPositionWithinLoop.size(); i++)
	{
		kmerPositionWithinLoop[i].resize(loopSequences.size(), std::numeric_limits<size_t>::max());
	}
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		iterateKmers(loopSequences[i], k, [&kmerPositionWithinLoop, &kmerToIndex, i](const uint64_t kmer, const size_t position)
		{
			if (kmerToIndex.count(kmer) == 0) return;
			assert(kmerPositionWithinLoop[kmerToIndex.at(kmer)][i] == std::numeric_limits<size_t>::max());
			kmerPositionWithinLoop[kmerToIndex.at(kmer)][i] = position;
		});
	}
	for (size_t i = 0; i < kmerPositionWithinLoop.size(); i++)
	{
		for (size_t j = 0; j < kmerPositionWithinLoop[i].size(); j++)
		{
			assert(kmerPositionWithinLoop[i][j] != std::numeric_limits<size_t>::max());
		}
	}
	std::cerr << "cluster with " << loopSequences.size() << " reads has " << kmerPositionWithinLoop.size() << " all-present kmers before chain filter" << std::endl;
	std::cerr << "positions before chain checking:";
	for (auto pair : kmerPositionsInConsensus)
	{
		std::cerr << " " << pair.first;
	}
	std::cerr << std::endl;
	while (true)
	{
		std::vector<size_t> indexConflictCount;
		indexConflictCount.resize(kmerPositionWithinLoop.size(), false);
		for (size_t i = 1; i < kmerPositionWithinLoop.size(); i++)
		{
			assert(kmerPositionWithinLoop[i].size() == loopSequences.size());
			assert(kmerPositionWithinLoop[i-1].size() == loopSequences.size());
			for (size_t j = 0; j < loopSequences.size(); j++)
			{
				assert(kmerPositionWithinLoop[i-1][j] != kmerPositionWithinLoop[i][j]);
				if (kmerPositionWithinLoop[i-1][j] < kmerPositionWithinLoop[i][j]) continue;
				indexConflictCount[i] += 1;
				indexConflictCount[i-1] += 1;
			}
		}
		std::vector<bool> removeIndex;
		removeIndex.resize(kmerPositionWithinLoop.size(), false);
		for (size_t i = 0; i < kmerPositionWithinLoop.size(); i++)
		{
			if (indexConflictCount[i] == 0) continue;
			if (i > 0 && indexConflictCount[i-1] > indexConflictCount[i]) continue;
			if (i+1 < kmerPositionWithinLoop.size() && indexConflictCount[i+1] > indexConflictCount[i]) continue;
			removeIndex[i] = true;
		}
		bool removedAny = false;
		for (size_t i = kmerPositionWithinLoop.size()-1; i < kmerPositionWithinLoop.size(); i--)
		{
			if (!removeIndex[i]) continue;
			kmerPositionWithinLoop.erase(kmerPositionWithinLoop.begin() + i);
			removedAny = true;
		}
		if (!removedAny) break;
	}
	return kmerPositionWithinLoop;
}

void clusterBubbleAlleles(std::vector<std::vector<size_t>>& result, const std::vector<std::string>& loopSequences, const std::vector<size_t>& previousMatchIndices, const std::vector<size_t>& currentMatchIndices, const size_t k)
{
	const size_t minCoverage = 5;
	const double maxDivergence = 0.1;
	size_t maxClusterDistance = 10;
	if (currentMatchIndices[0] < previousMatchIndices[0] + k)
	{
		bool allIdentical = true;
		for (size_t i = 1; i < loopSequences.size(); i++)
		{
			if (currentMatchIndices[i] - previousMatchIndices[i] != currentMatchIndices[0] - previousMatchIndices[0])
			{
				allIdentical = false;
			}
		}
		if (allIdentical) return;
	}
	std::vector<std::string> substrings;
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		assert(loopSequences[i].size() >= currentMatchIndices[i]+k);
		assert(currentMatchIndices[i] > previousMatchIndices[i]);
		substrings.emplace_back(loopSequences[i].substr(previousMatchIndices[i], currentMatchIndices[i]-previousMatchIndices[i]+k));
		assert(substrings.back().size() >= k+1);
		maxClusterDistance = std::max<size_t>(maxClusterDistance, substrings[i].size() * maxDivergence);
	}
	std::vector<std::vector<size_t>> distanceMatrix;
	distanceMatrix.emplace_back();
	size_t maxDistance = 0;
	for (size_t i = 1; i < loopSequences.size(); i++)
	{
		distanceMatrix.emplace_back();
		for (size_t j = 0; j < i; j++)
		{
			distanceMatrix.back().emplace_back(getEditDistanceWfa(substrings[i], substrings[j]));
			maxDistance = std::max(maxDistance, distanceMatrix.back().back());
		}
	}
	if (maxDistance < maxClusterDistance) return;
	std::vector<size_t> parent;
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		parent.emplace_back(i);
	}
	for (size_t i = 1; i < loopSequences.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (distanceMatrix[i][j] > maxClusterDistance) continue;
			merge(parent, i, j);
		}
	}
	phmap::flat_hash_set<size_t> distinctClusters;
	std::vector<std::vector<size_t>> clusters;
	clusters.resize(loopSequences.size());
	for (size_t i = 0; i < parent.size(); i++)
	{
		clusters[find(parent, i)].emplace_back(i);
		distinctClusters.insert(find(parent, i));
	}
	if (distinctClusters.size() < 2) return;
	bool allClustersAreCliques = true;
	bool hasSmallCluster = false;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		if (clusters[i].size() > 0 && clusters[i].size() < minCoverage)
		{
			hasSmallCluster = true;
			break;
		}
		for (size_t j = 1; j < clusters[i].size(); j++)
		{
			for (size_t k = 0; k < j; k++)
			{
				size_t smaller = std::min(clusters[i][j], clusters[i][k]);
				size_t bigger = std::max(clusters[i][j], clusters[i][k]);
				if (distanceMatrix[bigger][smaller] > maxClusterDistance)
				{
					allClustersAreCliques = false;
					break;
				}
			}
			if (!allClustersAreCliques) break;
		}
		if (!allClustersAreCliques) break;
	}
	if (hasSmallCluster) return;
	if (!allClustersAreCliques) return;
	// distance between clusters is always at least maxClusterDistance+1 because of union-find
	for (size_t i = 0; i < clusters.size(); i++)
	{
		for (size_t read : clusters[i])
		{
			result[read].emplace_back(i);
		}
	}
}

std::vector<std::vector<size_t>> splitByBubbleAlleles(const std::vector<OntLoop>& loops, const std::vector<std::string>& loopSequences, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t numThreads)
{
	const size_t k = 15;
	auto path = getConsensusPath(loops, graph);
	std::string consensusSeq = getConsensusSequence(path, graph, pathStartClip, pathEndClip);
	consensusSeq = getPolishedSequence(consensusSeq, loopSequences, numThreads);
	std::vector<std::vector<size_t>> kmerChain = getAllPresentKmerChain(consensusSeq, loopSequences, k);
	std::cerr << "cluster with " << loopSequences.size() << " reads has " << kmerChain.size() << " all-present kmers in chain" << std::endl;
	std::vector<std::vector<size_t>> allelesPerRead;
	allelesPerRead.resize(loopSequences.size());
	for (size_t i = 1; i < kmerChain.size(); i++)
	{
		clusterBubbleAlleles(allelesPerRead, loopSequences, kmerChain[i-1], kmerChain[i], k);
	}
	std::cerr << "cluster with " << loopSequences.size() << " reads has " << allelesPerRead[0].size() << " alleles" << std::endl;
	std::map<std::vector<size_t>, size_t> allelesToCluster;
	std::vector<std::vector<size_t>> result;
	for (size_t i = 0; i < allelesPerRead.size(); i++)
	{
		if (allelesToCluster.count(allelesPerRead[i]) == 0)
		{
			allelesToCluster[allelesPerRead[i]] = result.size();
			result.emplace_back();
		}
		result[allelesToCluster.at(allelesPerRead[i])].emplace_back(i);
	}
	return result;
}

std::vector<std::vector<std::tuple<size_t, size_t, std::string>>> getEditsForPhasing(const std::vector<OntLoop>& loops, const std::vector<std::string>& loopSequences, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t numThreads)
{
	auto path = getConsensusPath(loops, graph);
	std::string consensusSeq = getConsensusSequence(path, graph, pathStartClip, pathEndClip);
	consensusSeq = getPolishedSequence(consensusSeq, loopSequences, numThreads);
	size_t firstMatchPos = 0;
	size_t lastMatchPos = consensusSeq.size();;
	std::vector<std::vector<std::tuple<size_t, size_t, std::string>>> editsPerRead;
	editsPerRead.resize(loops.size());
	std::mutex resultMutex;
	iterateEdits(consensusSeq, loopSequences, numThreads, [&firstMatchPos, &lastMatchPos, &resultMutex, &editsPerRead](const size_t threadId, const size_t readId, const size_t firstMatchRef, const size_t lastMatchRef, const size_t firstMatchRead, const size_t lastMatchRead, const std::tuple<size_t, size_t, std::string, size_t>& edit)
	{
		auto filterEdit = std::make_tuple(std::get<0>(edit), std::get<1>(edit), std::get<2>(edit));
		editsPerRead[readId].emplace_back(filterEdit);
		std::lock_guard<std::mutex> lock { resultMutex };
		firstMatchPos = std::max(firstMatchPos, firstMatchRef);
		lastMatchPos = std::min(lastMatchPos, lastMatchRef);
	});
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (size_t j = editsPerRead[i].size()-1; j < editsPerRead[i].size(); j++)
		{
			if (std::get<0>(editsPerRead[i][j]) < firstMatchPos || std::get<1>(editsPerRead[i][j]) > lastMatchPos)
			{
				std::swap(editsPerRead[i][j], editsPerRead[i].back());
				editsPerRead[i].pop_back();
			}
		}
	}
	return editsPerRead;
}

std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>> getSNPMSA(const std::vector<std::vector<std::tuple<size_t, size_t, std::string>>>& editsPerRead)
{
	const size_t minCoverage = 3;
	const size_t minCoverageWithDels = 5;
	std::vector<std::vector<std::vector<size_t>>> readsWhichHaveAlt;
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (size_t j = 0; j < editsPerRead[i].size(); j++)
		{
			while (readsWhichHaveAlt.size() < std::get<1>(editsPerRead[i][j])) readsWhichHaveAlt.emplace_back();
			const std::tuple<size_t, size_t, std::string>& t = editsPerRead[i][j];
			if (std::get<1>(t) == std::get<0>(t)) continue;
			if (std::get<2>(t).size() >= 1 && std::get<1>(t)-std::get<0>(t) != std::get<2>(t).size()) continue;
			for (size_t k = std::get<0>(t); k < std::get<1>(t); k++)
			{
				size_t allele = 5;
				if (k-std::get<0>(t) < std::get<2>(t).size())
				{
					switch(std::get<2>(t)[k-std::get<0>(t)])
					{
						case 'A':
							allele = 0;
							break;
						case 'C':
							allele = 1;
							break;
						case 'G':
							allele = 2;
							break;
						case 'T':
							allele = 3;
							break;
						default:
							assert(false);
							break;
					}
				}
				while (readsWhichHaveAlt[k].size() <= allele) readsWhichHaveAlt[k].emplace_back(0);
				readsWhichHaveAlt[k][allele].emplace_back(i);
			}
		}
	}
	std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>> result;
	for (size_t i = 0; i < readsWhichHaveAlt.size(); i++)
	{
		size_t totalCoverage = 0;
		bool hasSmallAllele = false;
		size_t totalNonDeletionAlleles = 0;
		bool hasDeletion = readsWhichHaveAlt[i].size() > 5 && readsWhichHaveAlt[i][5].size() >= 1;
		for (size_t allele = 0; allele < readsWhichHaveAlt[i].size(); allele++)
		{
			totalCoverage += readsWhichHaveAlt[i][allele].size();
			if (readsWhichHaveAlt[i][allele].size() == 0) continue;
			if (allele != 5) totalNonDeletionAlleles += 1;
			if (hasDeletion && readsWhichHaveAlt[i][allele].size() < minCoverageWithDels && allele != 5)
			{
				hasSmallAllele = true;
			}
			if (readsWhichHaveAlt[i][allele].size() < minCoverage && allele != 5)
			{
				hasSmallAllele = true;
			}
		}
		if (hasSmallAllele) continue;
		if (totalCoverage+minCoverage > editsPerRead.size()) continue;
		if (hasDeletion && totalCoverage+minCoverageWithDels > editsPerRead.size()) continue;
		if (totalNonDeletionAlleles < 1) continue;
		result.emplace_back();
		result.back().second = i;
		result.back().first.resize(7);
		std::vector<bool> readHasAlt;
		readHasAlt.resize(editsPerRead.size(), false);
		for (size_t allele = 0; allele < readsWhichHaveAlt[i].size(); allele++)
		{
			for (size_t read : readsWhichHaveAlt[i][allele])
			{
				result.back().first[allele].emplace_back(read);
				readHasAlt[read] = true;
			}
		}
		for (size_t j = 0; j < readHasAlt.size(); j++)
		{
			if (readHasAlt[j]) continue;
			result.back().first[6].emplace_back(j);
		}
		assert(result.back().first[6].size() >= minCoverage);
		assert(!hasDeletion || result.back().first[6].size() >= minCoverageWithDels);
		for (size_t j = result.back().first.size()-1; j < result.back().first.size(); j--)
		{
			if (result.back().first[j].size() != 0) continue;
			std::swap(result.back().first[j], result.back().first.back());
			result.back().first.pop_back();
		}
		assert(result.back().first.size() >= 2);
		assert((!hasDeletion && result.back().first.size() == totalNonDeletionAlleles+1) || (hasDeletion && result.back().first.size() == totalNonDeletionAlleles+2));
	}
	return result;
}

std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>> getPhasableVariantInfoBiallelicAltRefSNPsBigIndels(const std::vector<std::vector<std::tuple<size_t, size_t, std::string>>>& editsPerRead)
{
	const size_t minCoverage = 3;
	std::map<std::tuple<size_t, size_t, std::string>, size_t> editCounts;
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (size_t j = 0; j < editsPerRead[i].size(); j++)
		{
			editCounts[editsPerRead[i][j]] += 1;
		}
	}
	std::vector<std::tuple<size_t, size_t, std::string>> distinctEdits;
	for (auto pair : editCounts)
	{
		distinctEdits.emplace_back(pair.first);
	}
	std::sort(distinctEdits.begin(), distinctEdits.end());
	std::set<std::tuple<size_t, size_t, std::string>> discardedEdits;
	size_t currentClusterStart = 0;
	size_t currentClusterLastPos = std::get<1>(distinctEdits[0]);
	std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>> result;
	std::map<std::tuple<size_t, size_t, std::string>, std::pair<size_t, size_t>> validEditIndices;
	for (size_t i = 1; i <= distinctEdits.size(); i++)
	{
		if (i == distinctEdits.size() || std::get<0>(distinctEdits[i]) > currentClusterLastPos)
		{
			if (i-currentClusterStart <= 1)
			{
				size_t countSNPs = 0;
				size_t countBigIndels = 0;
				for (size_t j = currentClusterStart; j < i; j++)
				{
					if (isSNP(distinctEdits[j])) countSNPs += 1;
					if (isBigIndel(distinctEdits[j])) countBigIndels += 1;
				}
				if (countSNPs >= 1 || countBigIndels >= 1)
				{
					result.emplace_back();
					result.back().first.emplace_back(); // ref allele
					size_t minPos = std::get<0>(distinctEdits[currentClusterStart]);
					size_t maxPos = std::get<1>(distinctEdits[currentClusterStart]);
					for (size_t j = currentClusterStart; j < i; j++)
					{
						result.back().first.emplace_back();
						validEditIndices[distinctEdits[j]] = std::make_pair(result.size()-1, result.back().first.size()-1);
						minPos = std::min(minPos, std::get<0>(distinctEdits[j]));
						maxPos = std::max(maxPos, std::get<1>(distinctEdits[j]));
					}
					assert(result.back().first.size() >= 2);
					result.back().second = (minPos+maxPos)/2;
				}
			}
			currentClusterStart = i;
		}
		if (i != distinctEdits.size()) currentClusterLastPos = std::max(currentClusterLastPos, std::get<1>(distinctEdits[i]));
	}
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (const auto& t : editsPerRead[i])
		{
			if (validEditIndices.count(t) == 0) continue;
			auto index = validEditIndices.at(t);
			result[index.first].first[index.second].emplace_back(i);
		}
	}
	for (size_t i = result.size()-1; i < result.size(); i--)
	{
		bool hasMultipleAllelesFromSameRead = false;
		std::vector<bool> hasAlt;
		hasAlt.resize(editsPerRead.size(), false);
		for (size_t allele = 1; allele < result[i].first.size(); allele++)
		{
			for (size_t j : result[i].first[allele])
			{
				if (hasAlt[j])
				{
					hasMultipleAllelesFromSameRead = true;
					break;
				}
				hasAlt[j] = true;
			}
			if (hasMultipleAllelesFromSameRead) break;
		}
		if (hasMultipleAllelesFromSameRead)
		{
			std::swap(result[i], result.back());
			result.pop_back();
			continue;
		}
		for (size_t j = 0; j < hasAlt.size(); j++)
		{
			if (hasAlt[j]) continue;
			result[i].first[0].emplace_back(j);
		}
	}
	for (size_t i = result.size()-1; i < result.size(); i--)
	{
		bool hasLowCoverage = false;
		for (size_t j = 0; j < result[i].first.size(); j++)
		{
			if (result[i].first[j].size() >= minCoverage) continue;
			hasLowCoverage = true;
			break;
		}
		if (hasLowCoverage)
		{
			std::swap(result[i], result.back());
			result.pop_back();
		}
	}
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) { return left.second < right.second; });
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].first.size(); j++)
		{
			std::sort(result[i].first[j].begin(), result[i].first[j].end());
		}
	}
	return result;
}

std::pair<std::vector<std::vector<OntLoop>>, std::vector<std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>>>> splitLoopsAndPhaseInfoToClusters(const std::vector<OntLoop>& previous, const std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>>& phasableVariantInfo, const std::vector<std::vector<size_t>>& splitted)
{
	const size_t minCoverage = 3;
	std::vector<size_t> clusterOfLoop;
	std::vector<size_t> readIndexWithinCluster;
	readIndexWithinCluster.resize(previous.size(), std::numeric_limits<size_t>::max());
	clusterOfLoop.resize(previous.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < splitted.size(); i++)
	{
		for (size_t j = 0; j < splitted[i].size(); j++)
		{
			assert(splitted[i][j] < previous.size());
			assert(clusterOfLoop[splitted[i][j]] == std::numeric_limits<size_t>::max());
			assert(readIndexWithinCluster[splitted[i][j]] == std::numeric_limits<size_t>::max());
			clusterOfLoop[splitted[i][j]] = i;
			readIndexWithinCluster[splitted[i][j]] = j;
		}
	}
	for (size_t i = 0; i < readIndexWithinCluster.size(); i++)
	{
		assert(clusterOfLoop[i] != std::numeric_limits<size_t>::max());
		assert(readIndexWithinCluster[i] != std::numeric_limits<size_t>::max());
	}
	std::vector<std::vector<OntLoop>> splittedLoops;
	std::vector<std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>>> splittedPhasableVariantInfo;
	splittedLoops.resize(splitted.size());
	splittedPhasableVariantInfo.resize(splitted.size());
	for (size_t i = 0; i < splitted.size(); i++)
	{
		splittedLoops[i].resize(splitted[i].size());
	}
	for (size_t j = 0; j < previous.size(); j++)
	{
		assert(j < clusterOfLoop.size());
		assert(clusterOfLoop[j] < splittedLoops.size());
		assert(j < readIndexWithinCluster.size());
		assert(readIndexWithinCluster[j] < splittedLoops[clusterOfLoop[j]].size());
		splittedLoops[clusterOfLoop[j]][readIndexWithinCluster[j]] = previous[j];
	}
	for (size_t site = 0; site < phasableVariantInfo.size(); site++)
	{
		for (size_t j = 0; j < splitted.size(); j++)
		{
			splittedPhasableVariantInfo[j].emplace_back();
			splittedPhasableVariantInfo[j][site].first.resize(phasableVariantInfo[site].first.size());
			splittedPhasableVariantInfo[j][site].second = phasableVariantInfo[site].second;
		}
		for (size_t allele = 0; allele < phasableVariantInfo[site].first.size(); allele++)
		{
			for (size_t read : phasableVariantInfo[site].first[allele])
			{
				splittedPhasableVariantInfo[clusterOfLoop[read]][site].first[allele].emplace_back(readIndexWithinCluster[read]);
			}
		}
	}
	for (size_t i = 0; i < splittedPhasableVariantInfo.size(); i++)
	{
		bool removedAny = false;
		for (size_t j = splittedPhasableVariantInfo[i].size()-1; j < splittedPhasableVariantInfo[i].size(); j--)
		{
			bool hasLowCoverage = false;
			for (size_t k = 0; k < splittedPhasableVariantInfo[i][j].first.size(); k++)
			{
				if (splittedPhasableVariantInfo[i][j].first[k].size() < minCoverage) hasLowCoverage = true;
			}
			if (hasLowCoverage)
			{
				std::swap(splittedPhasableVariantInfo[i][j], splittedPhasableVariantInfo[i].back());
				splittedPhasableVariantInfo[i].pop_back();
				removedAny = true;
			}
		}
		if (removedAny)
		{
			std::sort(splittedPhasableVariantInfo[i].begin(), splittedPhasableVariantInfo[i].end(), [](const auto& left, const auto& right) { return left.second < right.second; });
		}
	}
	for (size_t i = 0; i < splittedPhasableVariantInfo.size(); i++)
	{
		for (size_t j = 0; j < splittedPhasableVariantInfo[i].size(); j++)
		{
			for (size_t k = 0; k < splittedPhasableVariantInfo[i][j].first.size(); k++)
			{
				for (size_t read : splittedPhasableVariantInfo[i][j].first[k])
				{
					assert(read < splitted[i].size());
				}
			}
		}
	}
	return std::make_pair(splittedLoops, splittedPhasableVariantInfo);
}

bool vectorsMatch(const std::vector<size_t>& left, const std::vector<size_t>& right)
{
	if (left.size() != right.size()) return false;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i] != right[i]) return false;
	}
	return true;
}

std::vector<std::vector<size_t>> trySplitTwoSites(const std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>>& phasableVariantInfo, const size_t numReads)
{
	const size_t minDistance = 100;
	const size_t minCoverage = 3;
	std::vector<bool> siteUsed;
	siteUsed.resize(phasableVariantInfo.size(), false);
	std::vector<std::vector<size_t>> phaseClustersPerRead;
	phaseClustersPerRead.resize(numReads);
	for (size_t i = 0; i < phasableVariantInfo.size(); i++)
	{
		for (size_t j = i+1; j < phasableVariantInfo.size(); j++)
		{
			if (siteUsed[i]) break;
			if (siteUsed[j]) continue;
			assert(phasableVariantInfo[i].second < phasableVariantInfo[j].second);
			if (phasableVariantInfo[j].second < phasableVariantInfo[i].second + minDistance) continue;
			std::vector<std::pair<size_t, size_t>> allelesPerRead;
			allelesPerRead.resize(numReads, std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
			for (size_t k = 0; k < phasableVariantInfo[i].first.size(); k++)
			{
				for (size_t read : phasableVariantInfo[i].first[k])
				{
					assert(allelesPerRead[read].first == std::numeric_limits<size_t>::max());
					allelesPerRead[read].first = k;
				}
			}
			for (size_t k = 0; k < phasableVariantInfo[j].first.size(); k++)
			{
				for (size_t read : phasableVariantInfo[j].first[k])
				{
					assert(allelesPerRead[read].second == std::numeric_limits<size_t>::max());
					allelesPerRead[read].second = k;
				}
			}
			phmap::flat_hash_set<std::pair<size_t, size_t>> alleles;
			phmap::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> parent;
			for (auto t : allelesPerRead)
			{
				alleles.insert(t);
				parent[t] = t;
			}
			for (auto t : alleles)
			{
				for (auto t2 : alleles)
				{
					if (t2 == t) continue;
					size_t edits = 0;
					if (t.first != t2.first) edits += 1;
					if (t.second != t2.second) edits += 1;
					assert(edits > 0);
					if (edits == 1) merge(parent, t, t2);
				}
			}
			phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> clusters;
			for (auto t : alleles)
			{
				auto p = find(parent, t);
				if (clusters.count(p) == 1) continue;
				size_t cluster = clusters.size();
				clusters[p] = cluster;
			}
			if (clusters.size() < 2) continue;
			phmap::flat_hash_map<size_t, size_t> clusterCoverage;
			for (auto t : allelesPerRead)
			{
				clusterCoverage[clusters.at(find(parent, t))] += 1;
			}
			bool hasSmallCluster = false;
			for (auto pair : clusterCoverage)
			{
				if (pair.second < minCoverage) hasSmallCluster = true;
			}
			if (hasSmallCluster) continue;
			std::cerr << "split two sites " << phasableVariantInfo[i].second << " " << phasableVariantInfo[j].second << std::endl;
			for (size_t k = 0; k < allelesPerRead.size(); k++)
			{
				size_t cluster = clusters.at(find(parent, allelesPerRead[k]));
				phaseClustersPerRead[k].emplace_back(cluster);
			}
			siteUsed[i] = true;
			siteUsed[j] = true;
		}
	}
	if (phaseClustersPerRead[0].size() == 0) return std::vector<std::vector<size_t>> {};
	std::vector<std::vector<size_t>> result;
	std::map<std::vector<size_t>, size_t> clusters;
	for (size_t i = 0; i < phaseClustersPerRead.size(); i++)
	{
		if (clusters.count(phaseClustersPerRead[i]) == 0)
		{
			clusters[phaseClustersPerRead[i]] = result.size();
			result.emplace_back();
		}
		result[clusters.at(phaseClustersPerRead[i])].emplace_back(i);
	}
	return result;
}

std::vector<std::vector<size_t>> trySplitThreeSites(const std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>>& phasableVariantInfo)
{
	const size_t minDistance = 100;
	const size_t minCoverage = 3;
	std::vector<std::vector<size_t>> allelesPerRead;
	for (size_t i = 0; i < phasableVariantInfo.size(); i++)
	{
		for (size_t j = 0; j < phasableVariantInfo[i].first.size(); j++)
		{
			for (auto read : phasableVariantInfo[i].first[j])
			{
				while (allelesPerRead.size() <= read)
				{
					allelesPerRead.emplace_back();
					allelesPerRead.back().resize(phasableVariantInfo.size(), std::numeric_limits<size_t>::max());
				}
				assert(allelesPerRead[read][i] == std::numeric_limits<size_t>::max());
				allelesPerRead[read][i] = j;
			}
		}
	}
	for (size_t i = 0; i < allelesPerRead.size(); i++)
	{
		for (size_t j = 0; j < allelesPerRead[i].size(); j++)
		{
			assert(allelesPerRead[i][j] != std::numeric_limits<size_t>::max());
		}
	}
	for (size_t blocki = 0; blocki*10 < (phasableVariantInfo.size()+9)/10; blocki++)
	{
		for (size_t blockj = blocki; blockj*10 < (phasableVariantInfo.size()+9)/10; blockj++)
		{
			for (size_t blockk = blockj; blockk*10 < (phasableVariantInfo.size()+9)/10; blockk++)
			{
				std::set<std::vector<size_t>> distinctBlocks;
				for (size_t i = 0; i < allelesPerRead.size(); i++)
				{
					std::vector<size_t> thisBlock;
					for (size_t j = 0; j < 10; j++)
					{
						if (blocki*10+j < phasableVariantInfo.size())
						{
							thisBlock.emplace_back(allelesPerRead[i][blocki*10+j]);
						}
						else
						{
							thisBlock.emplace_back(std::numeric_limits<size_t>::max());
						}
					}
					for (size_t j = 0; j < 10; j++)
					{
						if (blockj*10+j < phasableVariantInfo.size())
						{
							thisBlock.emplace_back(allelesPerRead[i][blockj*10+j]);
						}
						else
						{
							thisBlock.emplace_back(std::numeric_limits<size_t>::max());
						}
					}
					for (size_t j = 0; j < 10; j++)
					{
						if (blockk*10+j < phasableVariantInfo.size())
						{
							thisBlock.emplace_back(allelesPerRead[i][blockk*10+j]);
						}
						else
						{
							thisBlock.emplace_back(std::numeric_limits<size_t>::max());
						}
					}
					assert(thisBlock.size() == 30);
					distinctBlocks.emplace(thisBlock);
				}
				for (size_t offi = 0; offi < 10 && blocki*10+offi < phasableVariantInfo.size(); offi++)
				{
					for (size_t offj = (blocki == blockj ? offi+1 : 0); offj < 10 && blockj*10+offj < phasableVariantInfo.size(); offj++)
					{
						for (size_t offk = (blockj == blockk ? offj+1 : 0); offk < 10 && blockk*10+offk < phasableVariantInfo.size(); offk++)
						{
							size_t i = blocki*10+offi;
							size_t j = blockj*10+offj;
							size_t k = blockk*10+offk;
							assert(i < j);
							assert(j < k);
							assert(phasableVariantInfo[i].second < phasableVariantInfo[j].second);
							assert(phasableVariantInfo[j].second < phasableVariantInfo[k].second);
							if (phasableVariantInfo[j].second < phasableVariantInfo[i].second + minDistance) continue;
							if (phasableVariantInfo[k].second < phasableVariantInfo[j].second + minDistance) continue;
							std::set<std::tuple<size_t, size_t, size_t>> alleles;
							for (const auto& block : distinctBlocks)
							{
								alleles.emplace(block[offi], block[10+offj], block[20+offk]);
							}
							std::map<std::tuple<size_t, size_t, size_t>, std::tuple<size_t, size_t, size_t>> parent;
							for (auto t : alleles)
							{
								if (parent.count(t) == 0) parent[t] = t;
								for (auto t2 : alleles)
								{
									if (t == t2) continue;
									size_t matches = 0;
									if (std::get<0>(t) == std::get<0>(t2)) matches += 1;
									if (std::get<1>(t) == std::get<1>(t2)) matches += 1;
									if (std::get<2>(t) == std::get<2>(t2)) matches += 1;
									if (matches < 2) continue;
									while (parent[parent[t]] != parent[t]) parent[t] = parent[parent[t]];
									while (parent[parent[t2]] != parent[t2]) parent[t2] = parent[parent[t2]];
									parent[parent[t2]] = parent[t];
								}
							}
							std::map<std::tuple<size_t, size_t, size_t>, size_t> newClusters;
							for (auto t : alleles)
							{
								while (parent[parent[t]] != parent[t]) parent[t] = parent[parent[t]];
								if (newClusters.count(parent[t]) == 1) continue;
								size_t cluster = newClusters.size();
								newClusters[parent[t]] = cluster;
							}
							if (newClusters.size() < 2) continue;
							assert(i < phasableVariantInfo.size());
							assert(j < phasableVariantInfo.size());
							assert(k < phasableVariantInfo.size());
							std::vector<std::vector<size_t>> result;
							result.resize(newClusters.size());
							for (size_t m = 0; m < allelesPerRead.size(); m++)
							{
								result[newClusters.at(parent.at(std::make_tuple(allelesPerRead[m][i], allelesPerRead[m][j], allelesPerRead[m][k])))].emplace_back(m);
							}
							bool allClustersCovered = true;
							for (size_t m = 0; m < result.size(); m++)
							{
								if (result[m].size() < minCoverage) allClustersCovered = false;
							}
							if (allClustersCovered)
							{
								std::cerr << "split three sites " << phasableVariantInfo[i].second << " " << phasableVariantInfo[j].second << " " << phasableVariantInfo[k].second << std::endl;
								return result;
							}
						}
					}
				}
			}
		}
	}
	return std::vector<std::vector<size_t>> {};
}

void splitAndAddRecursively(std::vector<std::vector<size_t>>& result, const std::vector<OntLoop>& previous, const std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>>& phasableVariantInfo, const std::vector<size_t>& indicesHere)
{
	std::cerr << "recurse phase of cluster with size " << previous.size() << std::endl;
	std::cerr << "count sites " << phasableVariantInfo.size() << std::endl;
	std::cerr << "try split two sites" << std::endl;
	auto startTime = getTime();
	std::vector<std::vector<size_t>> splitted = trySplitTwoSites(phasableVariantInfo, previous.size());
	auto endTime = getTime();
	std::cerr << "two sites took " << formatTime(startTime, endTime) << std::endl;
	if (splitted.size() >= 2)
	{
		std::cerr << "two sites splitted to clusters of sizes";
		for (size_t i = 0; i < splitted.size(); i++)
		{
			std::cerr << " " << splitted[i].size();
		}
		std::cerr << std::endl;
		std::vector<std::vector<OntLoop>> splittedLoops;
		std::vector<std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>>> splittedPhasableVariantInfo;
		std::tie(splittedLoops, splittedPhasableVariantInfo) = splitLoopsAndPhaseInfoToClusters(previous, phasableVariantInfo, splitted);
		for (size_t i = 0; i < splitted.size(); i++)
		{
			std::vector<size_t> splittedIndices;
			for (size_t j : splitted[i]) splittedIndices.emplace_back(indicesHere[j]);
			splitAndAddRecursively(result, splittedLoops[i], splittedPhasableVariantInfo[i], splittedIndices);
		}
		return;
	}
	std::cerr << "try split three sites" << std::endl;
	startTime = getTime();
	splitted = trySplitThreeSites(phasableVariantInfo);
	endTime = getTime();
	std::cerr << "three sites took " << formatTime(startTime, endTime) << std::endl;
	if (splitted.size() >= 2)
	{
		std::cerr << "three sites splitted to clusters of sizes";
		for (size_t i = 0; i < splitted.size(); i++)
		{
			std::cerr << " " << splitted[i].size();
		}
		std::cerr << std::endl;
		std::vector<std::vector<OntLoop>> splittedLoops;
		std::vector<std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>>> splittedPhasableVariantInfo;
		std::tie(splittedLoops, splittedPhasableVariantInfo) = splitLoopsAndPhaseInfoToClusters(previous, phasableVariantInfo, splitted);
		for (size_t i = 0; i < splitted.size(); i++)
		{
			std::vector<size_t> splittedIndices;
			for (size_t j : splitted[i]) splittedIndices.emplace_back(indicesHere[j]);
			splitAndAddRecursively(result, splittedLoops[i], splittedPhasableVariantInfo[i], splittedIndices);
		}
		return;
	}
	std::cerr << "add cluster of size " << previous.size() << std::endl;
	result.emplace_back(indicesHere);
}

void filterUniqueRegionsToSpan(std::vector<std::pair<size_t, size_t>>& uniqueRegionsInRef, const size_t minPos, const size_t maxPos)
{
	for (size_t i = uniqueRegionsInRef.size()-1; i < uniqueRegionsInRef.size(); i--)
	{
		if (uniqueRegionsInRef[i].first < minPos) uniqueRegionsInRef[i].first = minPos;
		if (uniqueRegionsInRef[i].second > maxPos) uniqueRegionsInRef[i].second = maxPos;
		if (uniqueRegionsInRef[i].first >= uniqueRegionsInRef[i].second)
		{
			std::swap(uniqueRegionsInRef[i], uniqueRegionsInRef.back());
			uniqueRegionsInRef.pop_back();
		}
	}
	std::sort(uniqueRegionsInRef.begin(), uniqueRegionsInRef.end());
}

void filterMatchRegionsToUniqueRegions(const std::vector<std::pair<size_t, size_t>>& uniqueRegionsInRef, std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& matchRegionsInReads)
{
	for (size_t i = 0; i < matchRegionsInReads.size(); i++)
	{
		phmap::flat_hash_set<std::tuple<size_t, size_t, size_t>> newResult;
		for (size_t j = matchRegionsInReads[i].size()-1; j < matchRegionsInReads[i].size(); j--)
		{
			for (auto t : uniqueRegionsInRef)
			{
				if (t.first > std::get<1>(matchRegionsInReads[i][j])) continue;
				if (t.second < std::get<0>(matchRegionsInReads[i][j])) continue;
				size_t leftClip = 0;
				size_t rightClip = 0;
				if (t.first > std::get<0>(matchRegionsInReads[i][j])) leftClip = t.first - std::get<0>(matchRegionsInReads[i][j]);
				if (t.second < std::get<1>(matchRegionsInReads[i][j])) rightClip = std::get<1>(matchRegionsInReads[i][j]) - t.second;
				if (leftClip + rightClip < std::get<1>(matchRegionsInReads[i][j]) - std::get<0>(matchRegionsInReads[i][j]))
				{
					newResult.emplace(std::get<0>(matchRegionsInReads[i][j]) + leftClip, std::get<1>(matchRegionsInReads[i][j]) - rightClip, std::get<2>(matchRegionsInReads[i][j]) + leftClip);
				}
			}
		}
		matchRegionsInReads[i].clear();
		matchRegionsInReads[i].insert(matchRegionsInReads[i].end(), newResult.begin(), newResult.end());
		std::sort(matchRegionsInReads[i].begin(), matchRegionsInReads[i].end());
	}
}

std::vector<std::pair<size_t, size_t>> getAllmatchRegions(const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& matchRegionsInReads)
{
	std::vector<std::pair<size_t, size_t>> result;
	std::vector<size_t> indexPerRead;
	indexPerRead.resize(matchRegionsInReads.size(), 0);
	while (true)
	{
		size_t maxMatchStart = 0;
		size_t minMatchEnd = std::numeric_limits<size_t>::max();
		bool shouldBreak = false;
		for (size_t i = 0; i < matchRegionsInReads.size(); i++)
		{
			if (indexPerRead[i] == matchRegionsInReads[i].size())
			{
				shouldBreak = true;
				break;
			}
			assert(indexPerRead[i] < matchRegionsInReads[i].size());
			maxMatchStart = std::max(maxMatchStart, std::get<0>(matchRegionsInReads[i][indexPerRead[i]]));
			minMatchEnd = std::min(minMatchEnd, std::get<1>(matchRegionsInReads[i][indexPerRead[i]]));
		}
		if (shouldBreak) break;
		if (minMatchEnd < maxMatchStart)
		{
			for (size_t i = 0; i < matchRegionsInReads.size(); i++)
			{
				if (indexPerRead[i] == matchRegionsInReads[i].size()) continue;
				while (indexPerRead[i] < matchRegionsInReads[i].size() && std::get<1>(matchRegionsInReads[i][indexPerRead[i]]) < maxMatchStart) indexPerRead[i] += 1;
				if (indexPerRead[i] == matchRegionsInReads[i].size())
				{
					shouldBreak = true;
					break;
				}
			}
			continue;
		}
		if (shouldBreak) break;
		result.emplace_back(maxMatchStart, minMatchEnd);
		for (size_t i = 0; i < matchRegionsInReads.size(); i++)
		{
			if (indexPerRead[i] == matchRegionsInReads[i].size()) continue;
			if (std::get<1>(matchRegionsInReads[i][indexPerRead[i]]) == minMatchEnd) indexPerRead[i] += 1;
		}
	}
	return result;
}

std::vector<std::vector<size_t>> getReadSequenceLengthsBetweenAllMatchRegions(const std::vector<std::pair<size_t, size_t>>& allmatchRegions, const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& matchRegionsInReads)
{
	std::vector<std::vector<size_t>> result;
	result.resize(matchRegionsInReads.size());
	for (size_t read = 0; read < matchRegionsInReads.size(); read++)
	{
		size_t readIndex = 0;
		size_t lastMatchReadPos = std::numeric_limits<size_t>::max();
		for (size_t allmatchIndex = 0; allmatchIndex < allmatchRegions.size(); allmatchIndex++)
		{
			while (std::get<1>(matchRegionsInReads[read][readIndex]) < allmatchRegions[allmatchIndex].first)
			{
				readIndex += 1;
				assert(readIndex < matchRegionsInReads[read].size());
			}
			assert(std::get<0>(matchRegionsInReads[read][readIndex]) <= allmatchRegions[allmatchIndex].first);
			assert(std::get<1>(matchRegionsInReads[read][readIndex]) >= allmatchRegions[allmatchIndex].second);
			if (allmatchIndex > 0)
			{
				assert(lastMatchReadPos != std::numeric_limits<size_t>::max());
				result[read].emplace_back(std::get<2>(matchRegionsInReads[read][readIndex]) + allmatchRegions[allmatchIndex].first - std::get<0>(matchRegionsInReads[read][readIndex]) - lastMatchReadPos);
			}
			lastMatchReadPos = std::get<2>(matchRegionsInReads[read][readIndex]) + allmatchRegions[allmatchIndex].second - std::get<0>(matchRegionsInReads[read][readIndex]);
		}
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		assert(result[i].size() == allmatchRegions.size()-1);
	}
	return result;
}

std::vector<std::vector<size_t>> getReadLengthAlleles(const std::vector<std::pair<size_t, size_t>>& allmatchRegions, const std::vector<std::vector<size_t>>& readVariantLengths)
{
	std::vector<std::vector<size_t>> result;
	result.resize(readVariantLengths.size());
	for (size_t i = 0; i+1 < allmatchRegions.size(); i++)
	{
		size_t minSeparation = std::max<size_t>(20, (std::get<0>(allmatchRegions[i+1])-std::get<1>(allmatchRegions[i])) * 0.05);
		std::vector<size_t> lengths;
		for (size_t j = 0; j < readVariantLengths.size(); j++)
		{
			lengths.emplace_back(readVariantLengths[j][i]);
		}
		std::sort(lengths.begin(), lengths.end());
		std::vector<size_t> separators;
		for (size_t i = 1; i < lengths.size(); i++)
		{
			if (lengths[i] > lengths[i-1]+minSeparation)
			{
				separators.emplace_back((lengths[i] + lengths[i-1])/2);
			}
		}
		if (separators.size() == 0) continue;
		separators.emplace_back(lengths.back()+1);
		for (size_t j = 0; j < readVariantLengths.size(); j++)
		{
			size_t index = std::upper_bound(separators.begin(), separators.end(), readVariantLengths[j][i]) - separators.begin();
			assert(index < separators.size());
			result[j].emplace_back(index);
		}
	}
	return result;
}

std::vector<std::vector<size_t>> splitByBigIndels(const std::vector<std::vector<std::tuple<size_t, size_t, std::string>>>& editsPerRead)
{
	const size_t minLengthClusterCoverage = 10;
	if (editsPerRead.size() <= 10)
	{
		std::vector<std::vector<size_t>> result;
		result.emplace_back();
		for (size_t i = 0; i < editsPerRead.size(); i++)
		{
			result.back().emplace_back(i);
		}
		return result;
	}
	std::set<std::pair<size_t, size_t>> closeToIndelSpans;
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (auto t : editsPerRead[i])
		{
			if (isSNP(t)) continue;
			if (std::get<1>(t)-std::get<0>(t) == std::get<2>(t).size()) continue;
			closeToIndelSpans.emplace(std::get<0>(t), std::get<1>(t));
		}
	}
	std::vector<std::pair<size_t, size_t>> closeToIndelSpansVec { closeToIndelSpans.begin(), closeToIndelSpans.end() };
	std::sort(closeToIndelSpansVec.begin(), closeToIndelSpansVec.end());
	std::vector<std::pair<size_t, size_t>> nonIndelSpans;
	size_t currentIndelSpanEnd = 0;
	for (auto pair : closeToIndelSpansVec)
	{
		if (pair.first > currentIndelSpanEnd + 10)
		{
			nonIndelSpans.emplace_back(currentIndelSpanEnd, pair.first);
		}
		currentIndelSpanEnd = std::max(currentIndelSpanEnd, pair.second);
	}
	nonIndelSpans.emplace_back(currentIndelSpanEnd, currentIndelSpanEnd+2);
	assert(nonIndelSpans.size() >= 1);
	std::vector<size_t> indelSpanSeparators;
	for (size_t i = 1; i < nonIndelSpans.size(); i++)
	{
		indelSpanSeparators.emplace_back((nonIndelSpans[i].first + nonIndelSpans[i].second)/2);
	}
	std::vector<std::vector<int>> indelLengthSums;
	indelLengthSums.resize(editsPerRead.size());
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		indelLengthSums[i].resize(indelSpanSeparators.size(), 0);
		for (auto t : editsPerRead[i])
		{
			if (isSNP(t)) continue;
			if (std::get<1>(t)-std::get<0>(t) == std::get<2>(t).size()) continue;
			size_t index = std::upper_bound(indelSpanSeparators.begin(), indelSpanSeparators.end(), std::get<0>(t)) - indelSpanSeparators.begin();
			assert(index < indelLengthSums[i].size());
			indelLengthSums[i][index] += std::get<2>(t).size();
			indelLengthSums[i][index] -= std::get<1>(t);
			indelLengthSums[i][index] += std::get<0>(t);
		}
	}
	std::vector<std::vector<size_t>> indelAllelesPerRead;
	indelAllelesPerRead.resize(editsPerRead.size());
	for (size_t j = 0; j < indelSpanSeparators.size(); j++)
	{
		std::vector<int> indelsPerReads;
		for (size_t i = 0; i < editsPerRead.size(); i++)
		{
			indelsPerReads.emplace_back(indelLengthSums[i][j]);
		}
		std::sort(indelsPerReads.begin(), indelsPerReads.end());
		assert(j+1 < nonIndelSpans.size());
		assert(nonIndelSpans[j+1].first >= nonIndelSpans[j].second);
		std::vector<int> separators;
		int requiredLengthDifference = std::max<int>(15, (nonIndelSpans[j+1].first - nonIndelSpans[j].second) * 0.05);
		if (indelsPerReads.back() - indelsPerReads[0] < requiredLengthDifference) continue;
		bool hasSmallSeparator = false;
		for (size_t i = 1; i < indelsPerReads.size(); i++)
		{
			if (indelsPerReads[i] < indelsPerReads[i-1] + requiredLengthDifference) continue;
			size_t lastSeparatorIndex = 0;
			if (separators.size() > 0) lastSeparatorIndex = separators.back();
			separators.emplace_back((indelsPerReads[i] + indelsPerReads[i-1]) / 2);
			if (i - lastSeparatorIndex < minLengthClusterCoverage || i+minLengthClusterCoverage > indelsPerReads.size())
			{
				hasSmallSeparator = true;
				break;
			}
			for (size_t j = lastSeparatorIndex+1; j < i; j++)
			{
				if (indelsPerReads[j] > indelsPerReads[j-1]+10)
				{
					hasSmallSeparator = true;
					break;
				}
			}
			if (hasSmallSeparator) break;
		}
		if (hasSmallSeparator) continue;
		if (separators.size() == 0) continue;
		for (size_t i = 0; i < editsPerRead.size(); i++)
		{
			size_t allele = std::upper_bound(separators.begin(), separators.end(), indelLengthSums[i][j]) - separators.begin();
			assert(allele <= separators.size());
			indelAllelesPerRead[i].emplace_back(allele);
		}
	}
	std::cerr << "indel sites " << indelAllelesPerRead[0].size() << std::endl;
	std::map<std::vector<size_t>, size_t> allelesToCluster;
	std::vector<std::vector<size_t>> result;
	for (size_t i = 0; i < indelAllelesPerRead.size(); i++)
	{
		if (allelesToCluster.count(indelAllelesPerRead[i]) == 0)
		{
			allelesToCluster[indelAllelesPerRead[i]] = result.size();
			result.emplace_back();
		}
		result[allelesToCluster.at(indelAllelesPerRead[i])].emplace_back(i);
	}
	return result;
}

bool minorAllelesPerfectlyLinked(const phmap::flat_hash_map<std::tuple<size_t, size_t, std::string>, std::vector<size_t>>& readsWithEdit, const std::tuple<size_t, size_t, std::string>& editLeft, const std::tuple<size_t, size_t, std::string>& editRight, const size_t countReads)
{
	const size_t minCoverage = 10;
	std::vector<std::pair<bool, bool>> allelesPerRead;
	allelesPerRead.resize(countReads);
	for (size_t read : readsWithEdit.at(editLeft))
	{
		allelesPerRead[read].first = true;
	}
	for (size_t read : readsWithEdit.at(editRight))
	{
		allelesPerRead[read].second = true;
	}
	std::map<std::pair<bool, bool>, size_t> alleleCounts;
	for (size_t i = 0; i < allelesPerRead.size(); i++)
	{
		alleleCounts[allelesPerRead[i]] += 1;
	}
	bool hasSmallAllele = false;
	for (auto pair : alleleCounts)
	{
		if (pair.second < minCoverage) hasSmallAllele = true;
	}
	if (hasSmallAllele) return false;
	assert(alleleCounts.size() >= 2);
	assert(alleleCounts.size() <= 4);
	if (alleleCounts.size() == 4) return false;
	if (alleleCounts.size() == 2)
	{
		for (auto t1 : alleleCounts)
		{
			for (auto t2 : alleleCounts)
			{
				if (t1.first == t2.first) continue;
				size_t edits = 0;
				if (t1.first.first != t2.first.first) edits += 1;
				if (t1.first.second != t2.first.second) edits += 1;
				assert(edits == 1 || edits == 2);
				if (edits == 2) std::cerr << "minor alleles linked 2 edits" << std::endl;
				return edits == 2;
			}
		}
	}
	assert(alleleCounts.size() == 3);
	bool leftIsMinor = readsWithEdit.at(editLeft).size()*2 < countReads;
	bool rightIsMinor = readsWithEdit.at(editRight).size()*2 < countReads;
	if (leftIsMinor ^ rightIsMinor)
	{
		if (alleleCounts.count(std::make_pair(false, false)) == 0) return true;
		if (alleleCounts.count(std::make_pair(true, true)) == 0) return true;
		return false;
	}
	if (alleleCounts.count(std::make_pair(true, false)) == 0) return true;
	if (alleleCounts.count(std::make_pair(false, true)) == 0) return true;
	return false;
}

std::vector<std::vector<size_t>> splitByLinkedMinorAlleles(const std::vector<std::vector<std::tuple<size_t, size_t, std::string>>>& editsPerRead)
{
	const size_t minCoverage = std::max<size_t>(10, editsPerRead.size()*0.1);
	const size_t minDistance = 100;
	if (editsPerRead.size() < 2*minCoverage)
	{
		std::vector<std::vector<size_t>> result;
		result.emplace_back();
		for (size_t i = 0; i < editsPerRead.size(); i++)
		{
			result.back().emplace_back(i);
		}
		return result;
	}
	phmap::flat_hash_map<std::tuple<size_t, size_t, std::string>, std::vector<size_t>> readsWithEdit;
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (const auto& t : editsPerRead[i])
		{
			if (!isSNP(t) && !isBigIndel(t)) continue;
			readsWithEdit[t].emplace_back(i);
		}
	}
	std::vector<std::tuple<size_t, size_t, std::string>> coveredEdits;
	for (const auto& pair : readsWithEdit)
	{
		if (pair.second.size() < minCoverage) continue;
		if (pair.second.size()+minCoverage > editsPerRead.size()) continue;
		coveredEdits.emplace_back(pair.first);
	}
	std::sort(coveredEdits.begin(), coveredEdits.end());
	std::vector<std::vector<size_t>> allelesPerRead;
	allelesPerRead.resize(editsPerRead.size());
	for (size_t i = 0; i < coveredEdits.size(); i++)
	{
		for (size_t j = i+1; j < coveredEdits.size(); j++)
		{
			if (std::get<0>(coveredEdits[j]) < std::get<1>(coveredEdits[i])+minDistance) continue;
			if (!minorAllelesPerfectlyLinked(readsWithEdit, coveredEdits[i], coveredEdits[j], editsPerRead.size())) continue;
			std::cerr << "minor alleles perfectly linked poses " << std::get<0>(coveredEdits[i]) << "-" << std::get<1>(coveredEdits[i]) << " " << std::get<0>(coveredEdits[j]) << "-" << std::get<1>(coveredEdits[j]) << std::endl;
			std::vector<bool> hasFirstEdit;
			hasFirstEdit.resize(editsPerRead.size(), false);
			for (size_t read : readsWithEdit.at(coveredEdits[i]))
			{
				hasFirstEdit[read] = true;
				allelesPerRead[read].emplace_back(1);
			}
			for (size_t read = 0; read < hasFirstEdit.size(); read++)
			{
				if (hasFirstEdit[read]) continue;
				allelesPerRead[read].emplace_back(0);
			}
			std::vector<bool> hasSecondEdit;
			hasSecondEdit.resize(editsPerRead.size(), false);
			for (size_t read : readsWithEdit.at(coveredEdits[j]))
			{
				hasSecondEdit[read] = true;
				allelesPerRead[read].emplace_back(1);
			}
			for (size_t read = 0; read < hasSecondEdit.size(); read++)
			{
				if (hasSecondEdit[read]) continue;
				allelesPerRead[read].emplace_back(0);
			}
			break;
		}
		if (allelesPerRead[0].size() > 0) break;
	}
	std::map<std::vector<size_t>, size_t> allelesToCluster;
	std::vector<std::vector<size_t>> result;
	for (size_t i = 0; i < allelesPerRead.size(); i++)
	{
		if (allelesToCluster.count(allelesPerRead[i]) == 0)
		{
			size_t cluster = result.size();
			allelesToCluster[allelesPerRead[i]] = cluster;
			result.emplace_back();
		}
		result[allelesToCluster.at(allelesPerRead[i])].emplace_back(i);
	}
	return result;
}

std::vector<std::vector<size_t>> splitByEditLinkage(const std::vector<std::vector<std::tuple<size_t, size_t, std::string>>>& editsPerRead)
{
	const size_t minEditsWithinCluster = 5;
	const size_t minReadsWithinCluster = 5;
	const size_t minDistanceBetweenEdits = 100;
	phmap::flat_hash_map<std::tuple<size_t, size_t, std::string>, size_t> editToIndex;
	std::vector<std::vector<size_t>> readsWithEdit;
	std::vector<size_t> editPosition;
	for (size_t read = 0; read < editsPerRead.size(); read++)
	{
		for (const auto& edit : editsPerRead[read])
		{
			if (editToIndex.count(edit) == 0)
			{
				editToIndex[edit] = readsWithEdit.size();
				readsWithEdit.emplace_back();
				readsWithEdit.emplace_back();
				editPosition.emplace_back((std::get<0>(edit) + std::get<1>(edit)) / 2);
			}
			readsWithEdit[editToIndex.at(edit)].emplace_back(read);
		}
	}
	assert(readsWithEdit.size() % 2 == 0);
	for (size_t i = 0; i < readsWithEdit.size(); i += 2)
	{
		std::vector<bool> hasEdit;
		hasEdit.resize(editsPerRead.size(), false);
		for (size_t read : readsWithEdit[i])
		{
			hasEdit[read] = true;
		}
		for (size_t read = 0; read < hasEdit.size(); read++)
		{
			if (hasEdit[read]) continue;
			readsWithEdit[i+1].emplace_back(read);
		}
	}
	for (size_t i = 0; i < readsWithEdit.size(); i++)
	{
		std::sort(readsWithEdit[i].begin(), readsWithEdit[i].end());
	}
	std::vector<size_t> order;
	for (size_t i = 0; i < readsWithEdit.size(); i++)
	{
		order.emplace_back(i);
	}
	std::sort(order.begin(), order.end(), [&readsWithEdit](const size_t left, const size_t right)
	{
		return readsWithEdit[left] < readsWithEdit[right];
	});
	size_t clusterStart = 0;
	std::vector<std::vector<size_t>> readAlleles;
	readAlleles.resize(editsPerRead.size());
	for (size_t i = 1; i <= order.size(); i++)
	{
		if (i < order.size() && readsWithEdit[order[i-1]] == readsWithEdit[order[i]]) continue;
		if (i-clusterStart >= minEditsWithinCluster)
		{
			if (readsWithEdit[order[clusterStart]].size() >= minReadsWithinCluster && readsWithEdit[order[clusterStart]].size()+minReadsWithinCluster <= editsPerRead.size())
			{
				std::vector<size_t> positions;
				for (size_t j = clusterStart; j < i; j++)
				{
					positions.emplace_back(editPosition[order[j]]);
				}
				std::sort(positions.begin(), positions.end());
				size_t reallyDifferentLocations = 1;
				for (size_t j = 1; j < positions.size(); j++)
				{
					if (positions[j] < positions[j-1] + minDistanceBetweenEdits) continue;
					reallyDifferentLocations += 1;
				}
				if (reallyDifferentLocations >= minEditsWithinCluster)
				{
					std::vector<bool> hasAlt;
					hasAlt.resize(editsPerRead.size(), false);
					for (size_t read : readsWithEdit[order[clusterStart]])
					{
						hasAlt[read] = true;
						readAlleles[read].emplace_back(0);
					}
					for (size_t j = 0; j < editsPerRead.size(); j++)
					{
						if (hasAlt[j]) continue;
						readAlleles[j].emplace_back(1);
					}
				}
			}
		}
		clusterStart = i;
	}
	std::cerr << "edit linkage clusters " << readAlleles[0].size() << std::endl;
	std::map<std::vector<size_t>, size_t> allelesToCluster;
	std::vector<std::vector<size_t>> result;
	for (size_t i = 0; i < readAlleles.size(); i++)
	{
		if (allelesToCluster.count(readAlleles[i]) == 0)
		{
			allelesToCluster[readAlleles[i]] = result.size();
			result.emplace_back();
		}
		result[allelesToCluster.at(readAlleles[i])].emplace_back(i);
	}
	return result;
}

std::pair<size_t, size_t> getBubble(const GfaGraph& graph, const Node startNode)
{
	const std::pair<size_t, size_t> noBubble { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
	// std::numeric_limits<size_t>::max() cannot be used for undefined because Node can't contain it
	// break convention and use graph.numNodes() as undefined
	const Node undefinedNode { graph.numNodes(), true };
	if (graph.edges.count(startNode) == 0) return noBubble;
	if (graph.edges.at(startNode).size() != 2) return noBubble;
	Node otherSide = undefinedNode;
	size_t alleleOne = std::numeric_limits<size_t>::max();
	size_t alleleTwo = std::numeric_limits<size_t>::max();
	for (auto edge : graph.edges.at(startNode))
	{
		Node target = std::get<0>(edge);
		if (graph.edges.count(target) == 0) return noBubble;
		if (graph.edges.at(target).size() != 1) return noBubble;
		assert(graph.edges.count(reverse(target)) == 1);
		if (graph.edges.at(reverse(target)).size() != 1) return noBubble;
		if (otherSide == undefinedNode)
		{
			assert(alleleOne == std::numeric_limits<size_t>::max());
			alleleOne = target.id();
			otherSide = std::get<0>(*graph.edges.at(target).begin());
		}
		else
		{
			assert(alleleOne != std::numeric_limits<size_t>::max());
			assert(alleleTwo == std::numeric_limits<size_t>::max());
			alleleTwo = target.id();
			if (std::get<0>(*graph.edges.at(target).begin()) != otherSide) return noBubble;
		}
	}
	assert(otherSide != undefinedNode);
	assert(graph.edges.count(reverse(otherSide)) == 1);
	if (graph.edges.at(reverse(otherSide)).size() != 2) return noBubble;
	if (alleleOne == startNode.id()) return noBubble;
	if (alleleOne == otherSide.id()) return noBubble;
	if (alleleTwo == startNode.id()) return noBubble;
	if (alleleTwo == otherSide.id()) return noBubble;
	if (alleleOne == alleleTwo) return noBubble;
	return std::make_pair(alleleOne, alleleTwo);
}

std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>> getDBGvariants(const std::vector<OntLoop>& cluster, const std::string MBGPath, const std::string tmpPath)
{
	size_t minCoverage = 5;
	std::string graphFile = tmpPath + "/tmpgraph.gfa";
	std::string pathsFile = tmpPath + "/tmppaths.gaf";
	std::string readsFile = tmpPath + "/tmpreads.fa";
	{
		std::ofstream file { readsFile };
		for (size_t i = 0; i < cluster.size(); i++)
		{
			file << ">" << i << std::endl;
			file << cluster[i].rawSequence << std::endl;
		}
	}
	std::string mbgCommand = MBGPath + " -i " + readsFile + " -o " + graphFile + " --output-sequence-paths " + pathsFile + " --error-masking=no -a 5 -u 5 --no-multiplex-cleaning -r 18446744073709551615 -R 18446744073709551615 -k 31 --only-local-resolve --do-unsafe-guesswork-resolutions 1> " + tmpPath + "/tmpstdout.txt 2> " + tmpPath + "/tmpstderr.txt";
	std::cerr << "MBG command:" << std::endl;
	std::cerr << mbgCommand << std::endl;
	int runresult = system(mbgCommand.c_str());
	if (runresult != 0)
	{
		std::cerr << "MBG did not run successfully" << std::endl;
		std::abort();
	}
	GfaGraph graph;
	graph.loadFromFile(graphFile);
	std::vector<ReadPath> readPaths = loadReadPaths(pathsFile, graph);
	std::vector<phmap::flat_hash_set<size_t>> readsTouchingNode;
	readsTouchingNode.resize(graph.numNodes());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		size_t read = std::stoull(readPaths[i].readName);
		for (Node node : readPaths[i].path)
		{
			readsTouchingNode[node.id()].insert(read);
		}
	}
	phmap::flat_hash_map<size_t, std::pair<size_t, size_t>> nodeBelongsToBubble;
	size_t countBubbles = 0;
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (readsTouchingNode[i].size() != cluster.size()) continue;
		std::pair<size_t, size_t> bubble = getBubble(graph, Node { i, true });
		if (bubble.first != std::numeric_limits<size_t>::max())
		{
			assert(bubble.second != std::numeric_limits<size_t>::max());
			assert(bubble.first != bubble.second);
			if (nodeBelongsToBubble.count(std::min(bubble.first, bubble.second)) == 1)
			{
				assert(nodeBelongsToBubble.count(std::max(bubble.first, bubble.second)) == 1);
			}
			else
			{
				nodeBelongsToBubble[std::min(bubble.first, bubble.second)] = std::make_pair(countBubbles, 0);
				nodeBelongsToBubble[std::max(bubble.first, bubble.second)] = std::make_pair(countBubbles, 1);
				countBubbles += 1;
			}
		}
		bubble = getBubble(graph, Node { i, false });
		if (bubble.first != std::numeric_limits<size_t>::max())
		{
			assert(bubble.second != std::numeric_limits<size_t>::max());
			assert(bubble.first != bubble.second);
			if (nodeBelongsToBubble.count(std::min(bubble.first, bubble.second)) == 1)
			{
				assert(nodeBelongsToBubble.count(std::max(bubble.first, bubble.second)) == 1);
			}
			else
			{
				nodeBelongsToBubble[std::min(bubble.first, bubble.second)] = std::make_pair(countBubbles, 0);
				nodeBelongsToBubble[std::max(bubble.first, bubble.second)] = std::make_pair(countBubbles, 1);
				countBubbles += 1;
			}
		}
	}
	std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>> result;
	result.resize(countBubbles);
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i].second = i*1000; // fake positions, graph always has good good sites and doesn't falsely split a single variant into multiple edits
		result[i].first.resize(3); // always simple biallelic bubbles plus one for reads with no bubble allele because of sequencing errors
	}
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (nodeBelongsToBubble.count(i) == 0) continue;
		size_t bubbleIndex = nodeBelongsToBubble.at(i).first;
		size_t alleleIndex = nodeBelongsToBubble.at(i).second;
		for (size_t read : readsTouchingNode[i])
		{
			result[bubbleIndex].first[alleleIndex].emplace_back(read);
		}
	}
	for (size_t i = result.size()-1; i < result.size(); i--)
	{
		std::vector<bool> hasAllele;
		hasAllele.resize(cluster.size());
		assert(result[i].first.size() == 3);
		assert(result[i].first[2].size() == 0);
		bool good = true;
		if (result[i].first[0].size() < minCoverage) good = false;
		if (result[i].first[1].size() < minCoverage) good = false;
		for (size_t j = 0; j < 2; j++)
		{
			for (size_t read : result[i].first[j])
			{
				if (hasAllele[read]) good = false;
				hasAllele[read] = true;
			}
		}
		for (size_t j = 0; j < hasAllele.size(); j++)
		{
			if (hasAllele[j]) continue;
			result[i].first[2].emplace_back(j);
		}
		if (result[i].first[2].size() >= 3) good = false;
		// equality valid for breaking because truncation rounds down, only allow missing values once coverage is 20
		if (result[i].first[0].size() + result[i].first[1].size() <= cluster.size()*0.95) good = false;
		if (!good)
		{
			std::swap(result[i], result.back());
			result.pop_back();
		}
	}
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) { return left.second < right.second; });
	return result;
}

template <typename F>
void iterateSuccessorKmers(const __uint128_t kmer, const size_t k, F callback)
{
	const __uint128_t mask = ((__uint128_t)1 << (__uint128_t)(2ull*k))-(__uint128_t)1;
	for (size_t i = 0; i < 4; i++)
	{
		callback(((kmer << 2) & mask) + i);
	}
}

template <typename F>
void iteratePredecessorKmers(const __uint128_t kmer, const size_t k, F callback)
{
	for (size_t i = 0; i < 4; i++)
	{
		callback((kmer >> 2) + ((__uint128_t)i << (__uint128_t)(2*(k-1))));
	}
}

void addPath(std::string& result, const phmap::flat_hash_set<__uint128_t>& keptKmers, const size_t lastMatchPos, const size_t thisMatchPos, const std::string& rawSequence, const size_t k)
{
	const __uint128_t mask = ((__uint128_t)1 << (__uint128_t)(2ull*k))-(__uint128_t)1;
	__uint128_t lastMatchKmer = 0;
	__uint128_t thisMatchKmer = 0;
	for (size_t i = 0; i < k; i++)
	{
		lastMatchKmer <<= 2;
		thisMatchKmer <<= 2;
		switch(rawSequence[lastMatchPos+i])
		{
			case 'A':
				lastMatchKmer += 0;
				break;
			case 'C':
				lastMatchKmer += 1;
				break;
			case 'G':
				lastMatchKmer += 2;
				break;
			case 'T':
				lastMatchKmer += 3;
				break;
			default:
				assert(false);
		}
		switch(rawSequence[thisMatchPos+i])
		{
			case 'A':
				thisMatchKmer += 0;
				break;
			case 'C':
				thisMatchKmer += 1;
				break;
			case 'G':
				thisMatchKmer += 2;
				break;
			case 'T':
				thisMatchKmer += 3;
				break;
			default:
				assert(false);
		}
	}
	if (lastMatchKmer == thisMatchKmer)
	{
		result += rawSequence.substr(lastMatchPos + k, thisMatchPos - lastMatchPos);
		return;
	}
	assert(keptKmers.count(lastMatchKmer) == 1);
	assert(keptKmers.count(thisMatchKmer) == 1);
	std::string unambiguousExtensionFromLast;
	__uint128_t wanderKmer = lastMatchKmer;
	phmap::flat_hash_set<__uint128_t> visited;
	while (wanderKmer != thisMatchKmer)
	{
		size_t uniqueSuccessor = 4;
		for (size_t i = 0; i < 4; i++)
		{
			__uint128_t testKmer = ((wanderKmer << (__uint128_t)2) & mask) + i;
			if (keptKmers.count(testKmer) == 0) continue;
			if (uniqueSuccessor == 4)
			{
				uniqueSuccessor = i;
			}
			else
			{
				uniqueSuccessor = 5;
			}
		}
		if (uniqueSuccessor >= 4)
		{
			break;
		}
		unambiguousExtensionFromLast += "ACGT"[uniqueSuccessor];
		size_t nextKmer = ((wanderKmer << 2) & mask) + uniqueSuccessor;
		wanderKmer = nextKmer;
		if (visited.count(wanderKmer) == 1) break;
		visited.insert(wanderKmer);
	}
	if (wanderKmer == thisMatchKmer)
	{
		assert(unambiguousExtensionFromLast.size() >= 1);
		size_t pathLength = visited.size();
		if (pathLength + 50 > thisMatchPos-lastMatchPos)
		{
			if (thisMatchPos-lastMatchPos+50 > pathLength)
			{
				result += unambiguousExtensionFromLast;
				return;
			}
		}
	}
	visited.clear();
	wanderKmer = thisMatchKmer;
	unambiguousExtensionFromLast.clear();
	while (wanderKmer != lastMatchKmer)
	{
		size_t uniquePredecessor = 4;
		for (size_t i = 0; i < 4; i++)
		{
			__uint128_t testKmer = ((wanderKmer >> 2)) + (i << (2ull*(k-1ull)));
			if (keptKmers.count(testKmer) == 0) continue;
			if (uniquePredecessor == 4)
			{
				uniquePredecessor = i;
			}
			else
			{
				uniquePredecessor = 5;
			}
		}
		if (uniquePredecessor >= 4)
		{
			break;
		}
		unambiguousExtensionFromLast += "ACGT"[wanderKmer & 3];
		size_t nextKmer = ((wanderKmer >> 2)) + (uniquePredecessor << (2ull*(k-1ull)));
		wanderKmer = nextKmer;
		if (visited.count(wanderKmer) == 1) break;
		visited.insert(wanderKmer);
	}
	if (wanderKmer == lastMatchKmer)
	{
		assert(unambiguousExtensionFromLast.size() >= 1);
		size_t pathLength = visited.size();
		if (pathLength + 20 > thisMatchPos-lastMatchPos)
		{
			if (thisMatchPos-lastMatchPos+20 > pathLength)
			{
				std::reverse(unambiguousExtensionFromLast.begin(), unambiguousExtensionFromLast.end());
				result += unambiguousExtensionFromLast;
				return;
			}
		}
	}
	result += rawSequence.substr(lastMatchPos + k, thisMatchPos - lastMatchPos);
	return;
}

bool isTipFw(const __uint128_t kmer, const phmap::flat_hash_set<__uint128_t>& keptKmers, const phmap::flat_hash_map<__uint128_t, size_t>& kmerCounts, const size_t k)
{
	size_t countNeighbors = 0;
	iterateSuccessorKmers(kmer, k, [&keptKmers, &countNeighbors](const __uint128_t neighbor)
	{
		if (keptKmers.count(neighbor) == 1) countNeighbors += 1;
	});
	return countNeighbors == 0;
}

bool isTipBw(const __uint128_t kmer, const phmap::flat_hash_set<__uint128_t>& keptKmers, const phmap::flat_hash_map<__uint128_t, size_t>& kmerCounts, const size_t k)
{
	size_t countNeighbors = 0;
	iteratePredecessorKmers(kmer, k, [&keptKmers, &countNeighbors](const __uint128_t neighbor)
	{
		if (keptKmers.count(neighbor) == 1) countNeighbors += 1;
	});
	return countNeighbors == 0;
}

void checkCanRemoveBw(phmap::flat_hash_set<__uint128_t>& removeThese, const __uint128_t kmer, const phmap::flat_hash_set<__uint128_t>& keptKmers, const phmap::flat_hash_map<__uint128_t, size_t>& kmerCounts, const size_t k)
{
	phmap::flat_hash_set<__uint128_t> visited;
	visited.insert(kmer);
	__uint128_t wanderKmer = kmer;
	for (size_t iteration = 0; iteration < k; iteration++)
	{
		std::vector<__uint128_t> allPredecessors;
		iteratePredecessorKmers(wanderKmer, k, [&keptKmers, &allPredecessors](const __uint128_t testKmer)
		{
			if (keptKmers.count(testKmer) == 0) return;
			allPredecessors.emplace_back(testKmer);
		});
		if (allPredecessors.size() == 0) break;
		if (allPredecessors.size() >= 2)
		{
			bool allPredecessorsHaveGoodAlternate = true;
			for (__uint128_t predecessor : allPredecessors)
			{
				bool thisPredecessorHasGoodAlternate = false;
				iterateSuccessorKmers(predecessor, k, [&keptKmers, &kmerCounts, &thisPredecessorHasGoodAlternate](const __uint128_t neighbor)
				{
					if (keptKmers.count(neighbor) == 0) return;
					if (kmerCounts.at(neighbor) >= 3) thisPredecessorHasGoodAlternate = true;
				});
				if (!thisPredecessorHasGoodAlternate) allPredecessorsHaveGoodAlternate = false;
			}
			if (!allPredecessorsHaveGoodAlternate) return;
			break;
		}
		assert(allPredecessors.size() == 1);
		bool predecessorHasOtherNeighbor = false;
		bool predecessorHasBetterNeigbor = false;
		iterateSuccessorKmers(allPredecessors[0], k, [&keptKmers, &kmerCounts, &predecessorHasOtherNeighbor, &predecessorHasBetterNeigbor, wanderKmer](const __uint128_t neighbor)
		{
			if (keptKmers.count(neighbor) == 0) return;
			if (neighbor == wanderKmer) return;
			predecessorHasOtherNeighbor = true;
			if (kmerCounts.at(neighbor) >= 3) predecessorHasBetterNeigbor = true;
		});
		if (predecessorHasOtherNeighbor && predecessorHasBetterNeigbor) break;
		if (predecessorHasOtherNeighbor && !predecessorHasBetterNeigbor) return;
		if (kmerCounts.at(allPredecessors[0]) >= 3) return;
		wanderKmer = allPredecessors[0];
		visited.insert(wanderKmer);
	}
	assert(visited.size() >= 1);
	for (__uint128_t kmer : visited)
	{
		assert(keptKmers.count(kmer) == 1);
		assert(kmerCounts.at(kmer) == 2);
		removeThese.insert(kmer);
	}
}

void checkCanRemoveFw(phmap::flat_hash_set<__uint128_t>& removeThese, const __uint128_t kmer, const phmap::flat_hash_set<__uint128_t>& keptKmers, const phmap::flat_hash_map<__uint128_t, size_t>& kmerCounts, const size_t k)
{
	phmap::flat_hash_set<__uint128_t> visited;
	visited.insert(kmer);
	__uint128_t wanderKmer = kmer;
	for (size_t iteration = 0; iteration < k; iteration++)
	{
		std::vector<__uint128_t> allSuccessors;
		iterateSuccessorKmers(wanderKmer, k, [&keptKmers, &allSuccessors](const __uint128_t testKmer)
		{
			if (keptKmers.count(testKmer) == 0) return;
			allSuccessors.emplace_back(testKmer);
		});
		if (allSuccessors.size() == 0) break;
		if (allSuccessors.size() >= 2)
		{
			bool allSuccessorsHaveGoodAlternate = true;
			for (__uint128_t successor : allSuccessors)
			{
				bool thisSuccessorHasGoodAlternate = false;
				iteratePredecessorKmers(successor, k, [&keptKmers, &kmerCounts, &thisSuccessorHasGoodAlternate](const __uint128_t neighbor)
				{
					if (keptKmers.count(neighbor) == 0) return;
					if (kmerCounts.at(neighbor) >= 3) thisSuccessorHasGoodAlternate = true;
				});
				if (!thisSuccessorHasGoodAlternate) allSuccessorsHaveGoodAlternate = false;
			}
			if (!allSuccessorsHaveGoodAlternate) return;
			break;
		}
		assert(allSuccessors.size() == 1);
		bool successorHasOtherNeighbor = false;
		bool successorHasBetterNeigbor = false;
		iteratePredecessorKmers(allSuccessors[0], k, [&keptKmers, &kmerCounts, &successorHasOtherNeighbor, &successorHasBetterNeigbor, wanderKmer](const __uint128_t neighbor)
		{
			if (keptKmers.count(neighbor) == 0) return;
			if (neighbor == wanderKmer) return;
			successorHasOtherNeighbor = true;
			if (kmerCounts.at(neighbor) >= 3) successorHasBetterNeigbor = true;
		});
		if (successorHasOtherNeighbor && successorHasBetterNeigbor) break;
		if (successorHasOtherNeighbor && !successorHasBetterNeigbor) return;
		if (kmerCounts.at(allSuccessors[0]) >= 3) return;
		wanderKmer = allSuccessors[0];
		visited.insert(wanderKmer);
	}
	assert(visited.size() >= 1);
	for (__uint128_t kmer : visited)
	{
		assert(keptKmers.count(kmer) == 1);
		assert(kmerCounts.at(kmer) == 2);
		removeThese.insert(kmer);
	}
}

void removeTips(phmap::flat_hash_set<__uint128_t>& keptKmers, const phmap::flat_hash_map<__uint128_t, size_t>& kmerCounts, const size_t k)
{
	phmap::flat_hash_set<__uint128_t> removeThese;
	for (__uint128_t kmer : keptKmers)
	{
		assert(kmerCounts.count(kmer) == 1);
		assert(kmerCounts.at(kmer) >= 2);
		if (kmerCounts.at(kmer) != 2) continue;
		if (isTipFw(kmer, keptKmers, kmerCounts, k))
		{
			checkCanRemoveBw(removeThese, kmer, keptKmers, kmerCounts, k);
		}
		if (isTipBw(kmer, keptKmers, kmerCounts, k))
		{
			checkCanRemoveFw(removeThese, kmer, keptKmers, kmerCounts, k);
		}
	}
	for (__uint128_t kmer : removeThese)
	{
		assert(keptKmers.count(kmer) == 1);
		keptKmers.erase(kmer);
	}
}

std::vector<std::string> getSelfCorrectedLoopsOneIteration(const std::vector<std::string>& cluster, const size_t k)
{
	phmap::flat_hash_map<__uint128_t, size_t> kmerCounts;
	for (size_t i = 0; i < cluster.size(); i++)
	{
		iterateKmers(cluster[i], k, [&kmerCounts](const __uint128_t kmer, const size_t position) { kmerCounts[kmer] += 1; });
	}
	phmap::flat_hash_set<__uint128_t> keptKmers;
	for (auto pair : kmerCounts)
	{
		if (pair.second == 1) continue;
		keptKmers.insert(pair.first);
	}
	removeTips(keptKmers, kmerCounts, k);
	std::vector<std::string> result;
	for (size_t i = 0; i < cluster.size(); i++)
	{
		std::string resultHere;
		size_t lastKmerMatch = std::numeric_limits<size_t>::max();
		iterateKmers(cluster[i], k, [&keptKmers, &cluster, &resultHere, &lastKmerMatch, i, k](const __uint128_t kmer, const size_t pos)
		{
			if (keptKmers.count(kmer) == 0) return;
			if (lastKmerMatch == std::numeric_limits<size_t>::max())
			{
				resultHere += cluster[i].substr(0, pos+k);
				lastKmerMatch = pos;
				return;
			}
			if (pos == lastKmerMatch+1)
			{
				resultHere += cluster[i][lastKmerMatch+k];
				lastKmerMatch = pos;
				return;
			}
			addPath(resultHere, keptKmers, lastKmerMatch, pos, cluster[i], k);
			lastKmerMatch = pos;
		});
		if (lastKmerMatch == std::numeric_limits<size_t>::max())
		{
			resultHere = cluster[i];
		}
		else
		{
			resultHere += cluster[i].substr(lastKmerMatch+k);
		}
		result.emplace_back(resultHere);
	}
	return result;
}

std::vector<std::string> getSelfCorrectedLoops(const std::vector<OntLoop>& cluster)
{
	std::vector<std::string> raws;
	for (size_t i = 0; i < cluster.size(); i++)
	{
		raws.emplace_back(cluster[i].rawSequence);
	}
	std::vector<std::string> partwiseCorrected;
	partwiseCorrected.resize(cluster.size());
	for (size_t i = 0; i < 10; i++)
	{
		std::vector<std::string> parts;
		for (size_t j = 0; j < cluster.size(); j++)
		{
			parts.emplace_back(raws[j].substr(raws[j].size()*i/10, raws[j].size()*(i+1)/10 - raws[j].size()*i/10));
		}
		std::vector<std::string> correctedParts = getSelfCorrectedLoopsOneIteration(parts, 15);
		correctedParts = getSelfCorrectedLoopsOneIteration(correctedParts, 31);
		correctedParts = getSelfCorrectedLoopsOneIteration(correctedParts, 63);
		for (size_t j = 0; j < cluster.size(); j++)
		{
			partwiseCorrected[j] += correctedParts[j];
		}
	}
	phmap::flat_hash_map<uint64_t, size_t> kmerCounts;
	std::vector<std::string> fullCorrected = getSelfCorrectedLoopsOneIteration(partwiseCorrected, 31);
	fullCorrected = getSelfCorrectedLoopsOneIteration(partwiseCorrected, 63);
	return fullCorrected;
}

void callVariantsAndSplitRecursively(std::vector<std::vector<OntLoop>>& result, const std::vector<OntLoop>& cluster, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t numThreads, const std::string MBGPath, const std::string tmpPath, const bool correctBeforePhasing)
{
	if (cluster.size() <= 5)
	{
		std::cerr << "skip small cluster with size " << cluster.size() << std::endl;
		result.emplace_back(cluster);
		return;
	}
	std::cerr << "begin phasing cluster with size " << cluster.size() << std::endl;
	std::vector<std::string> sequences;
	if (correctBeforePhasing)
	{
		sequences = getSelfCorrectedLoops(cluster);
	}
	else
	{
		for (size_t i = 0; i < cluster.size(); i++)
		{
			sequences.emplace_back(cluster[i].rawSequence);
		}
	}
	auto startTime = getTime();
	auto edits = getEditsForPhasing(cluster, sequences, graph, pathStartClip, pathEndClip, numThreads);
	auto phasableVariantInfo = getPhasableVariantInfoBiallelicAltRefSNPsBigIndels(edits);
	auto endTime = getTime();
	std::cerr << "getting variants took " << formatTime(startTime, endTime) << std::endl;
	std::cerr << "got variant info of cluster with size " << cluster.size() << std::endl;
	std::vector<std::vector<size_t>> resultHere;
	std::vector<size_t> allIndices;
	for (size_t i = 0; i < cluster.size(); i++) allIndices.emplace_back(i);
	splitAndAddRecursively(resultHere, cluster, phasableVariantInfo, allIndices);
	assert(resultHere.size() >= 1);
	if (resultHere.size() == 1)
	{
		std::cerr << "split by bubble alleles size " << cluster.size() << std::endl;
		resultHere = splitByBubbleAlleles(cluster, sequences, graph, pathStartClip, pathEndClip, numThreads);
		std::cerr << "bubble alleles splitted to " << resultHere.size() << " clusters:";
		for (size_t i = 0; i < resultHere.size(); i++)
		{
			std::cerr << " " << resultHere[i].size();
		}
		std::cerr << std::endl;
	}
	if (resultHere.size() == 1)
	{
		std::cerr << "split by SNP MSA size " << cluster.size() << std::endl;
		auto SNPMSA = getSNPMSA(edits);
		std::cerr << SNPMSA.size() << " sites in SNP MSA" << std::endl;
		resultHere.clear();
		splitAndAddRecursively(resultHere, cluster, SNPMSA, allIndices);
		if (resultHere.size() > 1)
		{
			std::cerr << "SNP MSA splitted to " << resultHere.size() << " clusters:";
			for (size_t i = 0; i < resultHere.size(); i++)
			{
				std::cerr << " " << resultHere[i].size();
			}
			std::cerr << std::endl;
		}
	}
	if (resultHere.size() == 1)
	{
		std::cerr << "split by DBG variants size " << cluster.size() << std::endl;
		auto DBGvariants = getDBGvariants(cluster, MBGPath, tmpPath);
		std::cerr << DBGvariants.size() << " sites in DBG variants" << std::endl;
		resultHere.clear();
		splitAndAddRecursively(resultHere, cluster, DBGvariants, allIndices);
		if (resultHere.size() > 1)
		{
			std::cerr << "DBG variants splitted to " << resultHere.size() << " clusters:";
			for (size_t i = 0; i < resultHere.size(); i++)
			{
				std::cerr << " " << resultHere[i].size();
			}
			std::cerr << std::endl;
		}
	}
	if (resultHere.size() == 1)
	{
		std::cerr << "split by big indels size " << cluster.size() << std::endl;
		resultHere = splitByBigIndels(edits);
		std::cerr << "big indels splitted to " << resultHere.size() << " clusters:";
		for (size_t i = 0; i < resultHere.size(); i++)
		{
			std::cerr << " " << resultHere[i].size();
		}
		std::cerr << std::endl;
	}
	if (resultHere.size() == 1)
	{
		std::cerr << "split by edit linkage size " << cluster.size() << std::endl;
		resultHere = splitByEditLinkage(edits);
		std::cerr << "edit linkage splitted to " << resultHere.size() << " clusters:";
		for (size_t i = 0; i < resultHere.size(); i++)
		{
			std::cerr << " " << resultHere[i].size();
		}
		std::cerr << std::endl;
	}
	if (resultHere.size() == 1)
	{
		std::cerr << "split by linked minor alleles size " << cluster.size() << std::endl;
		resultHere = splitByLinkedMinorAlleles(edits);
		std::cerr << "linked minor alleles splitted to " << resultHere.size() << " clusters:";
		for (size_t i = 0; i < resultHere.size(); i++)
		{
			std::cerr << " " << resultHere[i].size();
		}
		std::cerr << std::endl;
	}
	if (resultHere.size() == 1)
	{
		if (correctBeforePhasing)
		{
			std::cerr << "final cluster with size " << cluster.size() << " with reads:";
			for (const auto& read : cluster)
			{
				std::cerr << " " << read.readName << "_" << read.approxStart << "_" << read.approxEnd;
			}
			std::cerr << std::endl;
			result.emplace_back(cluster);
			return;
		}
		else
		{
			callVariantsAndSplitRecursively(result, cluster, graph, pathStartClip, pathEndClip, numThreads, MBGPath, tmpPath, true);
			return;
		}
	}
	std::cerr << "split cluster with size " << cluster.size() << " corrected " << (correctBeforePhasing ? "yes" : "no")  << " into " << resultHere.size() << " clusters with reads:";
	for (const auto& read : cluster)
	{
		std::cerr << " " << read.readName << "_" << read.approxStart << "_" << read.approxEnd;
	}
	std::cerr << std::endl;
	for (size_t i = 0; i < resultHere.size(); i++)
	{
		std::vector<OntLoop> filteredClusters;
		for (size_t j : resultHere[i])
		{
			filteredClusters.emplace_back(cluster[j]);
		}
		callVariantsAndSplitRecursively(result, filteredClusters, graph, pathStartClip, pathEndClip, numThreads, MBGPath, tmpPath, false);
	}
}

void addRawSequencesToLoops(std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::string ontReadPath, const size_t numThreads)
{
	std::vector<std::vector<std::string>> sequences = getRawLoopSequences(clusters, ontReadPath, graph, pathStartClip, pathEndClip, numThreads);
	assert(sequences.size() == clusters.size());
	for (size_t i = 0; i < clusters.size(); i++)
	{
		assert(sequences[i].size() == clusters[i].size());
		for (size_t j = 0; j < clusters[i].size(); j++)
		{
			clusters[i][j].rawSequence = sequences[i][j];
		}
	}
}

std::vector<std::vector<OntLoop>> editSplitClusters(const std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t numThreads, const std::string MBGPath, const std::string tmpPath)
{
	std::vector<std::vector<OntLoop>> result;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		if (clusters[i].size() <= 5)
		{
			std::cerr << "skip low coverage cluster with size " << clusters[i].size() << std::endl;
			result.emplace_back(clusters[i]);
			continue;
		}
		callVariantsAndSplitRecursively(result, clusters[i], graph, pathStartClip, pathEndClip, numThreads, MBGPath, tmpPath, false);
	}
	return result;
}

void addRawSequenceNamesToLoops(std::vector<std::vector<OntLoop>>& clusters)
{
	for (size_t i = 0; i < clusters.size(); i++)
	{
		for (size_t j = 0; j < clusters[i].size(); j++)
		{
			clusters[i][j].rawLoopName = "morph_" + std::to_string(i) + "_rawloop_" + std::to_string(j) + "_" + clusters[i][j].readName + "_" + std::to_string(clusters[i][j].approxStart) + "_" + std::to_string(clusters[i][j].approxEnd);
		}
	}
}

void addSelfCorrectedOntLoopSequences(std::vector<std::vector<OntLoop>>& clusters)
{
	for (size_t i = 0; i < clusters.size(); i++)
	{
		std::cerr << "self-correct cluster " << i << " with size " << clusters[i].size() << std::endl;
		auto selfCorrectedLoops = getSelfCorrectedLoops(clusters[i]);
		assert(selfCorrectedLoops.size() == clusters[i].size());
		for (size_t j = 0; j < clusters[i].size(); j++)
		{
			clusters[i][j].selfCorrectedSequence = selfCorrectedLoops[j];
			std::cerr << "corrected read " << j << " from size " << clusters[i][j].rawSequence.size() << " to " << clusters[i][j].selfCorrectedSequence.size() << std::endl;
		}
	}
}

void DoClusterONTAnalysis(const ClusterParams& params)
{
	std::cerr << "reading allele graph" << std::endl;
	GfaGraph graph;
	graph.loadFromFile(params.basePath + "/processed-graph.gfa");
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
	std::unordered_set<size_t> anchorNodes;
	std::tie(borderNodes, anchorNodes, pathStartClip, pathEndClip) = getBorderNodes(heavyPath, graph);
	assert(borderNodes.size() > 0);
	auto loopSequences = extractLoopSequences(ontPaths, heavyPath, minLength, graph, borderNodes, anchorNodes, pathStartClip, pathEndClip);
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
	auto clusters = roughClusterLoopSequences(loopSequences, graph, pathStartClip, pathEndClip, coreNodes, params.maxClusterDifference);
	std::cerr << clusters.size() << " rough clusters" << std::endl;
	std::cerr << "getting exact locations of raw loop sequences" << std::endl;
	addRawSequencesToLoops(clusters, graph, pathStartClip, pathEndClip, params.ontReadPath, params.numThreads);
	std::cerr << "phase clusters by raw sequences" << std::endl;
	clusters = editSplitClusters(clusters, graph, pathStartClip, pathEndClip, params.numThreads, params.MBGPath, params.basePath + "/tmp");
	std::cerr << clusters.size() << " clusters" << std::endl;
	std::cerr << "cluster loops by density" << std::endl;
	clusters = densityClusterLoops(clusters, graph, pathStartClip, pathEndClip, coreNodes, params.maxClusterDifference, 5, params.minReclusterDistance);
	std::cerr << clusters.size() << " clusters" << std::endl;
	std::cerr << "phase clusters by raw sequences" << std::endl;
	clusters = editSplitClusters(clusters, graph, pathStartClip, pathEndClip, params.numThreads, params.MBGPath, params.basePath + "/tmp");
	std::cerr << clusters.size() << " clusters" << std::endl;
	// std::cerr << "phase clusters" << std::endl;
	// clusters = phaseClusters(clusters);
	// std::cerr << clusters.size() << " phased clusters" << std::endl;
	std::cerr << clusters.size() << " phased clusters" << std::endl;
	std::sort(clusters.begin(), clusters.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
	addRawSequenceNamesToLoops(clusters);
	std::cerr << "write raw ONT loop sequences" << std::endl;
	writeRawOntLoopSequences(params.basePath + "/raw_loops.fa", clusters);
	std::cerr << "get self-corrected ONT loop sequences" << std::endl;
	addSelfCorrectedOntLoopSequences(clusters);
	std::cerr << "write self-corrected ONT loop sequences" << std::endl;
	writeSelfCorrectedOntLoopSequences(params.basePath + "/loops_selfcorrected.fa", clusters);
	std::cerr << "getting morph consensuses" << std::endl;
	auto morphConsensuses = getMorphConsensuses(clusters, graph, pathStartClip, pathEndClip, params.namePrefix);
	std::cerr << "polishing morph consensuses" << std::endl;
	polishMorphConsensuses(morphConsensuses, clusters, params.MBGPath, params.basePath + "/tmp", params.numThreads);
	std::cerr << "write morph consensuses" << std::endl;
	writeMorphConsensuses(params.basePath + "/morphs.fa", params.basePath + "/morphs_preconsensus.fa", morphConsensuses);
	std::cerr << "write morph paths" << std::endl;
	writeMorphPaths(params.basePath + "/morphs.gaf", morphConsensuses, graph, pathStartClip, pathEndClip);
	std::cerr << "write morph graph and read paths" << std::endl;
	writeMorphGraphAndReadPaths(params.basePath + "/morphgraph.gfa", params.basePath + "/readpaths-morphgraph.gaf", morphConsensuses);
	if (params.winnowmapPath != "" && params.samtoolsPath != "")
	{
		std::cerr << "realign raw ONT loop sequences to morph consensuses" << std::endl;
		alignRawOntLoopsToMorphConsensuses(clusters, params.basePath + "/raw_loop_to_morphs_alignments.bam", params.basePath + "/tmp", params.numThreads, morphConsensuses, params.winnowmapPath, params.samtoolsPath);
		std::cerr << "realign self-corrected ONT loop sequences to morph consensuses" << std::endl;
		alignSelfCorrectedOntLoopsToMorphConsensuses(clusters, params.basePath + "/selfcorrected_loop_to_morphs_alignments.bam", params.basePath + "/tmp", params.numThreads, morphConsensuses, params.winnowmapPath, params.samtoolsPath);
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
	if (params.annotationFasta.size() > 0)
	{
		std::cerr << "lifting over annotations to morphs" << std::endl;
		liftoverAnnotationsToMorphs(params.basePath, morphConsensuses, params.annotationFasta, params.annotationGff3, params.basePath + "/tmp", params.liftoffPath);
	}
}
