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
#include "fastqloader.h"
#include "ClusterHandler.h"
#include "RibotinUtils.h"

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
				edges[fromnode].emplace(tonode, overlap, coverage);
				edges[reverse(tonode)].emplace(reverse(fromnode), overlap, coverage);
			}
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

Path getHeavyPath(const GfaGraph& graph)
{
	size_t maxCoverageNode = 0;
	for (size_t i = 1; i < graph.nodeCoverages.size(); i++)
	{
		if (graph.nodeCoverages[i] > graph.nodeCoverages[maxCoverageNode]) maxCoverageNode = i;
	}
	std::unordered_map<Node, std::pair<Node, double>> predecessor;
	std::priority_queue<std::tuple<Node, Node, double>, std::vector<std::tuple<Node, Node, double>>, CoverageComparer> checkStack;
	checkStack.emplace(Node { maxCoverageNode, true }, Node { maxCoverageNode, true }, 0);
	while (checkStack.size() > 0)
	{
		auto top = checkStack.top();
		checkStack.pop();
		if (graph.edges.count(std::get<0>(top)) == 0) continue;
		if (predecessor.count(std::get<0>(top)) == 1 && predecessor.at(std::get<0>(top)).second >= std::get<2>(top)) continue;
		predecessor[std::get<0>(top)] = std::make_pair(std::get<1>(top), std::get<2>(top));
		for (auto edge : graph.edges.at(std::get<0>(top)))
		{
			double coverage = std::min(graph.nodeCoverages.at(std::get<0>(edge).id()), std::get<2>(edge));
			if (predecessor.count(std::get<0>(edge)) == 1 && predecessor.at(std::get<0>(edge)).second > coverage) continue;
			checkStack.emplace(std::get<0>(edge), std::get<0>(top), coverage);
		}
	}
	Path result;
	result.nodes.emplace_back(maxCoverageNode, true);
	std::unordered_set<Node> visited;
	while (true)
	{
		assert(predecessor.count(result.nodes.back()) == 1);
		auto pre = predecessor.at(result.nodes.back());
		if (std::get<0>(pre).id() == maxCoverageNode)
		{
			break;
		}
		if (visited.count(reverse(pre.first)) == 1)
		{
			std::cerr << "The consensus sequence has a palindrome. Unable to build consensus." << std::endl;
			std::abort();
		}
		assert(visited.count(pre.first) == 0);
		visited.insert(pre.first);
		result.nodes.emplace_back(pre.first);
	}
	std::reverse(result.nodes.begin(), result.nodes.end());
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

void writePathSequence(const Path& path, const GfaGraph& graph, std::string outputFile)
{
	std::string pathseq = path.getSequence(graph.nodeSeqs);
	std::ofstream file { outputFile };
	file << ">heavy_path" << std::endl;
	file << pathseq;
}

void writePathGaf(const Path& path, const GfaGraph& graph, std::string outputFile)
{
	std::ofstream file { outputFile };
	file << getPathGaf(path, graph);
}

void runMBG(std::string basePath, std::string readPath, std::string MBGPath, size_t k)
{
	std::string mbgCommand;
	mbgCommand = MBGPath + " -o " + basePath + "/graph.gfa -i " + readPath + " -k " + std::to_string(k) + " -w " + std::to_string(k-30) + " -a 2 -u 3 -r 15000 -R 4000 --error-masking=msat --output-sequence-paths " + basePath + "/paths.gaf --only-local-resolve 1> " + basePath + "/mbg_stdout.txt 2> " + basePath + "/mbg_stderr.txt";
	std::cerr << "MBG command:" << std::endl;
	std::cerr << mbgCommand << std::endl;
	int result = system(mbgCommand.c_str());
	if (result != 0)
	{
		std::cerr << "MBG did not run successfully" << std::endl;
		std::abort();
	}
}

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
		size_t readstart = std::stoull(parts[2]);
		size_t readend = std::stoull(parts[3]);
		std::string pathstr = parts[5];
		size_t pathLength = std::stoull(parts[6]);
		size_t pathStart = std::stoull(parts[6]);
		size_t pathEnd = std::stoull(parts[6]);
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
	assert(startIndex != endIndex);
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
		std::string referenceSeq = getSequence(variants[variant].referencePath, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
		std::string variantSeq = getSequence(variants[variant].path, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
		assert(variants[variant].referenceEndPos < variants[variant].referenceStartPos || referenceSeq.size() == variants[variant].referenceEndPos - variants[variant].referenceStartPos);
		assert(variants[variant].referenceEndPos > variants[variant].referenceStartPos || referenceSeq.size() == pathLength + variants[variant].referenceEndPos - variants[variant].referenceStartPos);
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
	std::sort(variants.begin(), variants.end(), [](const Variant& left, const Variant& right) { return left.referenceStartPos < right.referenceStartPos; });
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

void liftoverAnnotationsToConsensus(const std::string& basePath, const std::string& consensusPath, const std::string& annotationFasta, const std::string& annotationGff3)
{
	{
		std::ofstream typefile { basePath + "/liftoff_types.txt" };
		typefile << "rRNA" << std::endl;
		typefile << "misc_RNA" << std::endl;
		typefile << "repeat_region" << std::endl;
		typefile << "gene" << std::endl;
		typefile << "pseudogene" << std::endl;
	}
	std::string command = "liftoff -f " + basePath + "/liftoff_types.txt -g " + annotationGff3 + " -o " + basePath + "/annotation.gff3 -u " + basePath + "/unmapped_features.txt -dir " + basePath + "/liftoff_intermediate_files/ " + consensusPath + " " + annotationFasta;
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
			if (getPathLength(looppath, graph.nodeSeqs, graph.edges) >= minLength)
			{
				result.emplace_back();
				result.back().originalReadLength = read.readLength;
				result.back().path = looppath;
				result.back().readName = read.readName;
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

size_t getEditDistanceWfa(const std::string& left, const std::string& right, const size_t maxEdits)
{
	if (left.size() > right.size()+maxEdits) return maxEdits+1;
	if (right.size() > left.size()+maxEdits) return maxEdits+1;
	std::vector<size_t> offsets;
	std::vector<size_t> nextOffsets;
	offsets.resize(1, 0);
	while (offsets[0] < left.size() && offsets[0] < right.size() && left[offsets[0]] == right[offsets[0]]) offsets[0] += 1;
	if (offsets[0] == left.size()) return right.size() - left.size();
	if (offsets[0] == right.size()) return left.size() - right.size();
	size_t zeroPos = 1;
	for (size_t score = 0; score < maxEdits; score++)
	{
		nextOffsets.resize(offsets.size()+2);
		for (size_t i = 0; i < nextOffsets.size(); i++)
		{
			nextOffsets[i] = 0;
			if (i >= 1 && i < nextOffsets.size()-1) nextOffsets[i] = offsets[i-1]+1;
			if (i >= 2) nextOffsets[i] = std::max(nextOffsets[i], offsets[i-2]);
			if (i < nextOffsets.size()-2) nextOffsets[i] = std::max(nextOffsets[i], offsets[i]+1);
			assert(nextOffsets[i] + i >= zeroPos);
			nextOffsets[i] = std::min(nextOffsets[i], right.size() - i + zeroPos);
			nextOffsets[i] = std::min(nextOffsets[i], left.size());
			while (nextOffsets[i] < left.size() && nextOffsets[i]+i-zeroPos < right.size() && left[nextOffsets[i]] == right[nextOffsets[i]+i-zeroPos]) nextOffsets[i] += 1;
			if (nextOffsets[i] == left.size() && nextOffsets[i]+i-zeroPos == right.size()) return score;
		}
		std::swap(offsets, nextOffsets);
		zeroPos += 1;
	}
	return maxEdits+1;
}

size_t getEditDistancePossiblyMemoized(const std::vector<Node>& left, const std::vector<Node>& right, const size_t leftStartClipBp, const size_t rightStartClipBp, const size_t leftEndClipBp, const size_t rightEndClipBp, const GfaGraph& graph, const size_t maxEdits, std::unordered_map<std::pair<std::vector<Node>, std::vector<Node>>, size_t>& memoizedEditDistances)
{
	if (left == right) return 0;
	auto key = canon(left, right);
	size_t add;
	if (memoizedEditDistances.count(key) == 0)
	{
		size_t startClip = 0;
		size_t endClip = 0;
		if (left[0] == right[0])
		{
			while (startClip+1 < left.size() && startClip+1 < right.size() && left[startClip+1] == right[startClip+1]) startClip += 1;
		}
		if (left.back() == right.back())
		{
			while (startClip+endClip+1 < left.size() && startClip+endClip+1 < right.size() && left[left.size()-2-endClip] == right[right.size()-2-endClip]) endClip += 1;
		}
		if (startClip > 0 || endClip > 0)
		{
			assert(startClip + endClip < left.size());
			assert(startClip + endClip < right.size());
			auto leftSubseq = getSequence(std::vector<Node> { left.begin() + startClip, left.end() - endClip }, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
			auto rightSubseq = getSequence(std::vector<Node> { right.begin() + startClip, right.end() - endClip }, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
			if (leftStartClipBp != 0 || leftEndClipBp != 0) leftSubseq = leftSubseq.substr(leftStartClipBp, left.size() - leftStartClipBp - leftEndClipBp);
			if (rightStartClipBp != 0 || rightEndClipBp != 0) rightSubseq = rightSubseq.substr(rightStartClipBp, right.size() - rightStartClipBp - rightEndClipBp);
			add = getEditDistanceWfa(leftSubseq, rightSubseq, maxEdits);
			memoizedEditDistances[key] = add;
		}
		else
		{
			auto leftSubseq = getSequence(left, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
			auto rightSubseq = getSequence(right, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
			if (leftStartClipBp != 0 || leftEndClipBp != 0) leftSubseq = leftSubseq.substr(leftStartClipBp, left.size() - leftStartClipBp - leftEndClipBp);
			if (rightStartClipBp != 0 || rightEndClipBp != 0) rightSubseq = rightSubseq.substr(rightStartClipBp, right.size() - rightStartClipBp - rightEndClipBp);
			add = getEditDistanceWfa(leftSubseq, rightSubseq, maxEdits);
			memoizedEditDistances[key] = add;
		}
	}
	else
	{
		add = memoizedEditDistances.at(key);
	}
	return add;
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
	std::vector<size_t> maxIncreasing;
	maxIncreasing.resize(unfilteredMatches.size(), 1);
	for (size_t i = 1; i < unfilteredMatches.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (unfilteredMatches[i].first <= unfilteredMatches[j].first) continue;
			if (unfilteredMatches[i].second <= unfilteredMatches[j].second) continue;
			maxIncreasing[i] = std::max(maxIncreasing[i], maxIncreasing[j]+1);
		}
	}
	if (maxIncreasing.back() == unfilteredMatches.size()) return unfilteredMatches;
	size_t maxpos = 0;
	for (size_t i = 1; i < maxIncreasing.size(); i++)
	{
		if (maxIncreasing[i] > maxIncreasing[maxpos]) maxpos = i;
	}
	std::vector<std::pair<size_t, size_t>> result;
	while (true)
	{
		result.push_back(unfilteredMatches[maxpos]);
		if (maxIncreasing[maxpos] == 1) break;
		size_t next = maxpos;
		for (size_t j = maxpos-1; j < maxpos; j--)
		{
			if (maxIncreasing[j] == maxIncreasing[maxpos]-1)
			{
				next = j;
				break;
			}
		}
		assert(next < maxpos);
		maxpos = next;
	}
	std::reverse(result.begin(), result.end());
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
		add = getEditDistancePossiblyMemoized(leftPath, rightPath, 0, 0, 0, 0, graph, maxEdits-result, memoizedEditDistances);
		result += add;
		if (result >= maxEdits) return maxEdits+1;
	}
	std::vector<Node> leftPath { left.begin(), left.begin()+nodeMatches[0].first+1 };
	std::vector<Node> rightPath { right.begin(), right.begin()+nodeMatches[0].second+1 };
	add = getEditDistancePossiblyMemoized(leftPath, rightPath, pathStartClip.at(leftPath[0]), pathStartClip.at(rightPath[0]), 0, 0, graph, maxEdits-result, memoizedEditDistances);
	result += add;
	if (result >= maxEdits) return maxEdits+1;
	leftPath = std::vector<Node> { left.begin()+nodeMatches.back().first, left.end() };
	rightPath = std::vector<Node> { right.begin()+nodeMatches.back().second, right.end() };
	add = getEditDistancePossiblyMemoized(leftPath, rightPath, 0, 0, pathEndClip.at(leftPath.back()), pathEndClip.at(rightPath.back()), graph, maxEdits-result, memoizedEditDistances);
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

std::vector<std::vector<OntLoop>> clusterLoopSequences(const std::vector<OntLoop>& loops, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::unordered_set<size_t>& coreNodes, const size_t maxEdits)
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
		loopLengths.emplace_back(getPathLength(loops[i].path, graph.nodeSeqs, graph.edges));
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
			if (edits > maxEdits) continue;
			parent[j] = parent[i];
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
	for (const auto& read : rawPaths)
	{
		const auto& path = read.path;
		size_t lastCore = 0;
		for (size_t i = 0; i < path.size(); i++)
		{
			if (coreNodes.count(path[i].id()) == 0) continue;
			std::vector<Node> subpath { path.begin() + lastCore, path.begin() + i + 1 };
			assert(subpath.size() >= 2 || (i == 0 && subpath.size() == 1));
			alleleCounts[subpath] += 1;
			lastCore = i;
		}
		// assert(coreIndex == alleleCounts.size()-1);
		std::vector<Node> subpath { path.begin() + lastCore, path.end()};
		assert(subpath.size() >= 2 || (lastCore == path.size()-1 && subpath.size() == 1));
		alleleCounts[subpath] += 1;
	}
	auto corePath = getMajorityPath(coreNodes, rawPaths);
	std::unordered_map<size_t, size_t> coreNodePositionInPath;
	for (size_t i = 0; i < corePath.size(); i++)
	{
		coreNodePositionInPath[corePath[i]] = i;
	}
	std::vector<std::pair<std::vector<Node>, size_t>> bestAlleles;
	bestAlleles.resize(corePath.size()+1, std::make_pair(std::vector<Node>{}, 0));
	for (auto pair : alleleCounts)
	{
		size_t index = std::numeric_limits<size_t>::max();
		if ((coreNodePositionInPath.count(pair.first[0].id()) == 0 || pair.first.size() == 1) && coreNodePositionInPath.count(pair.first.back().id()) == 1 && coreNodePositionInPath.at(pair.first.back().id()) == 0)
		{
			index = 0;
		}
		else if (pair.first.size() >= 2 && coreNodePositionInPath.count(pair.first[0].id()) == 1 && coreNodePositionInPath.count(pair.first.back().id()) == 1 && coreNodePositionInPath.at(pair.first.back().id()) == coreNodePositionInPath.at(pair.first[0].id())+1)
		{
			index = coreNodePositionInPath.at(pair.first.back().id());
		}
		else if (coreNodePositionInPath.count(pair.first[0].id()) == 1 && (coreNodePositionInPath.count(pair.first.back().id()) == 0 || pair.first.size() == 1) && coreNodePositionInPath.at(pair.first[0].id()) == corePath.size()-1)
		{
			index = corePath.size();
		}
		if (index == std::numeric_limits<size_t>::max()) continue;
		if (pair.second > bestAlleles[index].second)
		{
			bestAlleles[index] = pair;
		}
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
	return consensusPath;
}

std::string getConsensusSequence(const std::vector<Node>& consensusPath, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
{
	std::string consensusSeq = getSequence(consensusPath, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
	consensusSeq = consensusSeq.substr(pathStartClip.at(consensusPath[0]), consensusSeq.size() - pathStartClip.at(consensusPath[0]) - pathEndClip.at(consensusPath.back()));
	return consensusSeq;
}

std::vector<MorphConsensus> getMorphConsensuses(const std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
{
	std::vector<MorphConsensus> result;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		result.emplace_back();
		// todo handle cases where plurality path is not consensus
		result.back().path = getConsensusPath(clusters[i], graph);
		result.back().sequence = getConsensusSequence(result.back().path, graph, pathStartClip, pathEndClip);
		result.back().coverage = clusters[i].size();
		result.back().ontLoops = clusters[i];
		result.back().name = "morphconsensus" + std::to_string(i) + "_coverage" + std::to_string(result.back().coverage);
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

void orderLoopsByLength(std::vector<OntLoop>& loops, const GfaGraph& graph)
{
	std::vector<std::pair<size_t, size_t>> indexAndLength;
	indexAndLength.reserve(loops.size());
	for (size_t i = 0; i < loops.size(); i++)
	{
		indexAndLength.emplace_back(i, getPathLength(loops[i].path, graph.nodeSeqs, graph.edges));
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
	// runMBG(params.basePath, params.hifiReadPath, params.MBGPath, params.k);
	std::cerr << "reading graph" << std::endl;
	GfaGraph graph;
	graph.loadFromFile(params.basePath + "/graph.gfa");
	std::cerr << "getting consensus" << std::endl;
	Path heavyPath = getHeavyPath(graph);
	if (params.orientReferencePath.size() > 0)
	{
		std::cerr << "orienting consensus" << std::endl;
		heavyPath = orientPath(graph, heavyPath, params.orientReferencePath, 101);
	}
	std::cerr << "writing consensus" << std::endl;
	writePathSequence(heavyPath, graph, params.basePath + "/consensus.fa");
	writePathGaf(heavyPath, graph, params.basePath + "/consensus_path.gaf");
	std::cerr << "reading read paths" << std::endl;
	std::vector<ReadPath> readPaths = loadReadPaths(params.basePath + "/paths.gaf", graph);
	std::cerr << "getting variants" << std::endl;
	std::vector<Variant> variants = getVariants(graph, heavyPath, readPaths, 3);
	nameVariants(variants, graph, heavyPath);
	std::cerr << "writing variants" << std::endl;
	writeVariants(heavyPath, graph, variants, params.basePath + "/variants.txt");
	std::cerr << "writing variant graph" << std::endl;
	writeVariantGraph(params.basePath + "/graph.gfa", graph, heavyPath, variants, params.basePath + "/variant-graph.gfa");
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
		size_t lastMatch = 0;
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
		size_t lastMatch = 0;
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

void DoClusterONTAnalysis(const ClusterParams& params)
{
	std::cerr << "reading graph" << std::endl;
	GfaGraph graph;
	graph.loadFromFile(params.basePath + "/graph.gfa");
	std::cerr << "getting consensus" << std::endl;
	Path heavyPath = getHeavyPath(graph);
	if (params.orientReferencePath.size() > 0)
	{
		std::cerr << "orienting consensus" << std::endl;
		heavyPath = orientPath(graph, heavyPath, params.orientReferencePath, 101);
	}
	std::cerr << "extract corrected ultralong paths" << std::endl;
	size_t heavyPathLength = heavyPath.getSequence(graph.nodeSeqs).size();
	size_t minLength = heavyPathLength * 0.5;
	std::cerr << "consensus path length " << heavyPathLength << ", using " << minLength << " as minimum morph length" << std::endl;
	auto ontPaths = extractCorrectedONTPaths(params.basePath + "/ont-alns.gaf", heavyPath, minLength, graph);
	std::cerr << ontPaths.size() << " corrected paths" << std::endl;
	std::cerr << "extract morph paths from ONTs" << std::endl;
	std::unordered_map<Node, size_t> pathStartClip;
	std::unordered_map<Node, size_t> pathEndClip;
	std::unordered_set<size_t> borderNodes;
	std::tie(borderNodes, pathStartClip, pathEndClip) = getBorderNodes(heavyPath, graph);
	assert(borderNodes.size() > 0);
	auto loopSequences = extractLoopSequences(ontPaths, heavyPath, minLength, graph, borderNodes, pathStartClip, pathEndClip);
	auto coreNodes = getCoreNodes(loopSequences);
	std::cerr << loopSequences.size() << " morph paths in ONTs" << std::endl;
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
	std::cerr << "cluster morphs" << std::endl;
	orderLoopsByLength(loopSequences, graph);
	auto unphasedClusters = clusterLoopSequences(loopSequences, graph, pathStartClip, pathEndClip, coreNodes, 300);
	std::cerr << unphasedClusters.size() << " morph clusters" << std::endl;
	std::cerr << "phase morph clusters" << std::endl;
	auto clusters = phaseClusters(unphasedClusters);
	std::cerr << clusters.size() << " phased morph clusters" << std::endl;
	std::sort(clusters.begin(), clusters.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
	std::cerr << "getting morph consensuses" << std::endl;
	auto morphConsensuses = getMorphConsensuses(clusters, graph, pathStartClip, pathEndClip);
	std::cerr << "write morph consensuses" << std::endl;
	writeMorphConsensuses(params.basePath + "/morphs.fa", morphConsensuses);
	std::cerr << "write morph paths" << std::endl;
	writeMorphPaths(params.basePath + "/morphs.gaf", morphConsensuses, graph, pathStartClip, pathEndClip);
	std::cerr << "write morph graph and read paths" << std::endl;
	writeMorphGraphAndReadPaths(params.basePath + "/morphgraph.gfa", params.basePath + "/readpaths-morphgraph.gaf", morphConsensuses);
}
