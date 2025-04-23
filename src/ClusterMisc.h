#ifndef ClusterMisc_h
#define ClusterMisc_h

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <string>
#include <cstddef>
#include <vector>
#include "phmap.h"
#include "TwobitString.h"

class Node
{
public:
	Node() = default;
	Node(const Node& other) = default;
	Node(Node&& other) = default;
	Node& operator=(const Node& other) = default;
	Node& operator=(Node&& other) = default;
	Node(size_t id, bool forward);
	size_t id() const;
	bool forward() const;
	bool operator<(const Node& other) const;
	bool operator==(const Node& other) const;
	bool operator!=(const Node& other) const;
	size_t rawValue() const;
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
	enum MorphType
	{
		Inner,
		BorderOne,
		BorderTwo,
		Isolated
	};
	std::vector<Node> path;
	std::string sequence;
	std::string preCorrectionSequence;
	std::string name;
	std::vector<OntLoop> ontLoops;
	size_t coverage;
	MorphType type;
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
	void loadFromFile(std::string gfaFile);
	size_t numNodes() const;
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
	Path();
	std::vector<Node> nodes;
	std::vector<size_t> overlaps;
	size_t leftClip;
	size_t rightClip;
	std::string getSequence(const std::vector<std::string>& nodeSeqs) const;
};

class PathSequenceView
{
public:
	PathSequenceView(const std::vector<Node>& nodes, const std::vector<TwobitString>& nodeSeqs, const std::vector<TwobitString>& revCompNodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges, size_t leftClip, size_t rightClip);
	char operator[](size_t index) const;
	size_t size() const;
private:
	std::vector<const TwobitString*> nodePointers;
	std::vector<size_t> nodeStartPoses;
	std::vector<size_t> nodeStartOverlap;
	size_t leftClip;
	size_t rightClip;
	friend size_t getMatchLength(const PathSequenceView& left, const PathSequenceView& right, const size_t leftStart, const size_t rightStart);
};

std::string revcomp(std::string seq);
Node reverse(Node node);
std::string reverse(std::string node);
std::vector<Node> reverse(const std::vector<Node>& path);
std::vector<std::string> reverse(const std::vector<std::string>& path);
std::pair<std::vector<Node>, std::vector<Node>> canon(const std::vector<Node>& left, const std::vector<Node>& right);
std::pair<std::vector<std::string>, std::vector<std::string>> canon(const std::vector<std::string>& left, const std::vector<std::string>& right);
std::vector<Node> canon(const std::vector<Node>& path);
std::vector<std::string> canon(const std::vector<std::string>& path);
std::pair<Node, Node> canon(Node left, Node right);
std::pair<std::string, std::string> canon(const std::string& left, const std::string& right);
size_t find(std::vector<size_t>& parent, size_t key);
void merge(std::vector<size_t>& parent, size_t left, size_t right);
size_t getPathLength(const std::vector<Node>& nodes, const std::vector<std::string>& nodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges);
size_t getOverlap(const Node& from, const Node& to, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges);
PathSequenceView getSequenceView(const std::vector<Node>& nodes, const std::vector<TwobitString>& nodeSeqs, const std::vector<TwobitString>& revCompNodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges, size_t leftClip, size_t rightClip);
Path orientPath(const GfaGraph& graph, const Path& rawPath, const std::string& orientReferencePath, size_t k);
std::unordered_set<size_t> getCoreNodes(const std::vector<OntLoop>& paths);
bool orientPath(std::vector<Node>& path, const std::unordered_map<size_t, bool>& referenceOrientations);

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

#endif
