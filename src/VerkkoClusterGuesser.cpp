#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include "fastqloader.h"
#include "KmerMatcher.h"
#include "VerkkoClusterGuesser.h"

std::string revnode(std::string node)
{
	assert(node.size() >= 2);
	assert(node[0] == '>' || node[0] == '<');
	return (node[0] == '>' ? "<" : ">") + node.substr(1);
}

std::unordered_set<std::string> matchNodes(const KmerMatcher& matcher, const std::string& gfaPath)
{
	std::unordered_set<std::string> result;
	std::ifstream file { gfaPath };
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (!file.good()) break;
		if (line[0] != 'S') continue;
		std::string dummy, nodename, nodeseq;
		std::stringstream sstr { line };
		sstr >> dummy >> nodename >> nodeseq;
		if (nodeseq.size() >= 100000) continue;
		size_t matchLength = matcher.getMatchLength(nodeseq);
		if (matchLength < 2000) continue;
		result.insert(nodename);
	}
	return result;
}

std::string find(std::unordered_map<std::string, std::string>& parent, std::string key)
{
	while (parent.at(key) != parent.at(parent.at(key))) parent[key] = parent[parent[key]];
	return parent[key];
}

void merge(std::unordered_map<std::string, std::string>& parent, std::string left, std::string right)
{
	left = find(parent, left);
	right = find(parent, right);
	assert(parent.at(left) == left);
	assert(parent.at(right) == right);
	parent[right] = left;
}

std::vector<std::vector<std::string>> extendClusters(const std::unordered_set<std::string>& validNodes, const std::string& gfaPath)
{
	std::unordered_map<std::string, std::string> parent;
	for (auto node : validNodes)
	{
		parent[node] = node;
	}
	{
		std::ifstream file { gfaPath };
		while (file.good())
		{
			std::string line;
			getline(file, line);
			if (!file.good()) break;
			if (line[0] != 'L') continue;
			std::stringstream sstr { line };
			std::string dummy, fromnode, tonode;
			sstr >> dummy >> fromnode >> dummy >> tonode;
			if (validNodes.count(fromnode) == 0) continue;
			if (validNodes.count(tonode) == 0) continue;
			merge(parent, fromnode, tonode);
		}
	}
	for (auto node : validNodes)
	{
		find(parent, node);
	}
	std::unordered_map<std::string, size_t> clusterSize;
	for (auto node : validNodes)
	{
		clusterSize[find(parent, node)] += 1;
	}
	std::unordered_map<std::string, size_t> clusterNumber;
	size_t numClusters = 0;
	for (auto pair : clusterSize)
	{
		if (pair.second < 10) continue;
		clusterNumber[pair.first] = numClusters;
		numClusters += 1;
	}
	std::vector<std::vector<std::string>> result;
	result.resize(numClusters);
	for (auto node : validNodes)
	{
		auto cluster = find(parent, node);
		if (clusterNumber.count(cluster) == 0) continue;
		result[clusterNumber.at(cluster)].push_back(node);
	}
	return result;
}

std::string homopolymerCompress(std::string original)
{
	std::string result;
	result.push_back(original[0]);
	for (size_t i = 1; i < original.size(); i++)
	{
		if (original[i] != original[i-1]) result.push_back(original[i]);
	}
	return result;
}

bool hasCycle(const std::vector<std::string>& clusterNodes, const std::unordered_map<std::string, std::vector<std::string>>& edges)
{
	std::unordered_map<std::string, std::unordered_set<std::string>> reachable;
	for (auto node : clusterNodes)
	{
		if (edges.count(">" + node) == 1)
		{
			for (auto edge : edges.at(">" + node))
			{
				reachable[">" + node].insert(edge);
				if (edge == ">" + node) return true;
			}
		}
		if (edges.count("<" + node) == 1)
		{
			for (auto edge : edges.at("<" + node))
			{
				reachable["<" + node].insert(edge);
				if (edge == "<" + node) return true;
			}
		}
	}
	while (true)
	{
		bool changed = false;
		for (auto node : clusterNodes)
		{
			auto oldReachable = reachable[">" + node];
			for (auto node2 : oldReachable)
			{
				if (edges.count(node2) == 0) continue;
				for (auto edge : edges.at(node2))
				{
					if (reachable[">" + node].count(edge) == 1) continue;
					reachable[">" + node].insert(edge);
					changed = true;
					if (edge == ">" + node) return true;
				}
			}
			oldReachable = reachable["<" + node];
			for (auto node2 : oldReachable)
			{
				if (edges.count(node2) == 0) continue;
				for (auto edge : edges.at(node2))
				{
					if (reachable["<" + node].count(edge) == 1) continue;
					reachable["<" + node].insert(edge);
					changed = true;
					if (edge == "<" + node) return true;
				}
			}
		}
		if (!changed) break;
	}
	return false;
}

void filterOutAcyclic(std::vector<std::vector<std::string>>& clusterNodes, const std::string& gfaPath)
{
	std::unordered_map<std::string, std::vector<std::string>> edges;
	{
		std::ifstream file { gfaPath };
		while (file.good())
		{
			std::string line;
			getline(file, line);
			if (!file.good()) break;
			if (line[0] != 'L') continue;
			std::stringstream sstr { line };
			std::string dummy, fromnode, tonode;
			std::string fromorient, toorient;
			sstr >> dummy >> fromnode >> fromorient >> tonode >> toorient;
			assert(fromorient == "+" || fromorient == "-");
			assert(toorient == "+" || toorient == "-");
			fromnode = (fromorient == "+" ? ">" : "<") + fromnode;
			tonode = (toorient == "+" ? ">" : "<") + tonode;
			edges[fromnode].push_back(tonode);
			edges[revnode(tonode)].push_back(revnode(fromnode));
		}
	}
	for (size_t i = clusterNodes.size()-1; i < clusterNodes.size(); i--)
	{
		if (hasCycle(clusterNodes[i], edges)) continue;
		std::swap(clusterNodes[i], clusterNodes.back());
		clusterNodes.pop_back();
	}
}

std::vector<std::vector<std::string>> guessVerkkoRDNAClusters(std::string verkkoBasePath, const std::vector<std::string>& referencePath)
{
	KmerMatcher matcher { 101 };
	for (auto file : referencePath)
	{
		FastQ::streamFastqFromFile(file, false, [&matcher](FastQ& fastq)
		{
			matcher.addReferenceKmers(homopolymerCompress(fastq.sequence));
		});
	}
	auto nodes = matchNodes(matcher, verkkoBasePath + "/assembly.homopolymer-compressed.gfa");
	auto clusters = extendClusters(nodes, verkkoBasePath + "/assembly.homopolymer-compressed.gfa");
	filterOutAcyclic(clusters, verkkoBasePath + "/assembly.homopolymer-compressed.gfa");
	return clusters;
}
