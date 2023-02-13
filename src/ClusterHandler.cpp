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
#include "ClusterHandler.h"

std::string revcomp(std::string seq)
{
	std::reverse(seq.begin(), seq.end());
	for (size_t i = 0; i < seq.size(); i++)
	{
		switch(seq[i])
		{
			case 'a':
			case 'A':
				seq[i] = 'T';
				break;
			case 'c':
			case 'C':
				seq[i] = 'G';
				break;
			case 'g':
			case 'G':
				seq[i] = 'C';
				break;
			case 't':
			case 'T':
				seq[i] = 'A';
				break;
			default:
				assert(false);
				break;
		}
	}
	return seq;
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

std::vector<std::string> reverse(const std::vector<std::string>& path)
{
	std::vector<std::string> bw { path.rbegin(), path.rend() };
	for (size_t i = 0; i < bw.size(); i++)
	{
		bw[i] = reverse(bw[i]);
	}
	return bw;
}

std::vector<std::string> canon(const std::vector<std::string>& path)
{
	auto bw = reverse(path);
	if (bw < path) return bw;
	return path;
}

class ReadPath
{
public:
	std::string readName;
	size_t readStart;
	size_t readEnd;
	std::vector<std::string> path;
	size_t pathStartClip;
	size_t pathEndClip;
};

class Variant
{
public:
	std::vector<std::string> path;
	std::vector<std::string> referencePath;
	size_t coverage;
	size_t referenceCoverage;
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
				nodeSeqs[node] = seq;
				nodeCoverages[node] = std::stof(coveragetag.substr(5));
			}
			else if (test == "L")
			{
				std::string fromnode;
				std::string tonode;
				std::string fromorient, toorient;
				std::string overlapstr, edgecoveragetag;
				sstr >> fromnode >> fromorient >> tonode >> toorient >> overlapstr >> edgecoveragetag;
				if (fromorient == "+")
				{
					fromnode = ">" + fromnode;
				}
				else
				{
					fromnode = "<" + fromnode;
				}
				if (toorient == "+")
				{
					tonode = ">" + tonode;
				}
				else
				{
					tonode = "<" + tonode;
				}
				overlapstr.pop_back();
				size_t overlap = std::stoull(overlapstr);
				size_t coverage = std::stoull(edgecoveragetag.substr(5));
				edges[fromnode].emplace(tonode, overlap, coverage);
				edges[reverse(tonode)].emplace(reverse(fromnode), overlap, coverage);
			}
		}
	}
	std::unordered_map<std::string, size_t> nodeCoverages;
	std::unordered_map<std::string, std::string> nodeSeqs;
	std::unordered_map<std::string, std::set<std::tuple<std::string, size_t, size_t>>> edges;
};

class Path
{
public:
	std::vector<std::string> nodes;
	std::vector<size_t> overlaps;
	std::string getSequence(const std::unordered_map<std::string, std::string>& nodeSeqs) const
	{
		std::string result;
		for (size_t i = 0; i < nodes.size(); i++)
		{
			std::string add;
			add = nodeSeqs.at(nodes[i].substr(1));
			if (nodes[i][0] == '<')
			{
				add = revcomp(add);
			}
			else
			{
				assert(nodes[i][0] == '>');
			}
			add = add.substr(overlaps[i]);
			result += add;
		}
		return result;
	}
};

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
		result += node;
	}
	result += "\t";
	result += std::to_string(pathseq.size());
	result += std::to_string(path.overlaps[0]);
	result += "\t";
	result += std::to_string(path.overlaps[0]);
	result += "\t";
	result += std::to_string(pathseq.size());
	result += std::to_string(path.overlaps[0]);
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
	bool operator()(const std::tuple<std::string, std::string, double>& lhs, const std::tuple<std::string, std::string, double>& rhs) const
	{
		return std::get<2>(lhs) < std::get<2>(rhs);
	}
};

Path getHeavyPath(const GfaGraph& graph)
{
	std::string maxCoverageNode;
	for (auto pair : graph.nodeCoverages)
	{
		if (maxCoverageNode == "" || pair.second > graph.nodeCoverages.at(maxCoverageNode))
		{
			maxCoverageNode = pair.first;
		}
	}
	std::unordered_map<std::string, std::pair<std::string, double>> predecessor;
	std::priority_queue<std::tuple<std::string, std::string, double>, std::vector<std::tuple<std::string, std::string, double>>, CoverageComparer> checkStack;
	checkStack.emplace(">" + maxCoverageNode, ">" + maxCoverageNode, 0);
	while (checkStack.size() > 0)
	{
		auto top = checkStack.top();
		checkStack.pop();
		if (graph.edges.count(std::get<0>(top)) == 0) continue;
		if (predecessor.count(std::get<0>(top)) == 1 && predecessor.at(std::get<0>(top)).second >= std::get<2>(top)) continue;
		predecessor[std::get<0>(top)] = std::make_pair(std::get<1>(top), std::get<2>(top));
		for (auto edge : graph.edges.at(std::get<0>(top)))
		{
			double coverage = std::min(graph.nodeCoverages.at(std::get<0>(edge).substr(1)), std::get<2>(edge));
			if (predecessor.count(std::get<0>(edge)) == 1 && predecessor.at(std::get<0>(edge)).second > coverage) continue;
			checkStack.emplace(std::get<0>(edge), std::get<0>(top), coverage);
		}
	}
	Path result;
	result.nodes.emplace_back(">" + maxCoverageNode);
	std::unordered_set<std::string> visited;
	while (true)
	{
		assert(predecessor.count(result.nodes.back()) == 1);
		auto pre = predecessor.at(result.nodes.back());
		if (std::get<0>(pre).substr(1) == maxCoverageNode)
		{
			break;
		}
		assert(visited.count(pre.first) == 0);
		result.nodes.emplace_back(pre.first);
	}
	std::reverse(result.nodes.begin(), result.nodes.end());
	for (size_t i = 0; i < result.nodes.size(); i++)
	{
		std::string prev;
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

void runMBG(std::string basePath, std::string readPath, std::string MBGPath)
{
	std::string mbgCommand;
	mbgCommand = "/usr/bin/time -v " + MBGPath + " -o " + basePath + "/graph.gfa -i " +readPath + " -k 51 -w 20 -a 2 -u 3 -r 15000 -R 4000 --error-masking=msat --output-sequence-paths " + basePath + "/paths.gaf --only-local-resolve 1> " + basePath + "/mbg_stdout.txt 2> " + basePath + "/mbg_stderr.txt";
	std::cerr << "MBG command:" << std::endl;
	std::cerr << mbgCommand << std::endl;
	system(mbgCommand.c_str());
}

std::vector<std::string> parsePath(const std::string& pathstr)
{
	std::vector<std::string> result;
	size_t lastStart = 0;
	for (size_t i = 1; i < pathstr.size(); i++)
	{
		if (pathstr[i] != '<' && pathstr[i] != '>') continue;
		result.emplace_back(pathstr.substr(lastStart, i - lastStart));
		lastStart = i;
	}
	result.emplace_back(pathstr.substr(lastStart));
	return result;
}

std::vector<ReadPath> loadReadPaths(const std::string& filename)
{
	std::vector<ReadPath> result;
	std::ifstream file { filename };
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (!file.good()) break;
		std::stringstream sstr { line };
		ReadPath here;
		std::string dummy;
		std::string pathstr;
		size_t pathLength, pathStart, pathEnd;
		sstr >> here.readName >> dummy >> here.readStart >> here.readEnd >> dummy >> pathstr >> pathLength >> pathStart >> pathEnd;
		here.pathStartClip = pathStart;
		here.pathEndClip = pathLength - pathEnd;
		here.path = parsePath(pathstr);
		result.emplace_back(here);
	}
	return result;
}

std::vector<std::string> getReferenceAllele(const Path& heavyPath, const std::string startNode, const std::string endNode)
{
	size_t startIndex = heavyPath.nodes.size();
	size_t endIndex = heavyPath.nodes.size();
	for (size_t i = 0; i < heavyPath.nodes.size(); i++)
	{
		if (heavyPath.nodes[i] == startNode) startIndex = i;
		if (heavyPath.nodes[i] == endNode) endIndex = i;
	}
	assert(startIndex < heavyPath.nodes.size());
	assert(endIndex < heavyPath.nodes.size());
	assert(startIndex != endIndex);
	std::vector<std::string> result;
	if (endIndex > startIndex)
	{
		result.insert(result.end(), heavyPath.nodes.begin() + startIndex, heavyPath.nodes.begin() + endIndex + 1);
	}
	else
	{
		result.insert(result.end(), heavyPath.nodes.begin() + startIndex, heavyPath.nodes.end());
		result.insert(result.end(), heavyPath.nodes.begin(), heavyPath.nodes.begin() + endIndex + 1);
	}
	return result;
}

bool isSubstring(const std::vector<std::string>& small, const std::vector<std::string>& big)
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

size_t getReferenceCoverage(const std::vector<ReadPath>& readPaths, const std::vector<std::string>& referenceAllele)
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
	std::unordered_set<std::string> partOfHeavyPath;
	std::unordered_map<std::string, bool> nodeOrientation;
	for (auto node : heavyPath.nodes)
	{
		partOfHeavyPath.insert(node.substr(1));
		nodeOrientation[node.substr(1)] = node[0] == '>';
	}
	std::map<std::vector<std::string>, size_t> bubbleCoverages;
	for (const auto& path : readPaths)
	{
		size_t lastStart = std::numeric_limits<size_t>::max();
		for (size_t i = 0; i < path.path.size(); i++)
		{
			if (partOfHeavyPath.count(path.path[i].substr(1)) == 0) continue;
			if (lastStart != std::numeric_limits<size_t>::max())
			{
				std::vector<std::string> part { path.path.begin() + lastStart, path.path.begin() + i + 1 };
				part = canon(part);
				bubbleCoverages[part] += 1;
			}
			lastStart = i;
		}
	}
	std::vector<Variant> result;
	for (const auto& pair : bubbleCoverages)
	{
		if (pair.second < minCoverage) continue;
		std::vector<std::string> path = pair.first;
		if (nodeOrientation[path[0].substr(1)] != (path[0][0] == '>'))
		{
			path = reverse(path);
		}
		std::vector<std::string> referenceAllele = getReferenceAllele(heavyPath, path[0], path.back());
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
		for (const auto& node : variant.path)
		{
			file << node;
		}
		file << "\t";
		for (const auto& node : variant.referencePath)
		{
			file << node;
		}
		file << "\t" << variant.coverage << "\t" << variant.referenceCoverage << std::endl;
	}
}

void writeVariantGraph(std::string graphFileName, const Path& heavyPath, const std::vector<Variant>& variants, std::string variantGraphFileName)
{
	std::unordered_set<std::string> allowedNodes;
	std::set<std::pair<std::string, std::string>> allowedEdges;
	for (size_t i = 0; i < heavyPath.nodes.size(); i++)
	{
		std::string prev;
		if (i == 0)
		{
			prev = heavyPath.nodes.back();
		}
		else
		{
			prev = heavyPath.nodes[i-1];
		}
		std::string curr = heavyPath.nodes[i];
		allowedEdges.emplace(prev, curr);
		allowedEdges.emplace(reverse(curr), reverse(prev));
		allowedNodes.emplace(curr.substr(1));
	}
	for (const auto& variant : variants)
	{
		for (size_t i = 0; i < variant.path.size(); i++)
		{
			std::string prev;
			if (i == 0)
			{
				prev = variant.path.back();
			}
			else
			{
				prev = variant.path[i-1];
			}
			std::string curr = variant.path[i];
			allowedEdges.emplace(prev, curr);
			allowedEdges.emplace(reverse(curr), reverse(prev));
			allowedNodes.emplace(curr.substr(1));
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
			std::string node;
			sstr >> node;
			if (allowedNodes.count(node) == 1)
			{
				out << line << std::endl;
			}
		}
		else if (test == "L")
		{
			std::string fromnode, fromorient, tonode, toorient;
			sstr >> fromnode >> fromorient >> tonode >> toorient;
			fromnode = (fromorient == "+" ? ">" : "<") + fromnode;
			tonode = (toorient == "+" ? ">" : "<") + tonode;
			if (allowedEdges.count(std::make_pair(fromnode, tonode)) == 1)
			{
				out << line << std::endl;
			}
		}
	}
}

void HandleCluster(std::string basePath, std::string readPath, std::string MBGPath)
{
	std::cerr << "running MBG" << std::endl;
	runMBG(basePath, readPath, MBGPath);
	std::cerr << "reading graph" << std::endl;
	GfaGraph graph;
	graph.loadFromFile(basePath + "/graph.gfa");
	std::cerr << "getting consensus" << std::endl;
	Path heavyPath = getHeavyPath(graph);
	std::cerr << "writing consensus" << std::endl;
	writePathSequence(heavyPath, graph, basePath + "/consensus.fa");
	writePathGaf(heavyPath, graph, basePath + "/consensus_path.gaf");
	std::cerr << "reading read paths" << std::endl;
	std::vector<ReadPath> readPaths = loadReadPaths(basePath + "/paths.gaf");
	std::cerr << "getting variants" << std::endl;
	std::vector<Variant> variants = getVariants(graph, heavyPath, readPaths, 3);
	std::sort(variants.begin(), variants.end(), [](const Variant& left, const Variant& right) { return left.coverage > right.coverage; });
	std::cerr << "writing variants" << std::endl;
	writeVariants(heavyPath, graph, variants, basePath + "/variants.txt");
	std::cerr << "making variant graph" << std::endl;
	writeVariantGraph(basePath + "/graph.gfa", heavyPath, variants, basePath + "/variant-graph.gfa");
}
