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
	Path result;
	result.nodes.emplace_back(">" + maxCoverageNode);
	result.overlaps.emplace_back(0);
	std::unordered_set<std::string> visited;
	while (true)
	{
		std::string bestNeighbor;
		size_t bestNeighborOverlap = 0;
		double bestNeighborCoverage = 0;
		for (auto edge : graph.edges.at(result.nodes.back()))
		{
			std::string neighbor = std::get<0>(edge);
			double neighborCoverage = std::min(graph.nodeCoverages.at(neighbor.substr(1)), std::get<2>(edge));
			if (bestNeighbor == "" || neighborCoverage > bestNeighborCoverage)
			{
				bestNeighbor = neighbor;
				bestNeighborOverlap = std::get<1>(edge);
				bestNeighborCoverage = neighborCoverage;
			}
		}
		if (bestNeighbor.substr(1) == maxCoverageNode)
		{
			result.overlaps[0] = bestNeighborOverlap;
			break;
		}
		assert(visited.count(bestNeighbor) == 0);
		visited.insert(bestNeighbor);
		result.nodes.emplace_back(bestNeighbor);
		result.overlaps.emplace_back(bestNeighborOverlap);
	}
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

void HandleCluster(std::string basePath, std::string readPath, std::string MBGPath)
{
	std::cerr << "running MBG" << std::endl;
	runMBG(basePath, readPath, MBGPath);
	std::cerr << "getting consensus" << std::endl;
	GfaGraph graph;
	graph.loadFromFile(basePath + "/graph.gfa");
	Path heavyPath = getHeavyPath(graph);
	writePathSequence(heavyPath, graph, basePath + "/consensus.fa");
	writePathGaf(heavyPath, graph, basePath + "/consensus_path.gaf");
}
