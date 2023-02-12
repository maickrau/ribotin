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

std::string getPathString(const std::unordered_map<std::string, std::string>& nodeSeqs, const std::vector<std::pair<std::string, size_t>>& path)
{
	std::string result;
	for (size_t i = 0; i < path.size(); i++)
	{
		std::string add;
		add = nodeSeqs.at(path[i].first.substr(1));
		if (path[i].first[0] == '<')
		{
			add = revcomp(add);
		}
		else
		{
			assert(path[i].first[0] == '>');
		}
		add = add.substr(path[i].second);
		result += add;
	}
	return result;
}

std::string getPathGaf(const std::string& pathseq, const std::vector<std::pair<std::string, size_t>>& path)
{
	std::string result = "heavy_path\t";
	result += std::to_string(pathseq.size());
	result += "\t0\t";
	result += std::to_string(pathseq.size());
	result += "\t+\t";
	for (const auto& pair : path)
	{
		result += pair.first;
	}
	result += "\t";
	result += std::to_string(pathseq.size());
	result += std::to_string(path[0].second);
	result += "\t";
	result += std::to_string(path[0].second);
	result += "\t";
	result += std::to_string(pathseq.size());
	result += std::to_string(path[0].second);
	result += "\t";
	result += std::to_string(pathseq.size());
	result += "\t";
	result += std::to_string(pathseq.size());
	result += "\t60";
	return result;
}

void getHeavyPath(std::string gfaFile, std::string outputPathFile, std::string outputPathGaf)
{
	std::unordered_map<std::string, size_t> nodeCoverages;
	std::unordered_map<std::string, std::string> nodeSeqs;
	std::unordered_map<std::string, std::set<std::tuple<std::string, size_t, size_t>>> edges;
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
	std::string maxCoverageNode;
	for (auto pair : nodeCoverages)
	{
		if (maxCoverageNode == "" || pair.second > nodeCoverages.at(maxCoverageNode))
		{
			maxCoverageNode = pair.first;
		}
	}
	std::vector<std::pair<std::string, size_t>> path;
	path.emplace_back(">" + maxCoverageNode, 0);
	std::unordered_set<std::string> visited;
	while (true)
	{
		std::string bestNeighbor;
		size_t bestNeighborOverlap = 0;
		double bestNeighborCoverage = 0;
		for (auto edge : edges[path.back().first])
		{
			std::string neighbor = std::get<0>(edge);
			double neighborCoverage = std::min(nodeCoverages.at(neighbor.substr(1)), std::get<2>(edge));
			if (bestNeighbor == "" || neighborCoverage > bestNeighborCoverage)
			{
				bestNeighbor = neighbor;
				bestNeighborOverlap = std::get<1>(edge);
				bestNeighborCoverage = neighborCoverage;
			}
		}
		if (bestNeighbor.substr(1) == maxCoverageNode)
		{
			path[0].second = bestNeighborOverlap;
			break;
		}
		assert(visited.count(bestNeighbor) == 0);
		visited.insert(bestNeighbor);
		path.emplace_back(bestNeighbor, bestNeighborOverlap);
	}
	std::string pathseq = getPathString(nodeSeqs, path);
	{
		std::ofstream file { outputPathFile };
		file << ">heavy_path" << std::endl;
		file << pathseq;
	}
	{
		std::ofstream file { outputPathGaf };
		file << getPathGaf(pathseq, path);
	}
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
	getHeavyPath(basePath + "/graph.gfa", basePath + "/consensus.fa", basePath + "/consensus_path.gaf");
}
