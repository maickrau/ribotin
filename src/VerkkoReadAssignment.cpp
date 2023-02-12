#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <cassert>
#include <sstream>
#include "VerkkoReadAssignment.h"

bool getRukkiEnabled(std::string configfile)
{
	// not a good way but works for this one narrow use case
	std::ifstream file { configfile };
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (line.size() >= 10 && line.substr(0, 10) == "ruk_enable")
		{
			bool hasTrue = line.find("True") != std::string::npos;
			bool hasFalse = line.find("False") != std::string::npos;
			assert(hasTrue || hasFalse);
			assert(!(hasTrue && hasFalse));
			return hasTrue;
		}
	}
	assert(false);
	return false;
}

std::unordered_map<std::string, std::unordered_set<std::string>> getReadAssignmentToPieces(std::string layoutFile)
{
	std::unordered_map<std::string, std::unordered_set<std::string>> result;
	std::ifstream file { layoutFile };
	std::string currentTig;
	while (file.good())
	{
		std::string line;
		getline(file, line);
		std::stringstream sstr { line };
		std::string test;
		sstr >> test;
		if (test == "tig")
		{
			sstr >> currentTig;
			continue;
		}
		if (test == "len") continue;
		if (test == "rds") continue;
		if (test == "end") continue;
		result[test].insert(currentTig);
	}
	return result;
}

std::unordered_map<std::string, std::unordered_set<std::string>> getPieceAssignmentToPaths(std::string scfMapFile)
{
	std::unordered_map<std::string, std::unordered_set<std::string>> result;
	std::ifstream file { scfMapFile };
	std::string currentPath;
	while (file.good())
	{
		std::string line;
		getline(file, line);
		std::stringstream sstr { line };
		std::string test;
		sstr >> test;
		if (test == "path")
		{
			sstr >> test;
			sstr >> currentPath;
			continue;
		}
		if (test == "end") continue;
		if (test.size() > 0 && test[0] == '[') continue;
		result[test].insert(currentPath);
	}
	return result;
}

std::vector<std::string> parseNodes(std::string pathstr)
{
	std::vector<std::string> result;
	size_t lastStart = 1;
	//skip first because it's one of < > [
	for (size_t i = 1; i < pathstr.size(); i++)
	{
		if (pathstr[i] == ']')
		{
			lastStart = i+1;
			continue;
		}
		if (pathstr[i] == '[' || pathstr[i] == '<' || pathstr[i] == '>')
		{
			if (i > lastStart)
			{
				result.emplace_back(pathstr.substr(lastStart, i-lastStart));
			}
			lastStart = i+1;
		}
	}
	if (lastStart < pathstr.size())
	{
		result.emplace_back(pathstr.substr(lastStart));
	}
	return result;
}

std::unordered_map<std::string, std::unordered_set<std::string>> getNodeAssignmentToPaths(std::string pathsFile)
{
	std::unordered_map<std::string, std::unordered_set<std::string>> result;
	std::ifstream file { pathsFile };
	while (file.good())
	{
		std::string line;
		getline(file, line);
		std::stringstream sstr { line };
		std::string pathname;
		std::string pathstr;
		sstr >> pathname >> pathstr;
		if (pathname == "name") continue;
		auto nodes = parseNodes(pathstr);
		for (auto node : nodes)
		{
			result[node].insert(pathname);
		}
	}
	return result;
}

std::vector<std::vector<std::string>> getUniqueReadAssignmentToCluster(const std::unordered_map<std::string, std::unordered_set<std::string>>& piecesPerRead, const std::unordered_map<std::string, std::unordered_set<std::string>>& pathsPerPiece, const std::unordered_map<std::string, size_t>& uniqueClusterPerPath, size_t numClusters)
{
	std::unordered_map<std::string, size_t> clusterPerRead;
	for (auto pair : piecesPerRead)
	{
		for (auto piece : pair.second)
		{
			if (pathsPerPiece.count(piece) == 0) continue;
			for (auto path : pathsPerPiece.at(piece))
			{
				if (uniqueClusterPerPath.count(path) == 0) continue;
				size_t cluster = uniqueClusterPerPath.at(path);
				if (clusterPerRead.count(pair.first) == 1)
				{
					if (clusterPerRead.at(pair.first) != cluster)
					{
						clusterPerRead[pair.first] = std::numeric_limits<size_t>::max();
					}
				}
				else
				{
					clusterPerRead[pair.first] = cluster;
				}
			}
		}
	}
	std::vector<std::vector<std::string>> result;
	result.resize(numClusters);
	for (auto pair : clusterPerRead)
	{
		if (pair.second < numClusters)
		{
			result[pair.second].push_back(pair.first);
		}
	}
	return result;
}

std::unordered_set<std::string> getUniqueNodesPerCluster(const std::vector<std::vector<std::string>>& nodesPerCluster)
{
	std::unordered_set<std::string> unique;
	std::unordered_set<std::string> notUnique;
	for (const auto& cluster : nodesPerCluster)
	{
		for (const auto& node : cluster)
		{
			if (notUnique.count(node) == 1) continue;
			if (unique.count(node) == 1)
			{
				notUnique.insert(node);
				unique.erase(node);
			}
			else
			{
				unique.insert(node);
			}
		}
	}
	return unique;
}

std::unordered_set<std::string> getReadsTouchingNodes(const std::vector<std::string>& nodes, std::string nodemapFile, std::string mbgPathsFile)
{
	std::unordered_map<std::string, std::vector<std::string>> nodemap;
	{
		std::ifstream file { nodemapFile };
		while (file.good())
		{
			std::string line;
			getline(file, line);
			if (line.size() == 0) continue;
			std::stringstream sstr { line };
			std::string key;
			std::string pathstr;
			sstr >> key >> pathstr;
			assert(pathstr.find(":") != std::string::npos);
			pathstr.erase(pathstr.begin() + pathstr.find(":"), pathstr.end());
			auto nodes = parseNodes(pathstr);
			assert(nodemap.count(key) == 0);
			nodemap[key] = nodes;
		}
	}
	std::unordered_set<std::string> validNodes { nodes.begin(), nodes.end() };
	while (true)
	{
		std::unordered_set<std::string> newValidNodes;
		bool changed = false;
		for (auto node : validNodes)
		{
			if (nodemap.count(node) == 0)
			{
				newValidNodes.insert(node);
				continue;
			}
			else
			{
				newValidNodes.insert(nodemap.at(node).begin(), nodemap.at(node).end());
				changed = true;
			}
		}
		if (changed) validNodes = newValidNodes;
		if (!changed) break;
	}
	std::ifstream file { mbgPathsFile };
	std::unordered_set<std::string> result;
	while (file.good())
	{
		std::string readname;
		std::string pathstr;
		std::string dummy;
		std::string line;
		getline(file, line);
		if (line.size() == 0) continue;
		std::stringstream sstr { line };
		sstr >> readname >> dummy >> dummy >> dummy >> dummy >> pathstr;
		auto nodes = parseNodes(pathstr);
		for (auto node : nodes)
		{
			if (validNodes.count(node) == 1)
			{
				result.insert(readname);
				break;
			}
		}
	}
	return result;
}

std::unordered_map<std::string, size_t> getUniquePathAssignmentToCluster(const std::unordered_map<std::string, std::unordered_set<std::string>>& pathsPerNode, const std::vector<std::vector<std::string>>& nodesPerCluster)
{
	std::unordered_map<std::string, size_t> result;
	for (size_t i = 0; i < nodesPerCluster.size(); i++)
	{
		for (auto node : nodesPerCluster[i])
		{
			if (pathsPerNode.count(node) == 0) continue;
			for (auto path : pathsPerNode.at(node))
			{
				if (result.count(path) == 1)
				{
					if (result.at(path) != i)
					{
						result[path] = std::numeric_limits<size_t>::max();
					}
				}
				else
				{
					result[path] = i;
				}
			}
		}
	}
	return result;
}

std::vector<std::vector<std::string>> getReadNamesPerCluster(std::string verkkoBaseFolder, const std::vector<std::vector<std::string>>& nodesPerCluster)
{
	std::string scfMapFile = verkkoBaseFolder + "/6-layoutContigs/unitig-popped.layout.scfmap";
	std::string layoutFile = verkkoBaseFolder + "/6-layoutContigs/unitig-popped.layout";
	std::string pathsFile;
	if (getRukkiEnabled(verkkoBaseFolder + "/verkko.yml"))
	{
		pathsFile = verkkoBaseFolder + "/6-rukki/unitig-popped-unitig-normal-connected-tip.paths.gaf";
	}
	else
	{
		pathsFile = verkkoBaseFolder + "/6-layoutContigs/consensus_paths.txt";
	}
	auto piecesPerRead = getReadAssignmentToPieces(layoutFile);
	auto pathsPerPiece = getPieceAssignmentToPaths(scfMapFile);
	auto pathsPerNode = getNodeAssignmentToPaths(pathsFile);
	auto uniqueClusterPerPath = getUniquePathAssignmentToCluster(pathsPerNode, nodesPerCluster);
	auto uniqueReadPerCluster = getUniqueReadAssignmentToCluster(piecesPerRead, pathsPerPiece, uniqueClusterPerPath, nodesPerCluster.size());
	std::vector<std::vector<std::string>> result;
	result.resize(nodesPerCluster.size());
	for (size_t i = 0; i < nodesPerCluster.size(); i++)
	{
		auto potentialReads = getReadsTouchingNodes(nodesPerCluster[i], verkkoBaseFolder + "/6-layoutContigs/combined-nodemap.txt", verkkoBaseFolder + "/1-buildGraph/paths.gaf");
		for (auto read : uniqueReadPerCluster[i])
		{
			if (potentialReads.count(read) == 0) continue;
			result[i].push_back(read);
		}
	}
	return result;
}
