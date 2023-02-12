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

std::unordered_map<std::string, std::string> getUniqueReadAssignmentToPiece(std::string layoutFile)
{
	std::unordered_map<std::string, std::string> result;
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
		if (result.count(test) == 1)
		{
			if (result.at(test) != currentTig)
			{
				result[test] = "";
			}
		}
		else
		{
			result[test] = currentTig;
		}
	}
	return result;
}

std::unordered_map<std::string, std::string> getUniquePieceAssignmentToPath(std::string scfMapFile)
{
	std::unordered_map<std::string, std::string> result;
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
		if (result.count(test) == 1)
		{
			if (result.at(test) != currentPath)
			{
				result[test] = "";
			}
		}
		else
		{
			result[test] = currentPath;
		}
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

std::unordered_map<std::string, std::string> getUniqueNodeAssignmentToPath(std::string pathsFile)
{
	std::unordered_map<std::string, std::string> result;
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
			if (result.count(node) == 1)
			{
				if (result.at(node) != pathname)
				{
					result[node] = "";
				}
			}
			else
			{
				result[node] = pathname;
			}
		}
	}
	return result;
}

std::unordered_map<std::string, std::vector<std::string>> getUniqueReadAssignmentToPath(const std::unordered_map<std::string, std::string>& uniquePiecePerRead, const std::unordered_map<std::string, std::string>& uniquePathPerPiece)
{
	std::unordered_map<std::string, std::vector<std::string>> readsPerPath;
	for (auto pair : uniquePiecePerRead)
	{
		if (pair.second == "") continue;
		if (uniquePathPerPiece.count(pair.second) == 0) continue;
		auto path = uniquePathPerPiece.at(pair.second);
		if (path == "") continue;
		readsPerPath[path].push_back(pair.first);
	}
	return readsPerPath;
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
	auto uniquePiecePerRead = getUniqueReadAssignmentToPiece(layoutFile);
	auto uniquePathPerPiece = getUniquePieceAssignmentToPath(scfMapFile);
	auto uniquePathPerNode = getUniqueNodeAssignmentToPath(pathsFile);
	auto uniqueReadPerPath = getUniqueReadAssignmentToPath(uniquePiecePerRead, uniquePathPerPiece);
	auto uniqueNodesPerCluster = getUniqueNodesPerCluster(nodesPerCluster);
	std::vector<std::vector<std::string>> result;
	result.resize(nodesPerCluster.size());
	for (size_t i = 0; i < nodesPerCluster.size(); i++)
	{
		auto potentialReads = getReadsTouchingNodes(nodesPerCluster[i], verkkoBaseFolder + "/6-layoutContigs/combined-nodemap.txt", verkkoBaseFolder + "/1-buildGraph/paths.gaf");
		for (const auto& node : nodesPerCluster[i])
		{
			if (uniquePathPerNode.count(node) == 0) continue;
			if (uniquePathPerNode.at(node) == "") continue;
			if (uniqueNodesPerCluster.count(node) == 0) continue;
			auto path = uniquePathPerNode.at(node);
			if (uniqueReadPerPath.count(path) == 0) continue;
			for (auto read : uniqueReadPerPath.at(path))
			{
				if (potentialReads.count(read) == 0) continue;
				result[i].push_back(read);
			}
		}
	}
	return result;
}
