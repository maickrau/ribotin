#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <cassert>
#include <sstream>
#include "VerkkoReadAssignment.h"
#include "RibotinUtils.h"
#include "Logger.h"

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
	// assume missing value is false
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

std::vector<std::vector<std::string>> getUniqueReadAssignmentToTangle(const std::unordered_map<std::string, std::unordered_set<std::string>>& piecesPerRead, const std::unordered_map<std::string, std::unordered_set<std::string>>& pathsPerPiece, const std::unordered_map<std::string, size_t>& uniqueTanglePerPath, size_t numTangles)
{
	std::unordered_map<std::string, size_t> tanglePerRead;
	for (auto pair : piecesPerRead)
	{
		for (auto piece : pair.second)
		{
			if (pathsPerPiece.count(piece) == 0) continue;
			for (auto path : pathsPerPiece.at(piece))
			{
				if (uniqueTanglePerPath.count(path) == 0) continue;
				size_t tangle = uniqueTanglePerPath.at(path);
				if (tanglePerRead.count(pair.first) == 1)
				{
					if (tanglePerRead.at(pair.first) != tangle)
					{
						tanglePerRead[pair.first] = std::numeric_limits<size_t>::max();
					}
				}
				else
				{
					tanglePerRead[pair.first] = tangle;
				}
			}
		}
	}
	std::vector<std::vector<std::string>> result;
	result.resize(numTangles);
	for (auto pair : tanglePerRead)
	{
		if (pair.second < numTangles)
		{
			result[pair.second].push_back(pair.first);
		}
	}
	return result;
}

std::unordered_set<std::string> getUniqueNodesPerTangle(const std::vector<std::vector<std::string>>& nodesPerTangle)
{
	std::unordered_set<std::string> unique;
	std::unordered_set<std::string> notUnique;
	for (const auto& tangle : nodesPerTangle)
	{
		for (const auto& node : tangle)
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

std::unordered_map<std::string, size_t> getUniquePathAssignmentToTangle(const std::unordered_map<std::string, std::unordered_set<std::string>>& pathsPerNode, const std::vector<std::vector<std::string>>& nodesPerTangle)
{
	std::unordered_map<std::string, size_t> result;
	for (size_t i = 0; i < nodesPerTangle.size(); i++)
	{
		for (auto node : nodesPerTangle[i])
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

std::vector<std::vector<std::string>> getReadNamesPerTangle(std::string verkkoBaseFolder, const std::vector<std::vector<std::string>>& nodesPerTangle)
{
	std::string scfMapFile = verkkoBaseFolder + "/6-layoutContigs/unitig-popped.layout.scfmap";
	std::string layoutFile = verkkoBaseFolder + "/6-layoutContigs/unitig-popped.layout";
	std::string nodemapFile = verkkoBaseFolder + "/6-layoutContigs/combined-nodemap.txt";
	std::string hifipathFile = verkkoBaseFolder + "/1-buildGraph/paths.gaf";
	if (!fileExists(scfMapFile))
	{
		Logger::Log.log(Logger::LogLevel::Always) << "ERROR: could not find layout scfmap file in the verkko assembly folder!" << std::endl;
		std::abort();
	}
	if (!fileExists(layoutFile))
	{
		Logger::Log.log(Logger::LogLevel::Always) << "ERROR: could not find layout file in the verkko assembly folder!" << std::endl;
		std::abort();
	}
	if (!fileExists(nodemapFile))
	{
		Logger::Log.log(Logger::LogLevel::Always) << "ERROR: could not find node map file in the verkko assembly folder!" << std::endl;
		std::abort();
	}
	if (!fileExists(hifipathFile))
	{
		Logger::Log.log(Logger::LogLevel::Always) << "ERROR: could not find hifi path file in the verkko assembly folder!" << std::endl;
		std::abort();
	}
	std::string pathsFile;
	if (getRukkiEnabled(verkkoBaseFolder + "/verkko.yml"))
	{
		pathsFile = verkkoBaseFolder + "/6-rukki/unitig-popped-unitig-normal-connected-tip.paths.gaf";
		if (!fileExists(pathsFile))
		{
			pathsFile = verkkoBaseFolder + "/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.paths.gaf";
			if (!fileExists(pathsFile))
			{
				Logger::Log.log(Logger::LogLevel::Always) << "ERROR: could not find rukki paths file in the verkko assembly folder!" << std::endl;
				std::abort();
			}
		}
	}
	else
	{
		pathsFile = verkkoBaseFolder + "/6-layoutContigs/consensus_paths.txt";
	}
	if (!fileExists(pathsFile))
	{
		Logger::Log.log(Logger::LogLevel::Always) << "ERROR: could not find paths file in the verkko assembly folder!" << std::endl;
		std::abort();
	}
	auto piecesPerRead = getReadAssignmentToPieces(layoutFile);
	auto pathsPerPiece = getPieceAssignmentToPaths(scfMapFile);
	auto pathsPerNode = getNodeAssignmentToPaths(pathsFile);
	auto uniqueTanglePerPath = getUniquePathAssignmentToTangle(pathsPerNode, nodesPerTangle);
	auto uniqueReadPerTangle = getUniqueReadAssignmentToTangle(piecesPerRead, pathsPerPiece, uniqueTanglePerPath, nodesPerTangle.size());
	std::vector<std::vector<std::string>> result;
	result.resize(nodesPerTangle.size());
	for (size_t i = 0; i < nodesPerTangle.size(); i++)
	{
		auto potentialReads = getReadsTouchingNodes(nodesPerTangle[i], nodemapFile, hifipathFile);
		for (auto read : uniqueReadPerTangle[i])
		{
			if (potentialReads.count(read) == 0) continue;
			result[i].push_back(read);
		}
	}
	return result;
}
