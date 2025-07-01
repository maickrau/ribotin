#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <cassert>
#include <sstream>
#include <phmap.h>
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

std::vector<std::string> parseNodesWithOrientation(std::string pathstr)
{
	std::vector<std::string> result;
	size_t lastStart = 0;
	if (pathstr[0] == '[') lastStart = 1;
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
			lastStart = i;
			if (pathstr[i] == '[') lastStart = i+1;
		}
	}
	if (lastStart < pathstr.size())
	{
		result.emplace_back(pathstr.substr(lastStart));
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

std::unordered_set<std::string> getReadsTouchingNodes(const std::vector<std::string>& nodes, std::string nodemapFile, std::string mbgPathsFile, std::string nodeLengthFile, std::string edgeOverlapFile)
{
	phmap::flat_hash_map<std::string, size_t> nodeLength;
	{
		std::ifstream file { nodeLengthFile };
		while (file.good())
		{
			std::string line;
			getline(file, line);
			if (line.size() == 0) continue;
			std::stringstream sstr { line };
			std::string name;
			size_t length;
			sstr >> name >> length;
			assert(nodeLength.count(name) == 0 || nodeLength.at(name) == length);
			nodeLength[name] = length;
		}
	}
	phmap::flat_hash_map<std::pair<std::string, std::string>, size_t> edgeOverlaps;
	{
		std::ifstream file { edgeOverlapFile };
		while (file.good())
		{
			std::string line;
			getline(file, line);
			if (line.size() == 0) continue;
			std::stringstream sstr { line };
			std::string dummy, fromnodename, fromorient, tonodename, toorient, overlapstr;
			size_t overlap = 0;
			sstr >> dummy >> fromnodename >> fromorient >> tonodename >> toorient >> overlapstr;
			assert(overlapstr.size() >= 2);
			assert(overlapstr.back() == 'M');
			overlap = std::stoull(overlapstr.substr(0, overlapstr.size()-1));
			std::string fromnode = ((fromorient == "+") ? ">" : "<") + fromnodename;
			std::string tonode = ((toorient == "+") ? ">" : "<") + tonodename;
			auto key = std::make_pair(fromnode, tonode);
			assert(edgeOverlaps.count(key) == 0 || edgeOverlaps.at(key) == overlap);
			edgeOverlaps[key] = overlap;
			std::string revFromNode = ((fromorient == "+") ? "<" : ">") + fromnodename;
			std::string revTonode = ((toorient == "+") ? "<" : ">") + tonodename;
			key = std::make_pair(revTonode, revFromNode);
			assert(edgeOverlaps.count(key) == 0 || edgeOverlaps.at(key) == overlap);
			edgeOverlaps[key] = overlap;
		}
	}
	std::unordered_map<std::string, std::tuple<std::vector<std::string>, size_t, size_t>> nodemap;
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
			size_t firstColon = pathstr.find(":");
			assert(pathstr.find(":", firstColon+1) != std::string::npos);
			size_t secondColon = pathstr.find(":", firstColon+1);
			assert(secondColon > firstColon);
			auto nodes = parseNodesWithOrientation(pathstr.substr(0, firstColon));
			size_t leftClip = std::stoull(pathstr.substr(firstColon+1, secondColon-firstColon-1));
			size_t rightClip = std::stoull(pathstr.substr(secondColon+1));
			assert(nodemap.count(key) == 0);
			nodemap[key] = std::make_tuple(nodes, leftClip, rightClip);
		}
	}
	phmap::flat_hash_set<std::tuple<std::string, size_t, size_t>> validNodesWithClip;
	for (auto node : nodes)
	{
		if (node[0] == '>')
		{
			assert(nodeLength.count(node.substr(1)) == 1);
			assert(nodeLength.at(node.substr(1)) >= 50000);
			validNodesWithClip.emplace(node.substr(1), nodeLength.at(node.substr(1))-50000, 0);
		}
		else if (node[0] == '<')
		{
			assert(nodeLength.count(node.substr(1)) == 1);
			assert(nodeLength.at(node.substr(1)) >= 50000);
			validNodesWithClip.emplace(node.substr(1), 0, nodeLength.at(node.substr(1))-50000);
		}
		else
		{
			validNodesWithClip.emplace(node, 0, 0);
		}
	}
	while (true)
	{
		phmap::flat_hash_set<std::tuple<std::string, size_t, size_t>> newValidNodes;
		bool changed = false;
		for (auto t : validNodesWithClip)
		{
			if (nodemap.count(std::get<0>(t)) == 0)
			{
				newValidNodes.emplace(t);
				continue;
			}
			changed = true;
			const auto& vec = std::get<0>(nodemap.at(std::get<0>(t)));
			const size_t leftClip = std::get<1>(nodemap.at(std::get<0>(t))) + std::get<1>(t);
			const size_t rightClip = std::get<2>(nodemap.at(std::get<0>(t))) + std::get<2>(t);
			size_t pathLength = 0;
			for (size_t i = 0; i < vec.size(); i++)
			{
				assert(nodeLength.count(vec[i].substr(1)) == 1);
				pathLength += nodeLength.at(vec[i].substr(1));
				if (i > 0)
				{
					auto key = std::make_pair(vec[i-1], vec[i]);
					assert(edgeOverlaps.count(key) == 1);
					pathLength -= edgeOverlaps.at(key);
				}
			}
			assert(leftClip + rightClip < pathLength);
			size_t nodePathStartPos = 0;
			for (size_t i = 0; i < vec.size(); i++)
			{
				if (i > 0)
				{
					auto key = std::make_pair(vec[i-1], vec[i]);
					assert(edgeOverlaps.count(key) == 1);
					assert(nodePathStartPos >= edgeOverlaps.at(key));
					nodePathStartPos -= edgeOverlaps.at(key);
				}
				size_t nodePathEndPos = nodePathStartPos + nodeLength.at(vec[i].substr(1));
				assert(nodePathEndPos <= pathLength);
				if (nodePathEndPos > leftClip || pathLength-nodePathStartPos > rightClip)
				{
					size_t nodeLeftClip = 0;
					size_t nodeRightClip = 0;
					if (nodePathStartPos < leftClip) nodeLeftClip = leftClip - nodePathStartPos;
					if (pathLength-nodePathEndPos < rightClip) nodeRightClip = rightClip - (pathLength-nodePathEndPos);
					assert(vec[i][0] == '<' || vec[i][0] == '>');
					if (vec[i][0] == '<')
					{
						std::swap(nodeLeftClip, nodeRightClip);
					}
					if (nodeLeftClip + nodeRightClip < nodeLength.at(vec[i].substr(1)))
					{
						newValidNodes.emplace(vec[i].substr(1), nodeLeftClip, nodeRightClip);
					}
				}
				nodePathStartPos += nodeLength.at(vec[i].substr(1));
			}
		}
		if (changed) validNodesWithClip = newValidNodes;
		if (!changed) break;
	}
	phmap::flat_hash_map<std::string, std::vector<std::pair<size_t, size_t>>> validPositionsPerNode;
	for (auto t : validNodesWithClip)
	{
		validPositionsPerNode[std::get<0>(t)].emplace_back(std::get<1>(t), std::get<2>(t));
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
		size_t pathlength, pathstartpos, pathendpos;
		if (line.size() == 0) continue;
		std::stringstream sstr { line };
		sstr >> readname >> dummy >> dummy >> dummy >> dummy >> pathstr >> pathlength >> pathstartpos >> pathendpos;
		auto nodes = parseNodesWithOrientation(pathstr);
		assert(nodes.size() >= 1);
		if (nodes.size() == 1)
		{
			size_t wantedStartClip = pathstartpos;
			size_t wantedEndClip = pathlength-pathendpos;
			if (nodes[0][0] == '<')
			{
				std::swap(wantedStartClip, wantedEndClip);
			}
			if (validPositionsPerNode.count(nodes[0].substr(1)) == 1)
			{
				for (auto pos : validPositionsPerNode.at(nodes[0].substr(1)))
				{
					if (wantedStartClip < pathlength-pos.second && wantedEndClip < pathlength-pos.first)
					{
						result.insert(readname);
						break;
					}
				}
			}
		}
		else
		{
			for (size_t i = 1; i+1 < nodes.size(); i++)
			{
				if (validPositionsPerNode.count(nodes[i].substr(1)) == 1)
				{
					result.insert(readname);
					break;
				}
			}
			if (validPositionsPerNode.count(nodes[0].substr(1)) == 1)
			{
				size_t wantedStartClip = pathstartpos;
				size_t nodelen = nodeLength.at(nodes[0].substr(1));
				for (auto pos : validPositionsPerNode.at(nodes[0].substr(1)))
				{
					if (wantedStartClip < nodelen-pos.second)
					{
						result.insert(readname);
						break;
					}
				}
			}
			if (validPositionsPerNode.count(nodes.back().substr(1)) == 1)
			{
				size_t wantedEndClip = pathlength-pathendpos;
				size_t nodelen = nodeLength.at(nodes.back().substr(1));
				for (auto pos : validPositionsPerNode.at(nodes.back().substr(1)))
				{
					if (wantedEndClip < nodelen-pos.first)
					{
						result.insert(readname);
						break;
					}
				}
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
			if (node[0] == '>' || node[0] == '<') continue;
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
	for (size_t i = 0; i < nodesPerTangle.size(); i++)
	{
		for (auto node : nodesPerTangle[i])
		{
			if (node[0] != '>' && node[0] != '<') continue;
			std::string checknode = node.substr(1);
			if (pathsPerNode.count(checknode) == 0) continue;
			for (auto path : pathsPerNode.at(checknode))
			{
				if (result.count(path) == 1)
				{
					if (result.at(path) >= nodesPerTangle.size())
					{
						if (result.at(path) != i + nodesPerTangle.size())
						{
							result[path] = std::numeric_limits<size_t>::max();
						}
					}
				}
				else
				{
					result[path] = i + nodesPerTangle.size();
				}
			}
		}
	}
	for (auto& pair : result)
	{
		if (pair.second >= nodesPerTangle.size()) pair.second -= nodesPerTangle.size();
	}
	return result;
}

std::vector<std::vector<std::string>> getReadNamesPerTangle(std::string verkkoBaseFolder, const std::vector<std::vector<std::string>>& nodesPerTangle)
{
	std::string scfMapFile = verkkoBaseFolder + "/6-layoutContigs/unitig-popped.layout.scfmap";
	std::string layoutFile = verkkoBaseFolder + "/6-layoutContigs/unitig-popped.layout";
	std::string nodemapFile = verkkoBaseFolder + "/6-layoutContigs/combined-nodemap.txt";
	std::string nodeLengthFile = verkkoBaseFolder + "/6-layoutContigs/nodelens.txt";
	std::string edgeOverlapFile = verkkoBaseFolder + "/6-layoutContigs/combined-edges.gfa";
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
		auto potentialReads = getReadsTouchingNodes(nodesPerTangle[i], nodemapFile, hifipathFile, nodeLengthFile, edgeOverlapFile);
		for (auto read : uniqueReadPerTangle[i])
		{
			if (potentialReads.count(read) == 0) continue;
			result[i].push_back(read);
		}
	}
	return result;
}
