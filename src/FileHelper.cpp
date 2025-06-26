#include <fstream>
#include "FileHelper.h"
#include "RibotinUtils.h"

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

Path readHeavyPath(const GfaGraph& graph, const std::string& heavyPathGafFile)
{
	std::ifstream file { heavyPathGafFile };
	std::string line;
	getline(file, line);
	auto parts = split(line, '\t');
	Path result;
	size_t pathLength = std::stoull(parts[6]);
	size_t pathStart = std::stoull(parts[7]);
	size_t pathEnd = std::stoull(parts[8]);
	result.leftClip = pathStart;
	result.rightClip = pathLength - pathEnd;
	result.nodes = parsePath(parts[5], graph.nodeNameToId);
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
			if (prev == result.nodes[0])
			{
				assert(result.nodes.size() >= 3);
				prev = result.nodes[result.nodes.size()-2];
			}
		}
		assert(graph.edges.count(prev) == 1);
		for (auto edge : graph.edges.at(prev))
		{
			if (std::get<0>(edge) != result.nodes[i]) continue;
			result.overlaps.push_back(std::get<1>(edge));
		}
	}
	assert(result.getSequence(graph.nodeSeqs).size() == pathEnd - pathStart);
	return result;
}

void writePathSequence(const Path& path, const GfaGraph& graph, std::string outputFile, const std::string& namePrefix)
{
	std::string pathseq = path.getSequence(graph.nodeSeqs);
	std::ofstream file { outputFile };
	file << ">";
	if (namePrefix != "") file << namePrefix << "_";
	file << "heavy_path" << std::endl;
	file << pathseq << std::endl;
}

void writePathGaf(const Path& path, const GfaGraph& graph, std::string outputFile)
{
	std::ofstream file { outputFile };
	file << getPathGaf(path, graph);
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
		size_t readlength = std::stoull(parts[1]);
		size_t readstart = std::stoull(parts[2]);
		size_t readend = std::stoull(parts[3]);
		std::string pathstr = parts[5];
		size_t pathLength = std::stoull(parts[6]);
		size_t pathStart = std::stoull(parts[7]);
		size_t pathEnd = std::stoull(parts[8]);
		size_t mapq = std::stoull(parts[11]);
		here.readName = readname;
		here.readLength = readlength;
		here.readStart = readstart;
		here.readEnd = readend;
		here.pathStartClip = pathStart;
		here.pathEndClip = pathLength - pathEnd;
		here.path = parsePath(pathstr, graph.nodeNameToId);
		here.mapq = mapq;
		assert(here.path.size() >= 1);
		result.emplace_back(here);
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
		assert(result.count(str.id()) == 0 || result.at(str.id()) == str.forward());
		result[str.id()] = str.forward();
	}
	return result;
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
		if (reverse) std::swap(readPath.pathStartClip, readPath.pathEndClip);
		result.push_back(readPath);
	}
	return result;
}
