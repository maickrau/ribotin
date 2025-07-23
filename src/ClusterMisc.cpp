#include <sstream>
#include <map>
#include <fstream>
#include <set>
#include "ClusterMisc.h"
#include "fastqloader.h"
#include "Logger.h"

Node::Node(size_t id, bool forward) :
	value(id + (forward ? 0x8000000000000000 : 0))
{
}

size_t Node::id() const
{
	return value & 0x7FFFFFFFFFFFFFFF;
}

bool Node::forward() const
{
	return (value & 0x8000000000000000) == 0x8000000000000000;
}

bool Node::operator<(const Node& other) const
{
	return value < other.value;
}

bool Node::operator==(const Node& other) const
{
	return value == other.value;
}

bool Node::operator!=(const Node& other) const
{
	return value != other.value;
}

size_t Node::rawValue() const
{
	return value;
}

std::vector<char> getComplement()
{
	std::vector<char> result;
	result.resize(256, 0);
	result['a'] = 'T';
	result['A'] = 'T';
	result['c'] = 'G';
	result['C'] = 'G';
	result['g'] = 'C';
	result['G'] = 'C';
	result['t'] = 'A';
	result['T'] = 'A';
	return result;
}

std::vector<char> complement = getComplement();

std::string revcomp(std::string seq)
{
	std::reverse(seq.begin(), seq.end());
	for (size_t i = 0; i < seq.size(); i++)
	{
		seq[i] = complement[seq[i]];
	}
	return seq;
}

Node reverse(Node node)
{
	return Node { node.id(), !node.forward() };
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

std::vector<Node> reverse(const std::vector<Node>& path)
{
	std::vector<Node> bw { path.rbegin(), path.rend() };
	for (size_t i = 0; i < bw.size(); i++)
	{
		bw[i] = reverse(bw[i]);
	}
	return bw;
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

std::pair<std::vector<Node>, std::vector<Node>> canon(const std::vector<Node>& left, const std::vector<Node>& right)
{
	if (left < right) return std::make_pair(left, right);
	return std::make_pair(right, left);
}

std::pair<std::vector<std::string>, std::vector<std::string>> canon(const std::vector<std::string>& left, const std::vector<std::string>& right)
{
	if (left < right) return std::make_pair(left, right);
	return std::make_pair(right, left);
}

std::vector<Node> canon(const std::vector<Node>& path)
{
	auto bw = reverse(path);
	if (bw < path) return bw;
	return path;
}

std::vector<std::string> canon(const std::vector<std::string>& path)
{
	auto bw = reverse(path);
	if (bw < path) return bw;
	return path;
}

std::pair<Node, Node> canon(Node left, Node right)
{
	if (reverse(right) < left) return std::make_pair(reverse(right), reverse(left));
	if (left < reverse(right)) return std::make_pair(left, right);
	if (reverse(left) < right) return std::make_pair(reverse(right), reverse(left));
	return std::make_pair(left, right);
}

std::pair<std::string, std::string> canon(const std::string& left, const std::string& right)
{
	auto revleft = reverse(left);
	auto revright = reverse(right);
	if (revright < left) return std::make_pair(revright, revleft);
	if (left < revright) return std::make_pair(left, right);
	if (revleft < right) return std::make_pair(revright, revleft);
	return std::make_pair(left, right);
}

size_t find(std::vector<size_t>& parent, size_t key)
{
	while (parent.at(key) != parent.at(parent.at(key))) parent[key] = parent[parent[key]];
	return parent[key];
}

void merge(std::vector<size_t>& parent, size_t left, size_t right)
{
	left = find(parent, left);
	right = find(parent, right);
	assert(parent.at(left) == left);
	assert(parent.at(right) == right);
	parent[right] = left;
}

void GfaGraph::loadFromFile(std::string gfaFile)
{
	std::ifstream file { gfaFile };
	std::set<std::tuple<Node, Node, size_t, size_t>> fileedges;
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
			if (nodeNameToId.count(node) == 0)
			{
				size_t id = nodeNameToId.size();
				nodeNameToId[node] = id;
				nodeSeqs.emplace_back();
				revCompNodeSeqs.emplace_back();
				nodeCoverages.emplace_back(0);
				nodeNames.push_back(node);
			}
			nodeSeqs[nodeNameToId[node]] = seq;
			revCompNodeSeqs[nodeNameToId[node]] = revcomp(seq);
			nodeCoverages[nodeNameToId[node]] = std::stof(coveragetag.substr(5));
		}
		else if (test == "L")
		{
			std::string fromnodename;
			std::string tonodename;
			std::string fromorient, toorient;
			std::string overlapstr, edgecoveragetag;
			sstr >> fromnodename >> fromorient >> tonodename >> toorient >> overlapstr >> edgecoveragetag;
			if (nodeNameToId.count(fromnodename) == 0)
			{
				size_t id = nodeNameToId.size();
				nodeNameToId[fromnodename] = id;
				nodeSeqs.emplace_back();
				revCompNodeSeqs.emplace_back();
				nodeCoverages.emplace_back(0);
				nodeNames.push_back(fromnodename);
			}
			if (nodeNameToId.count(tonodename) == 0)
			{
				size_t id = nodeNameToId.size();
				nodeNameToId[tonodename] = id;
				nodeSeqs.emplace_back();
				revCompNodeSeqs.emplace_back();
				nodeCoverages.emplace_back(0);
				nodeNames.push_back(tonodename);
			}
			Node fromnode, tonode;
			if (fromorient == "+")
			{
				fromnode = Node { nodeNameToId[fromnodename], true };
			}
			else
			{
				fromnode = Node { nodeNameToId[fromnodename], false };
			}
			if (toorient == "+")
			{
				tonode = Node { nodeNameToId[tonodename], true };
			}
			else
			{
				tonode = Node { nodeNameToId[tonodename], false };
			}
			overlapstr.pop_back();
			size_t overlap = std::stoull(overlapstr);
			size_t coverage = std::stoull(edgecoveragetag.substr(5));
			fileedges.emplace(fromnode, tonode, overlap, coverage);
			fileedges.emplace(reverse(tonode), reverse(fromnode), overlap, coverage);
		}
	}
	for (auto t : fileedges)
	{
		edges[std::get<0>(t)].emplace(std::get<1>(t), std::get<2>(t), std::get<3>(t));
	}
	for (size_t i = 0; i < nodeSeqs.size(); i++)
	{
		nodeSeqsTwobit.emplace_back(nodeSeqs[i]);
		revCompNodeSeqsTwobit.emplace_back(revCompNodeSeqs[i]);
	}
}

size_t GfaGraph::numNodes() const
{
	return nodeSeqs.size();
}

Path::Path() :
	nodes(),
	overlaps(),
	leftClip(0),
	rightClip(0)
{
}

std::string Path::getSequence(const std::vector<std::string>& nodeSeqs) const
{
	std::string result;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		std::string add;
		add = nodeSeqs.at(nodes[i].id());
		if (!nodes[i].forward())
		{
			add = revcomp(add);
		}
		if (i > 0) add = add.substr(overlaps[i]);
		result += add;
	}
	result = result.substr(leftClip, result.size() - leftClip - rightClip);
	return result;
}

PathSequenceView::PathSequenceView(const std::vector<Node>& nodes, const std::vector<TwobitString>& nodeSeqs, const std::vector<TwobitString>& revCompNodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges, size_t leftClip, size_t rightClip)
{
	size_t realPos = 0;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		if (nodes[i].forward())
		{
			nodePointers.emplace_back(&(nodeSeqs[nodes[i].id()]));
		}
		else
		{
			nodePointers.emplace_back(&(revCompNodeSeqs[nodes[i].id()]));
		}
		size_t overlap = 0;
		if (i > 0)
		{
			overlap = std::numeric_limits<size_t>::max();
			for (auto edge : edges.at(nodes[i-1]))
			{
				if (std::get<0>(edge) != nodes[i]) continue;
				assert(overlap == std::numeric_limits<size_t>::max());
				overlap = std::get<1>(edge);
			}
			assert(overlap != std::numeric_limits<size_t>::max());
		}
		assert(overlap < nodePointers.back()->size());
		nodeStartOverlap.push_back(overlap);
		nodeStartPoses.push_back(realPos);
		realPos += nodePointers.back()->size() - overlap;
	}
	this->leftClip = leftClip;
	this->rightClip = rightClip;
}

char PathSequenceView::operator[](size_t index) const
{
	if (leftClip > 0) index += leftClip;
	size_t nodeIndex = 0;
	while (nodeIndex+1 < nodeStartPoses.size() && nodeStartPoses[nodeIndex+1] <= index) nodeIndex += 1;
	size_t nodeOffset = index - nodeStartPoses[nodeIndex] + nodeStartOverlap[nodeIndex];
	return (*(nodePointers[nodeIndex])).getChar(nodeOffset);
}

size_t PathSequenceView::size() const
{
	return nodeStartPoses.back() + nodePointers.back()->size() - nodeStartOverlap.back() - leftClip - rightClip;
}

size_t getOverlap(const Node& from, const Node& to, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges)
{
	size_t overlap = std::numeric_limits<size_t>::max();
	for (auto edge : edges.at(from))
	{
		if (std::get<0>(edge) != to) continue;
		assert(overlap == std::numeric_limits<size_t>::max());
		overlap = std::get<1>(edge);
	}
	assert(overlap != std::numeric_limits<size_t>::max());
	return overlap;
}

size_t getPathLength(const std::vector<Node>& nodes, const std::vector<std::string>& nodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges)
{
	size_t result = 0;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		size_t add = nodeSeqs.at(nodes[i].id()).size();
		if (i > 0)
		{
			size_t overlap = getOverlap(nodes[i-1], nodes[i], edges);
			add -= overlap;
		}
		result += add;
	}
	return result;
}

PathSequenceView getSequenceView(const std::vector<Node>& nodes, const std::vector<TwobitString>& nodeSeqs, const std::vector<TwobitString>& revCompNodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges, size_t leftClip, size_t rightClip)
{
	return PathSequenceView(nodes, nodeSeqs, revCompNodeSeqs, edges, leftClip, rightClip);
}

Path orientPath(const GfaGraph& graph, const Path& rawPath, const std::string& orientReferencePath, size_t k)
{
	Path result = rawPath;
	std::unordered_map<uint64_t, size_t> kmerReferencePosition;
	FastQ::streamFastqFromFile(orientReferencePath, false, [&kmerReferencePosition, k](FastQ& fastq)
	{
		for (size_t i = 0; i < fastq.sequence.size()-k; i++)
		{
			uint64_t hash = std::hash<std::string>{}(fastq.sequence.substr(i, k));
			if (kmerReferencePosition.count(hash) == 1)
			{
				kmerReferencePosition[hash] = std::numeric_limits<size_t>::max();
			}
			else
			{
				kmerReferencePosition[hash] = i;
			}
		}
	});
	{
		std::string pathSeq = result.getSequence(graph.nodeSeqs);
		std::string reverseSeq = revcomp(pathSeq);
		size_t fwMatches = 0;
		size_t bwMatches = 0;
		for (size_t i = 0; i < pathSeq.size() - k; i++)
		{
			uint64_t hash = std::hash<std::string>{}(pathSeq.substr(i, k));
			if (kmerReferencePosition.count(hash) == 1) fwMatches += 1; // no need to check if the hash is unique, it's still valid for orientation
		}
		for (size_t i = 0; i < reverseSeq.size() - k; i++)
		{
			uint64_t hash = std::hash<std::string>{}(reverseSeq.substr(i, k));
			if (kmerReferencePosition.count(hash) == 1) bwMatches += 1; // no need to check if the hash is unique, it's still valid for orientation
		}
		if (fwMatches == 0 && bwMatches == 0)
		{
			Logger::Log.log(Logger::LogLevel::DebugInfo) << "Can't be matched to the reference!" << std::endl;
			return rawPath;
		}
		if (bwMatches > fwMatches)
		{
			Logger::Log.log(Logger::LogLevel::DebugInfo) << "reverse complement" << std::endl;
			std::reverse(result.nodes.begin(), result.nodes.end());
			for (size_t i = 0; i < result.nodes.size(); i++)
			{
				result.nodes[i] = reverse(result.nodes[i]);
			}
			std::reverse(result.overlaps.begin(), result.overlaps.end());
			result.overlaps.insert(result.overlaps.begin(), result.overlaps.back());
			result.overlaps.pop_back();
		}
	}
	std::string pathSeq = result.getSequence(graph.nodeSeqs);
	// todo: this does not work in general due to spurious k-mer hits, but seems to be good enough for now?
	// should instead do chaining between matches and pick leftmost ref pos from them
	std::pair<size_t, size_t> minMatchPos { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
	for (size_t i = 0; i < pathSeq.size() - k; i++)
	{
		uint64_t hash = std::hash<std::string>{}(pathSeq.substr(i, k));
		if (kmerReferencePosition.count(hash) == 0) continue;
		auto refPos = kmerReferencePosition.at(hash);
		if (refPos == std::numeric_limits<size_t>::max()) continue;
		if (refPos < minMatchPos.second)
		{
			minMatchPos.first = i;
			minMatchPos.second = refPos;
		}
	}
	assert(minMatchPos.first != std::numeric_limits<size_t>::max());
	size_t rotatePosition = 0;
	if (minMatchPos.first >= minMatchPos.second)
	{
		rotatePosition = minMatchPos.first - minMatchPos.second;
	}
	else
	{
		assert(pathSeq.size() + minMatchPos.first >= minMatchPos.second);
		rotatePosition = pathSeq.size() + minMatchPos.first - minMatchPos.second;
	}
	size_t rotateNodes = 0;
	size_t pos = 0;
	size_t extraRotate = 0;
	for (size_t i = 1; i < result.nodes.size(); i++)
	{
		pos += graph.nodeSeqs.at(result.nodes[i-1].id()).size();
		pos -= result.overlaps[i];
		if (pos > rotatePosition)
		{
			break;
		}
		rotateNodes = i;
		assert(rotatePosition >= pos);
		extraRotate = rotatePosition - pos;
	}
	assert(rotateNodes != std::numeric_limits<size_t>::max());
	assert(rotateNodes < result.nodes.size());
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "rotate by " << rotateNodes << " nodes and " << extraRotate << " base pairs" << std::endl;
	if (rotateNodes != 0)
	{
		Path result2;
		result2.nodes.insert(result2.nodes.end(), result.nodes.begin() + rotateNodes, result.nodes.end());
		result2.nodes.insert(result2.nodes.end(), result.nodes.begin(), result.nodes.begin() + rotateNodes);
		result2.overlaps.insert(result2.overlaps.end(), result.overlaps.begin() + rotateNodes, result.overlaps.end());
		result2.overlaps.insert(result2.overlaps.end(), result.overlaps.begin(), result.overlaps.begin() + rotateNodes);
		result = result2;
	}
	result.leftClip = extraRotate;
	if (result.overlaps[0] <= result.leftClip)
	{
		result.nodes.push_back(result.nodes[0]);
		result.overlaps.push_back(result.overlaps[0]);
		result.rightClip = graph.nodeSeqs.at(result.nodes.back().id()).size() - extraRotate;
	}
	else
	{
		result.rightClip = result.overlaps[0] - result.leftClip;
	}
	return result;
}

std::unordered_set<size_t> getCoreNodes(const std::vector<OntLoop>& paths)
{
	std::unordered_set<size_t> existingNodes;
	std::unordered_set<size_t> notUniqueNodes;
	for (const auto& read : paths)
	{
		std::unordered_map<size_t, size_t> coverage;
		for (const auto& node : read.path)
		{
			coverage[node.id()] += 1;
		}
		for (auto pair : coverage)
		{
			if (pair.second == 1) existingNodes.insert(pair.first);
			if (pair.second != 1) notUniqueNodes.insert(pair.first);
		}
	}
	for (const auto& read : paths)
	{
		std::unordered_set<size_t> nodesHere;
		for (const auto& node : read.path) nodesHere.insert(node.id());
		for (auto node : existingNodes)
		{
			if (nodesHere.count(node) == 0) notUniqueNodes.insert(node);
		}
	}
	std::unordered_set<size_t> firstNodes;
	std::unordered_set<size_t> lastNodes;
	for (const auto& read : paths)
	{
		firstNodes.insert(read.path[0].id());
		lastNodes.insert(read.path.back().id());
	}
	std::unordered_set<size_t> potentialCoreNodes;
	for (auto node : existingNodes)
	{
		if (notUniqueNodes.count(node) == 1) continue;
		if (firstNodes.count(node) == 1 && lastNodes.count(node) == 1) continue;
		potentialCoreNodes.insert(node);
	}
	// static_assert(false, "todo: this did nothing before, check carefully");
	// std::unordered_map<size_t, size_t> uniqueCorePredecessor;
	// std::unordered_map<size_t, size_t> uniqueCoreSuccessor;
	// for (const auto& read : paths)
	// {
	// 	size_t lastCore = std::numeric_limits<size_t>::max();
	// 	for (const auto& node : read.path)
	// 	{
	// 		size_t nodename = node.id();
	// 		if (potentialCoreNodes.count(nodename) == 0) continue;
	// 		if (lastCore.size() == 0) continue;
	// 		if (uniqueCorePredecessor.count(nodename) == 1 && uniqueCorePredecessor.at(nodename) != lastCore) uniqueCorePredecessor[nodename] = std::numeric_limits<size_t>::max();
	// 		if (uniqueCoreSuccessor.count(lastCore) == 1 && uniqueCoreSuccessor.at(lastCore) != nodename) uniqueCoreSuccessor[lastCore] = std::numeric_limits<size_t>::max();
	// 		if (uniqueCorePredecessor.count(nodename) == 0) uniqueCorePredecessor[nodename] = lastCore;
	// 		if (uniqueCoreSuccessor.count(lastCore) == 0) uniqueCoreSuccessor[lastCore] = nodename;
	// 		lastCore = nodename;
	// 	}
	// }
	// for (const auto& pair : uniqueCoreSuccessor)
	// {
	// 	if (pair.second != std::numeric_limits<size_t>::max()) continue;
	// 	if (uniqueCorePredecessor.count(pair.first) == 0) continue;
	// 	if (uniqueCorePredecessor.at(pair.first) != std::numeric_limits<size_t>::max()) continue;
	// 	potentialCoreNodes.erase(pair.first);
	// }
	std::unordered_map<size_t, size_t> uniqueCoreIndex;
	for (const auto& read : paths)
	{
		size_t coreIndex = 0;
		for (const auto& node : read.path)
		{
			size_t nodename = node.id();
			if (potentialCoreNodes.count(nodename) == 0) continue;
			if (uniqueCoreIndex.count(nodename) == 1 && uniqueCoreIndex.at(nodename) != coreIndex) uniqueCoreIndex[nodename] = std::numeric_limits<size_t>::max();
			if (uniqueCoreIndex.count(nodename) == 0) uniqueCoreIndex[nodename] = coreIndex;
			coreIndex += 1;
		}
	}
	for (const auto& pair : uniqueCoreIndex)
	{
		if (pair.second != std::numeric_limits<size_t>::max()) continue;
		potentialCoreNodes.erase(pair.first);
	}
	return potentialCoreNodes;
}

bool orientPath(std::vector<Node>& path, const std::unordered_map<size_t, bool>& referenceOrientations)
{
	size_t fwMatches = 0;
	size_t bwMatches = 0;
	for (const auto& node : path)
	{
		if (referenceOrientations.count(node.id()) == 0) continue;
		bool orientHere = node.forward();
		bool refOrient = referenceOrientations.at(node.id());
		if (orientHere == refOrient)
		{
			fwMatches += 1;
		}
		else
		{
			bwMatches += 1;
		}
	}
	if (bwMatches > fwMatches)
	{
		path = reverse(path);
		return true;
	}
	return false;
}

phmap::flat_hash_set<size_t> getNodesInMainCycle(const GfaGraph& graph, const Path& heavyPath, const size_t minCoverage)
{
	std::vector<Node> checkStack;
	phmap::flat_hash_set<size_t> coveredNodes;
	phmap::flat_hash_set<Node> reachableNodes;
	for (size_t i = 0; i < heavyPath.nodes.size(); i++)
	{
		coveredNodes.insert(heavyPath.nodes[i].id());
		checkStack.emplace_back(Node { heavyPath.nodes[i].id(), true });
		checkStack.emplace_back(Node { heavyPath.nodes[i].id(), false });
	}
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (graph.nodeCoverages[i] < minCoverage) continue;
		coveredNodes.insert(i);
	}
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (reachableNodes.count(top) == 1) continue;
		reachableNodes.insert(top);
		if (graph.edges.count(top) == 1)
		{
			for (auto edge : graph.edges.at(top))
			{
				auto target = std::get<0>(edge);
				if (coveredNodes.count(target.id()) == 0) continue;
				checkStack.emplace_back(target);
			}
		}
	}
	phmap::flat_hash_set<size_t> nodesInMainCycle;
	for (size_t node : coveredNodes)
	{
		if (reachableNodes.count(Node { node, true }) == 0) continue;
		if (reachableNodes.count(Node { node, false }) == 0) continue;
		nodesInMainCycle.insert(node);
	}
	return nodesInMainCycle;
}

phmap::flat_hash_map<Node, size_t> getNodeDistancesToMainCycle(const GfaGraph& graph, const phmap::flat_hash_set<size_t>& keptNodesInCycles)
{
	phmap::flat_hash_map<Node, size_t> shortestDistanceToCyclicComponent;
	std::vector<std::pair<Node, size_t>> checkStack;
	for (size_t node : keptNodesInCycles)
	{
		checkStack.emplace_back(Node { node, true }, 0);
		checkStack.emplace_back(Node { node, false }, 0);
	}
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (shortestDistanceToCyclicComponent.count(top.first) == 1)
		{
			assert(shortestDistanceToCyclicComponent.at(top.first) <= top.second);
			continue;
		}
		shortestDistanceToCyclicComponent[top.first] = top.second;
		bool addedAny = false;
		if (graph.edges.count(top.first) == 1)
		{
			for (auto edge : graph.edges.at(top.first))
			{
				if (keptNodesInCycles.count(std::get<0>(edge).id()) == 1) continue;
				checkStack.emplace_back(std::get<0>(edge), top.second + graph.nodeSeqs[std::get<0>(edge).id()].size() - std::get<1>(edge));
				assert(checkStack.back().second > top.second);
				addedAny = true;
			}
		}
		if (addedAny) std::sort(checkStack.begin(), checkStack.end(), [](auto left, auto right) { return left.second > right.second; });
	}
	return shortestDistanceToCyclicComponent;
}
