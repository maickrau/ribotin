#include <map>
#include <limits>
#include "phmap.h"
#include "HeavyPath.h"
#include "TwobitString.h"
#include "Logger.h"

class PathRecWideCoverage
{
public:
	PathRecWideCoverage() = default;
	PathRecWideCoverage(const size_t coverage, const size_t length)
	{
		coverages.emplace_back(coverage, length);
	}
	bool operator<(const PathRecWideCoverage& other) const
	{
		for (size_t i = 0; i < coverages.size() && i < other.coverages.size(); i++)
		{
			if (coverages[i].first < other.coverages[i].first) return true;
			if (coverages[i].first > other.coverages[i].first) return false;
			if (coverages[i].second < other.coverages[i].second) return false;
			if (coverages[i].second > other.coverages[i].second) return true;
		}
		if (coverages.size() < other.coverages.size()) return true;
		return false;
	}
	void addCoverage(const size_t coverage, const size_t length)
	{
		coverages.emplace_back(coverage, length);
		std::sort(coverages.begin(), coverages.end());
		bool removedAny = false;
		for (size_t i = coverages.size()-1; i > 0; i--)
		{
			if (coverages[i-1].first != coverages[i].first) continue;
			coverages[i-1].second += coverages[i].second;
			std::swap(coverages[i], coverages.back());
			coverages.pop_back();
			removedAny = true;
		}
		if (removedAny)
		{
			std::sort(coverages.begin(), coverages.end());
		}
	}
	PathRecWideCoverage& operator+=(const PathRecWideCoverage& other)
	{
		for (size_t i = 0; i < other.coverages.size(); i++)
		{
			addCoverage(other.coverages[i].first, other.coverages[i].second);
		}
		return *this;
	}
	std::vector<std::pair<size_t, size_t>> coverages;
private:
};

bool localComponentIsAcyclic(const std::vector<std::vector<size_t>>& nodesInLocalComponent, const std::vector<size_t>& selfDistances, const size_t estimatedMinimalSelfDistance, const size_t component)
{
	if (nodesInLocalComponent.at(component).size() == 1)
	{
		if (selfDistances[nodesInLocalComponent.at(component)[0]] >= estimatedMinimalSelfDistance)
		{
			return true;
		}
	}
	return false;
}

void addAllShortestPaths(std::vector<std::vector<std::tuple<PathRecWideCoverage, std::vector<Node>, size_t>>>& componentInEdges, const std::vector<bool>& nodeExists, const GfaGraph& graph, const std::vector<bool>& nodeOrientation, const std::vector<size_t>& nodesInComponent, const std::vector<size_t>& nodeBelongsToStronglyConnectedLocalComponent)
{
	assert(nodesInComponent.size() >= 1);
	phmap::flat_hash_set<Node> nodesBefore;
	phmap::flat_hash_set<Node> nodesAfter;
	size_t thisComponent = nodeBelongsToStronglyConnectedLocalComponent[nodesInComponent[0]];
	for (size_t unorientedNode : nodesInComponent)
	{
		Node node { unorientedNode, nodeOrientation[unorientedNode] };
		assert(nodeBelongsToStronglyConnectedLocalComponent[node.id()] == thisComponent);
		assert(graph.edges.count(node) == 1);
		assert(graph.edges.count(reverse(node)) == 1);
		for (auto edge : graph.edges.at(node))
		{
			Node target = std::get<0>(edge);
			if (!nodeExists[target.id()]) continue;
			if (nodeBelongsToStronglyConnectedLocalComponent[node.id()] == nodeBelongsToStronglyConnectedLocalComponent[target.id()]) continue;
			assert(nodeBelongsToStronglyConnectedLocalComponent[target.id()] > thisComponent || nodeBelongsToStronglyConnectedLocalComponent[target.id()] == 0);
			nodesAfter.insert(target);
		}
		for (auto edge : graph.edges.at(reverse(node)))
		{
			Node target = reverse(std::get<0>(edge));
			if (!nodeExists[target.id()]) continue;
			if (nodeBelongsToStronglyConnectedLocalComponent[node.id()] == nodeBelongsToStronglyConnectedLocalComponent[target.id()]) continue;
			if (!(nodeBelongsToStronglyConnectedLocalComponent[target.id()] < thisComponent))
			{
				std::cerr << graph.nodeNames[unorientedNode] << " " << graph.nodeNames[target.id()] << std::endl;
				std::cerr << nodeBelongsToStronglyConnectedLocalComponent[target.id()] << " " << thisComponent << std::endl;
			}
			assert(nodeBelongsToStronglyConnectedLocalComponent[target.id()] < thisComponent);
			nodesBefore.insert(target);
		}
	}
	assert(nodesBefore.size() >= 1);
	assert(nodesAfter.size() >= 1);
	for (Node before : nodesBefore)
	{
		std::vector<std::pair<size_t, std::pair<Node, Node>>> checkStack;
		phmap::flat_hash_map<Node, Node> predecessor;
		assert(graph.edges.count(before) == 1);
		for (auto edge : graph.edges.at(before))
		{
			Node target = std::get<0>(edge);
			if (nodeBelongsToStronglyConnectedLocalComponent[target.id()] != thisComponent) continue;
			size_t overlap = std::get<1>(edge);
			checkStack.emplace_back(graph.nodeSeqs[target.id()].size() - overlap, std::make_pair(before, target));
		}
		std::sort(checkStack.begin(), checkStack.end(), [](auto left, auto right) { return left.first > right.first; });
		assert(checkStack.size() >= 1);
		while (checkStack.size() >= 1)
		{
			auto top = checkStack.back();
			assert(nodeBelongsToStronglyConnectedLocalComponent[top.second.second.id()] == thisComponent);
			checkStack.pop_back();
			if (predecessor.count(top.second.second) == 1) continue;
			predecessor[top.second.second] = top.second.first;
			assert(graph.edges.count(top.second.second) == 1);
			for (auto edge : graph.edges.at(top.second.second))
			{
				Node target = std::get<0>(edge);
				if (!nodeExists[target.id()]) continue;
				if (nodeBelongsToStronglyConnectedLocalComponent[target.id()] != thisComponent)
				{
					assert(nodesAfter.count(target) == 1);
					std::vector<Node> pathHere;
					pathHere.emplace_back(target);
					pathHere.emplace_back(top.second.second);
					while (predecessor.count(pathHere.back()) == 1) pathHere.emplace_back(predecessor.at(pathHere.back()));
					std::reverse(pathHere.begin(), pathHere.end());
					size_t pathLength = getPathLength(pathHere, graph.nodeSeqs, graph.edges);
					size_t sourceComponent = nodeBelongsToStronglyConnectedLocalComponent[before.id()];
					size_t targetComponent = nodeBelongsToStronglyConnectedLocalComponent[target.id()];
					assert(targetComponent < componentInEdges.size());
					assert(sourceComponent < targetComponent || targetComponent == 0);
					componentInEdges[targetComponent].emplace_back(PathRecWideCoverage { 1, pathLength }, pathHere, sourceComponent);
					continue;
				}
				size_t overlap = std::get<1>(edge);
				checkStack.emplace_back(graph.nodeSeqs[target.id()].size() - overlap, std::make_pair(top.second.second, target));
			}
		}
	}
}

// components in nodesInLocalComponent and nodeBelongsToStronglyConnectedLocalComponent must be sorted by topological order
Path getTopologicallyOrderedHeaviestPath(const GfaGraph& graph, const std::vector<ReadPath>& readPaths, const std::vector<bool>& nodeExists, const std::vector<size_t>& selfDistances, const size_t estimatedMinimalSelfDistance, const std::vector<bool>& nodeOrientation, const std::vector<size_t>& nodeBelongsToStronglyConnectedLocalComponent, const std::vector<std::vector<size_t>>& nodesInLocalComponent)
{
	assert(nodesInLocalComponent[0].size() == 1);
	std::vector<bool> componentIsCyclic;
	componentIsCyclic.resize(nodesInLocalComponent.size(), false);
	for (size_t i = 0; i < nodesInLocalComponent.size(); i++)
	{
		if (localComponentIsAcyclic(nodesInLocalComponent, selfDistances, estimatedMinimalSelfDistance, i)) continue;
		componentIsCyclic[i] = true;
	}
	assert(!componentIsCyclic[0]);
	std::vector<std::vector<std::tuple<PathRecWideCoverage, std::vector<Node>, size_t>>> componentInEdges;
	componentInEdges.resize(nodesInLocalComponent.size());
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (!nodeExists[i]) continue;
		size_t thisComponent = nodeBelongsToStronglyConnectedLocalComponent[i];
		if (componentIsCyclic[thisComponent]) continue;
		Node nodeHere { i, nodeOrientation[i] };
		assert(graph.edges.count(nodeHere) == 1);
		assert(graph.edges.count(reverse(nodeHere)) == 1);
		for (auto edge : graph.edges.at(nodeHere))
		{
			Node target = std::get<0>(edge);
			if (!nodeExists[target.id()]) continue;
			size_t targetComponent = nodeBelongsToStronglyConnectedLocalComponent.at(target.id());
			if (componentIsCyclic[targetComponent]) continue;
			size_t overlap = std::get<1>(edge);
			size_t coverage = std::get<2>(edge);
			assert(thisComponent < targetComponent || targetComponent == 0);
			componentInEdges[targetComponent].emplace_back(PathRecWideCoverage { coverage, 1 }, std::vector<Node> { nodeHere, target }, thisComponent);
			std::get<0>(componentInEdges[targetComponent].back()).addCoverage(graph.nodeCoverages[nodeHere.id()], graph.nodeSeqs[nodeHere.id()].size() - overlap);
		}
	}
	std::map<std::vector<Node>, size_t> cyclicComponentThroughGoerCoverages;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		size_t lastAcyclicIndex = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (!nodeExists[readPaths[i].path[j].id()])
			{
				lastAcyclicIndex = std::numeric_limits<size_t>::max();
				continue;
			}
			if (componentIsCyclic[nodeBelongsToStronglyConnectedLocalComponent[readPaths[i].path[j].id()]]) continue;
			if (lastAcyclicIndex == std::numeric_limits<size_t>::max())
			{
				lastAcyclicIndex = j;
				continue;
			}
			if (j == lastAcyclicIndex+1)
			{
				lastAcyclicIndex = j;
				continue;
			}
			std::vector<Node> path { readPaths[i].path.begin()+lastAcyclicIndex, readPaths[i].path.begin() + j + 1};
			assert(path.size() >= 3);
			assert(nodeExists[path[0].id()]);
			assert(nodeExists[path.back().id()]);
			assert(nodeBelongsToStronglyConnectedLocalComponent[path.back().id()] != nodeBelongsToStronglyConnectedLocalComponent[path[0].id()]);
			if (nodeBelongsToStronglyConnectedLocalComponent[path.back().id()] < nodeBelongsToStronglyConnectedLocalComponent[path[0].id()])
			{
				path = reverse(path);
			}
			cyclicComponentThroughGoerCoverages[path] += 1;
			lastAcyclicIndex = j;
		}
	}
	for (const auto& pair : cyclicComponentThroughGoerCoverages)
	{
		assert(pair.first.size() >= 3);
		for (size_t i = 0; i < pair.first.size(); i++)
		{
			assert(nodeExists[pair.first[i].id()]);
		}
		size_t preComponent = nodeBelongsToStronglyConnectedLocalComponent[pair.first[0].id()];
		size_t postComponent = nodeBelongsToStronglyConnectedLocalComponent[pair.first.back().id()];
		assert(!componentIsCyclic[preComponent]);
		assert(!componentIsCyclic[postComponent]);
		for (size_t i = 1; i+1 < pair.first.size(); i++)
		{
			assert(componentIsCyclic[nodeBelongsToStronglyConnectedLocalComponent[pair.first[i].id()]]);
		}
		size_t pathLength = getPathLength(pair.first, graph.nodeSeqs, graph.edges);
		assert(preComponent < postComponent || postComponent == 0);
		componentInEdges[postComponent].emplace_back(PathRecWideCoverage { pair.second, pathLength }, pair.first, preComponent);
	}
	for (size_t i = 0; i < nodesInLocalComponent.size(); i++)
	{
		if (!componentIsCyclic[i]) continue;
		addAllShortestPaths(componentInEdges, nodeExists, graph, nodeOrientation, nodesInLocalComponent[i], nodeBelongsToStronglyConnectedLocalComponent);
	}
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (!nodeExists[i]) continue;
	}
	for (size_t i = 0; i < componentInEdges.size(); i++)
	{
		for (auto edge : componentInEdges[i])
		{
			if (!(std::get<2>(edge) < i || i == 0))
			{
				std::cerr << std::get<2>(edge) << " " << i << std::endl;
			}
			assert(std::get<2>(edge) < i || i == 0);
		}
	}
	std::vector<std::tuple<PathRecWideCoverage, std::vector<Node>, size_t>> predecessor;
	predecessor.resize(nodesInLocalComponent.size());
	for (size_t i = 0; i < nodesInLocalComponent.size(); i++)
	{
		if (componentIsCyclic[i]) continue;
		assert(componentInEdges[i].size() >= 1);
		std::tuple<PathRecWideCoverage, std::vector<Node>, size_t> bestSoFar = componentInEdges[i][0];
		std::get<0>(bestSoFar) += std::get<0>(predecessor[std::get<2>(bestSoFar)]);
		for (size_t j = 1; j < componentInEdges[i].size(); j++)
		{
			auto optionHere = componentInEdges[i][j];
			std::get<0>(optionHere) += std::get<0>(predecessor[std::get<2>(optionHere)]);
			if (std::get<0>(bestSoFar) < std::get<0>(optionHere)) bestSoFar = optionHere;
		}
		predecessor[i] = bestSoFar;
	}
	Path result;
	size_t pos = nodesInLocalComponent.size()-1;
	while (true)
	{
		std::vector<Node> addHere = reverse(std::get<1>(predecessor[pos]));
		assert(result.nodes.size() == 0 || addHere[0] == result.nodes.back());
		result.nodes.insert(result.nodes.end(), addHere.begin()+1, addHere.end());
		pos = std::get<2>(predecessor[pos]);
		if (pos == nodesInLocalComponent.size()-1) break;
		assert(result.nodes.size() < nodeExists.size());
	}
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
		}
		assert(graph.edges.count(prev) == 1);
		for (auto edge : graph.edges.at(prev))
		{
			if (std::get<0>(edge) != result.nodes[i]) continue;
			result.overlaps.push_back(std::get<1>(edge));
		}
	}
	assert(result.overlaps.size() == result.nodes.size());
	result.leftClip = 0;
	result.rightClip = result.overlaps[0];
	return result;
}

std::vector<size_t> getLocalComponentTopologicalOrder(const std::vector<phmap::flat_hash_set<size_t>>& localComponentInEdges, const std::vector<phmap::flat_hash_set<size_t>>& localComponentOutEdges, const size_t startComponent)
{
	std::vector<bool> visited;
	visited.resize(localComponentInEdges.size(), false);
	visited[startComponent] = true;
	std::vector<size_t> checkStack;
	assert(localComponentInEdges[startComponent].size() >= 1);
	assert(localComponentOutEdges[startComponent].size() >= 1);
	for (size_t edge : localComponentOutEdges[startComponent])
	{
		checkStack.emplace_back(edge);
	}
	std::vector<size_t> result;
	std::cerr << "start component " << startComponent << std::endl;
	while (checkStack.size() > 0)
	{
		size_t top = checkStack.back();
		checkStack.pop_back();
		assert(localComponentInEdges[top].size() >= 1);
		if (visited[top]) continue;
		bool allPredecessorsVisited = true;
		for (size_t neighbor : localComponentInEdges[top])
		{
			if (!visited[neighbor]) allPredecessorsVisited = false;
		}
		if (!allPredecessorsVisited) continue;
		visited[top] = true;
		result.emplace_back(top);
		std::cerr << "visit " << top << std::endl;
		for (size_t neighbor : localComponentOutEdges[top])
		{
			checkStack.emplace_back(neighbor);
		}
	}
	std::cerr << result.size() << std::endl;
	std::cerr << localComponentInEdges.size() << std::endl;
	assert(result.size()+1 == localComponentInEdges.size());
	result.emplace_back(startComponent);
	assert(result.size() == visited.size());
	assert(result.size() == localComponentInEdges.size());
	for (size_t i = 0; i < visited.size(); i++)
	{
		assert(visited[i]);
		assert(i+1 == visited.size() || result[i] != startComponent);
	}
	return result;
}

std::vector<bool> getNodeOrientations(const GfaGraph& graph, const size_t startNode, const std::vector<bool>& nodeExists)
{
	assert(nodeExists[startNode]);
	std::vector<bool> result;
	result.resize(graph.numNodes(), true);
	std::vector<bool> visited;
	visited.resize(graph.numNodes(), false);
	std::vector<Node> checkStack;
	checkStack.emplace_back(Node { startNode, true });
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		assert(nodeExists[top.id()]);
		checkStack.pop_back();
		if (visited[top.id()])
		{
			assert(result[top.id()] == top.forward());
			continue;
		}
		visited[top.id()] = true;
		result[top.id()] = top.forward();
		if (graph.edges.count(top) == 1)
		{
			for (auto edge : graph.edges.at(top))
			{
				if (std::get<0>(edge).id() == top.id()) continue;
				if (!nodeExists[std::get<0>(edge).id()]) continue;
				checkStack.emplace_back(std::get<0>(edge));
			}
		}
		if (graph.edges.count(reverse(top)) == 1)
		{
			for (auto revedge : graph.edges.at(reverse(top)))
			{
				if (std::get<0>(revedge).id() == top.id()) continue;
				if (!nodeExists[std::get<0>(revedge).id()]) continue;
				checkStack.emplace_back(reverse(std::get<0>(revedge)));
			}
		}
	}
	return result;
}

size_t getShortestPathDistance(const GfaGraph& graph, const Node startNode, const Node endNode)
{
	std::vector<std::pair<size_t, Node>> checkStack;
	if (graph.edges.count(startNode) == 0) return std::numeric_limits<size_t>::max();
	if (graph.edges.count(reverse(endNode)) == 0) return std::numeric_limits<size_t>::max();
	for (auto edge : graph.edges.at(startNode))
	{
		Node target = std::get<0>(edge);
		size_t overlap = std::get<1>(edge);
		checkStack.emplace_back(graph.nodeSeqs[target.id()].size() - overlap, target);
	}
	std::sort(checkStack.begin(), checkStack.end(), [](auto left, auto right) { return left.first > right.first; });
	phmap::flat_hash_set<Node> visited;
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (visited.count(top.second) == 1) continue;
		visited.insert(top.second);
		if (top.second == endNode)
		{
			return top.first;
		}
		if (graph.edges.count(top.second) == 0) continue;
		for (auto edge : graph.edges.at(top.second))
		{
			Node target = std::get<0>(edge);
			size_t overlap = std::get<1>(edge);
			checkStack.emplace_back(top.first + graph.nodeSeqs[target.id()].size() - overlap, target);
		}
		std::sort(checkStack.begin(), checkStack.end(), [](auto left, auto right) { return left.first > right.first; });
	}
	return std::numeric_limits<size_t>::max();
}

size_t getSelfDistance(const GfaGraph& graph, const size_t startNode)
{
	return getShortestPathDistance(graph, Node { startNode, true }, Node { startNode, true });
}

// https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
void strongConnectRecurse(std::vector<size_t>& result, std::vector<size_t>& index, std::vector<size_t>& lowLink, std::vector<bool>& onStack, std::vector<size_t>& stack, size_t& nextIndex, size_t& nextStronglyConnectedComponent, const GfaGraph& graph, const std::vector<bool>& nodeOrientation, const std::vector<bool>& nodeExists, const size_t v, const size_t globalStartNode)
{
	assert(nodeExists[v]);
	assert(!onStack[v]);
	index[v] = nextIndex;
	lowLink[v] = nextIndex;
	nextIndex += 1;
	stack.push_back(v);
	onStack[v] = true;
	for (auto edge : graph.edges.at(Node { v, nodeOrientation[v] }))
	{
		size_t w = std::get<0>(edge).id();
		if (w == globalStartNode) continue;
		if (!nodeExists[w]) continue;
		if (index[w] == std::numeric_limits<size_t>::max())
		{
			strongConnectRecurse(result, index, lowLink, onStack, stack, nextIndex, nextStronglyConnectedComponent, graph, nodeOrientation, nodeExists, w, globalStartNode);
			lowLink[v] = std::min(lowLink[v], lowLink[w]);
		}
		else if (onStack[w])
		{
			lowLink[v] = std::min(lowLink[v], index[w]);
		}
	}
	if (lowLink[v] == index[v])
	{
		while (true)
		{
			assert(stack.size() >= 1);
			size_t w = stack.back();
			onStack[w] = false;
			stack.pop_back();
			result[w] = nextStronglyConnectedComponent;
			if (w == v) break;
		}
		nextStronglyConnectedComponent += 1;
	}
}

std::vector<size_t> getStronglyConnectedLocalComponents(const GfaGraph& graph, const std::vector<bool>& nodeOrientation, const std::vector<bool>& nodeExists, const size_t bestStartNode)
{
	std::vector<size_t> result;
	result.resize(graph.numNodes(), std::numeric_limits<size_t>::max());
	std::vector<size_t> index;
	index.resize(graph.numNodes(), std::numeric_limits<size_t>::max());
	std::vector<size_t> lowLink;
	lowLink.resize(graph.numNodes(), std::numeric_limits<size_t>::max());
	std::vector<bool> onStack;
	onStack.resize(graph.numNodes(), false);
	std::vector<size_t> stack;
	size_t nextIndex = 0;
	size_t nextStronglyConnectedComponent = 0;
	strongConnectRecurse(result, index, lowLink, onStack, stack, nextIndex, nextStronglyConnectedComponent, graph, nodeOrientation, nodeExists, bestStartNode, bestStartNode);
	assert(stack.size() == 0);
	return result;
}

Path getHeavyPath(const GfaGraph& graph, const std::vector<ReadPath>& readPaths, const size_t localResolveLength)
{
	std::vector<size_t> selfDistances;
	size_t distanceSum = 0;
	size_t distanceDivisor = 0;
	std::vector<std::pair<size_t, size_t>> selfDistanceAndLengths;
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		size_t selfDistance = getSelfDistance(graph, i);
		selfDistances.emplace_back(selfDistance);
		if (selfDistance == std::numeric_limits<size_t>::max()) continue;
		if (selfDistance < localResolveLength) continue;
		selfDistanceAndLengths.emplace_back(selfDistance, graph.nodeSeqs[i].size());
		distanceSum += selfDistance * graph.nodeSeqs[i].size();
		distanceDivisor += graph.nodeSeqs[i].size();
	}
	std::sort(selfDistanceAndLengths.begin(), selfDistanceAndLengths.end());
	size_t totalLength = 0;
	for (auto pair : selfDistanceAndLengths)
	{
		totalLength += pair.second;
	}
	assert(totalLength > 0);
	size_t estimatedMinimalSelfDistance = std::numeric_limits<size_t>::max();
	size_t sum = 0;
	for (auto pair : selfDistanceAndLengths)
	{
		sum += pair.second;
		if (sum*2 >= totalLength)
		{
			estimatedMinimalSelfDistance = pair.first/2;
			break;
		}
	}
	assert(estimatedMinimalSelfDistance != std::numeric_limits<size_t>::max());
	Logger::Log.log(Logger::LogLevel::Always) << "using " << estimatedMinimalSelfDistance << " as local self-repeat distance threshold" << std::endl;
	size_t bestStartNode = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (selfDistances[i] == std::numeric_limits<size_t>::max()) continue;
		if (selfDistances[i] < estimatedMinimalSelfDistance) continue;
		if (bestStartNode == std::numeric_limits<size_t>::max()) bestStartNode = i;
		if (graph.nodeCoverages[i] > graph.nodeCoverages[bestStartNode]) bestStartNode = i;
	}
	assert(bestStartNode != std::numeric_limits<size_t>::max());
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "best start node " << graph.nodeNames[bestStartNode] << std::endl;
	std::vector<bool> nodeExists;
	nodeExists.resize(graph.numNodes(), false);
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		nodeExists[i] = selfDistances[i] != std::numeric_limits<size_t>::max();
	}
	std::vector<bool> nodeOrientation = getNodeOrientations(graph, bestStartNode, nodeExists);
	nodeExists.resize(0);
	nodeExists.resize(graph.numNodes(), false);
	nodeExists[bestStartNode] = true;
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (i == bestStartNode) continue;
		if (selfDistances[i] == std::numeric_limits<size_t>::max()) continue;
		size_t distancePre = getShortestPathDistance(graph, Node { bestStartNode, true }, Node { i, nodeOrientation[i] });
		if (distancePre == std::numeric_limits<size_t>::max()) continue;
		size_t distancePost = getShortestPathDistance(graph, Node { i, nodeOrientation[i] }, Node { bestStartNode, true });
		if (distancePost == std::numeric_limits<size_t>::max()) continue;
		if (distancePre + distancePost > selfDistances[bestStartNode] + estimatedMinimalSelfDistance && distancePre + distancePost > selfDistances[i] + estimatedMinimalSelfDistance) continue;
		nodeExists[i] = true;
	}
	std::vector<size_t> nodeBelongsToStronglyConnectedLocalComponent = getStronglyConnectedLocalComponents(graph, nodeOrientation, nodeExists, bestStartNode);
	size_t countComponents = 0;
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (!nodeExists[i]) continue;
		assert(nodeBelongsToStronglyConnectedLocalComponent[i] != std::numeric_limits<size_t>::max());
		countComponents = std::max(countComponents, nodeBelongsToStronglyConnectedLocalComponent[i]);
	}
	assert(nodeBelongsToStronglyConnectedLocalComponent[bestStartNode] == countComponents);
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (!nodeExists[i]) continue;
		nodeBelongsToStronglyConnectedLocalComponent[i] = countComponents - nodeBelongsToStronglyConnectedLocalComponent[i];
	}
	countComponents += 1;
	std::vector<std::vector<size_t>> nodesInLocalComponent;
	nodesInLocalComponent.resize(countComponents);
	for (size_t i = 0; i < nodeBelongsToStronglyConnectedLocalComponent.size(); i++)
	{
		if (!nodeExists[i]) continue;
		nodesInLocalComponent[nodeBelongsToStronglyConnectedLocalComponent[i]].emplace_back(i);
	}
	Path result = getTopologicallyOrderedHeaviestPath(graph, readPaths, nodeExists, selfDistances, estimatedMinimalSelfDistance, nodeOrientation, nodeBelongsToStronglyConnectedLocalComponent, nodesInLocalComponent);
	return result;
}
