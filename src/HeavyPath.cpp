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
			nodesAfter.insert(target);
		}
		for (auto edge : graph.edges.at(reverse(node)))
		{
			Node target = reverse(std::get<0>(edge));
			if (!nodeExists[target.id()]) continue;
			if (nodeBelongsToStronglyConnectedLocalComponent[node.id()] == nodeBelongsToStronglyConnectedLocalComponent[target.id()]) continue;
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
					componentInEdges[targetComponent].emplace_back(PathRecWideCoverage { 1, pathLength }, pathHere, sourceComponent);
					continue;
				}
				size_t overlap = std::get<1>(edge);
				checkStack.emplace_back(graph.nodeSeqs[target.id()].size() - overlap, std::make_pair(top.second.second, target));
			}
		}
	}
}

Path getTopologicallyOrderedHeaviestPath(const GfaGraph& graph, const std::vector<ReadPath>& readPaths, const std::vector<bool>& nodeExists, const std::vector<size_t>& selfDistances, const size_t estimatedMinimalSelfDistance, const std::vector<bool>& nodeOrientation, const std::vector<size_t>& localComponentOrder, const std::vector<size_t>& nodeBelongsToStronglyConnectedLocalComponent, const std::vector<std::vector<size_t>>& nodesInLocalComponent)
{
	std::vector<size_t> componentOrderInTopologicalOrder;
	componentOrderInTopologicalOrder.resize(localComponentOrder.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < localComponentOrder.size(); i++)
	{
		assert(componentOrderInTopologicalOrder[localComponentOrder[i]] == std::numeric_limits<size_t>::max());
		componentOrderInTopologicalOrder[localComponentOrder[i]] = i;
	}
	for (size_t i = 0; i < componentOrderInTopologicalOrder.size(); i++)
	{
		assert(componentOrderInTopologicalOrder[i] != std::numeric_limits<size_t>::max());
	}
	assert(localComponentOrder.size() >= 1);
	assert(localComponentOrder.size() == nodesInLocalComponent.size());
	assert(nodesInLocalComponent[localComponentOrder[0]].size() == 1);
	std::vector<bool> componentIsCyclic;
	componentIsCyclic.resize(localComponentOrder.size(), false);
	for (size_t i = 0; i < localComponentOrder.size(); i++)
	{
		if (localComponentIsAcyclic(nodesInLocalComponent, selfDistances, estimatedMinimalSelfDistance, i)) continue;
		componentIsCyclic[i] = true;
	}
	assert(!componentIsCyclic[localComponentOrder[0]]);
	std::vector<std::vector<std::tuple<PathRecWideCoverage, std::vector<Node>, size_t>>> componentInEdges;
	componentInEdges.resize(localComponentOrder.size());
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
			assert(componentOrderInTopologicalOrder[thisComponent] < componentOrderInTopologicalOrder[targetComponent] || thisComponent == localComponentOrder.back());
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
			if (componentOrderInTopologicalOrder[nodeBelongsToStronglyConnectedLocalComponent[path.back().id()]] < componentOrderInTopologicalOrder[nodeBelongsToStronglyConnectedLocalComponent[path[0].id()]])
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
		assert(componentOrderInTopologicalOrder[preComponent] < componentOrderInTopologicalOrder[postComponent] || preComponent == localComponentOrder.back());
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
			assert(componentOrderInTopologicalOrder[std::get<2>(edge)] < componentOrderInTopologicalOrder[i] || std::get<2>(edge) == localComponentOrder.back());
		}
	}
	std::vector<std::tuple<PathRecWideCoverage, std::vector<Node>, size_t>> predecessor;
	predecessor.resize(localComponentOrder.size());
	for (size_t i = 0; i < localComponentOrder.size(); i++)
	{
		size_t component = localComponentOrder[i];
		if (componentIsCyclic[component]) continue;
		assert(componentInEdges[component].size() >= 1);
		std::tuple<PathRecWideCoverage, std::vector<Node>, size_t> bestSoFar = componentInEdges[component][0];
		std::get<0>(bestSoFar) += std::get<0>(predecessor[std::get<2>(bestSoFar)]);
		for (size_t j = 1; j < componentInEdges[component].size(); j++)
		{
			auto optionHere = componentInEdges[component][j];
			std::get<0>(optionHere) += std::get<0>(predecessor[std::get<2>(optionHere)]);
			if (std::get<0>(bestSoFar) < std::get<0>(optionHere)) bestSoFar = optionHere;
		}
		predecessor[component] = bestSoFar;
	}
	Path result;
	size_t pos = localComponentOrder.back();
	while (true)
	{
		std::vector<Node> addHere = reverse(std::get<1>(predecessor[pos]));
		assert(result.nodes.size() == 0 || addHere[0] == result.nodes.back());
		result.nodes.insert(result.nodes.end(), addHere.begin()+1, addHere.end());
		pos = std::get<2>(predecessor[pos]);
		if (pos == localComponentOrder.back()) break;
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
		for (size_t neighbor : localComponentOutEdges[top])
		{
			checkStack.emplace_back(neighbor);
		}
	}
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

Path getHeavyPath(const GfaGraph& graph, const std::vector<ReadPath>& readPaths, const size_t localResolveLength)
{
	std::vector<size_t> selfDistances;
	size_t distanceSum = 0;
	size_t distanceDivisor = 0;
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		size_t selfDistance = getSelfDistance(graph, i);
		selfDistances.emplace_back(selfDistance);
		if (selfDistance == std::numeric_limits<size_t>::max()) continue;
		if (selfDistance < localResolveLength) continue;
		distanceSum += selfDistance * graph.nodeSeqs[i].size();
		distanceDivisor += graph.nodeSeqs[i].size();
	}
	size_t estimatedMinimalSelfDistance = distanceSum/distanceDivisor/2;
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
	std::vector<size_t> nodeBelongsToStronglyConnectedLocalComponent;
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		nodeBelongsToStronglyConnectedLocalComponent.emplace_back(i);
	}
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (!nodeExists[i]) continue;
		if (selfDistances[i] == std::numeric_limits<size_t>::max()) continue;
		if (selfDistances[i] >= estimatedMinimalSelfDistance) continue;
		assert(graph.edges.count(Node { i, true }) == 1);
		assert(graph.edges.count(Node { i, false }) == 1);
		for (auto edge : graph.edges.at(Node { i, true }))
		{
			Node target = std::get<0>(edge);
			if (!nodeExists[target.id()]) continue;
			if (selfDistances[target.id()] < estimatedMinimalSelfDistance && getShortestPathDistance(graph, target, Node { i, true }) < estimatedMinimalSelfDistance)
			{
				merge(nodeBelongsToStronglyConnectedLocalComponent, i, target.id());
			}
		}
		for (auto edge : graph.edges.at(Node { i, false }))
		{
			Node target = std::get<0>(edge);
			if (!nodeExists[target.id()]) continue;
			if (selfDistances[target.id()] < estimatedMinimalSelfDistance && getShortestPathDistance(graph, target, Node { i, false }) < estimatedMinimalSelfDistance)
			{
				merge(nodeBelongsToStronglyConnectedLocalComponent, i, target.id());
			}
		}
	}
	assert(nodeBelongsToStronglyConnectedLocalComponent.at(bestStartNode) == bestStartNode);
	size_t countComponents = 0;
	{
		phmap::flat_hash_map<size_t, size_t> componentRenaming;
		for (size_t i = 0; i < nodeBelongsToStronglyConnectedLocalComponent.size(); i++)
		{
			if (!nodeExists[i]) continue;
			size_t component = find(nodeBelongsToStronglyConnectedLocalComponent, i);
			if (componentRenaming.count(component) == 1) continue;
			size_t newName = componentRenaming.size();
			componentRenaming[component] = newName;
		}
		for (size_t i = 0; i < nodeBelongsToStronglyConnectedLocalComponent.size(); i++)
		{
			if (!nodeExists[i]) continue;
			nodeBelongsToStronglyConnectedLocalComponent[i] = componentRenaming.at(nodeBelongsToStronglyConnectedLocalComponent[i]);
		}
		countComponents = componentRenaming.size();
	}
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (!nodeExists[i]) continue;
	}
	std::vector<std::vector<size_t>> nodesInLocalComponent;
	nodesInLocalComponent.resize(countComponents);
	for (size_t i = 0; i < nodeBelongsToStronglyConnectedLocalComponent.size(); i++)
	{
		if (!nodeExists[i]) continue;
		nodesInLocalComponent[nodeBelongsToStronglyConnectedLocalComponent[i]].emplace_back(i);
	}
	std::vector<phmap::flat_hash_set<size_t>> localComponentInEdges;
	std::vector<phmap::flat_hash_set<size_t>> localComponentOutEdges;
	localComponentInEdges.resize(countComponents);
	localComponentOutEdges.resize(countComponents);
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (!nodeExists[i]) continue;
		assert(graph.edges.count(Node { i, true }) == 1);
		assert(graph.edges.count(Node { i, false }) == 1);
		size_t thisComponent = nodeBelongsToStronglyConnectedLocalComponent[i];
		for (auto edge : graph.edges.at(Node { i, nodeOrientation[i] }))
		{
			Node target = std::get<0>(edge);
			if (!nodeExists[target.id()]) continue;
			size_t otherComponent = nodeBelongsToStronglyConnectedLocalComponent[target.id()];
			if (otherComponent == thisComponent) continue;
			localComponentOutEdges[thisComponent].emplace(otherComponent);
		}
		for (auto edge : graph.edges.at(Node { i, !nodeOrientation[i] }))
		{
			Node target = std::get<0>(edge);
			if (!nodeExists[target.id()]) continue;
			size_t otherComponent = nodeBelongsToStronglyConnectedLocalComponent[target.id()];
			if (otherComponent == thisComponent) continue;
			localComponentInEdges[thisComponent].emplace(otherComponent);
		}
	}
	std::vector<size_t> localComponentOrder = getLocalComponentTopologicalOrder(localComponentInEdges, localComponentOutEdges, nodeBelongsToStronglyConnectedLocalComponent.at(bestStartNode));
	assert(localComponentOrder.size() == localComponentInEdges.size());
	Path result = getTopologicallyOrderedHeaviestPath(graph, readPaths, nodeExists, selfDistances, estimatedMinimalSelfDistance, nodeOrientation, localComponentOrder, nodeBelongsToStronglyConnectedLocalComponent, nodesInLocalComponent);
	return result;
}
