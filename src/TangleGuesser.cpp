#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include "TangleGuesser.h"

std::string revnode(std::string node)
{
	assert(node.size() >= 2);
	assert(node[0] == '>' || node[0] == '<');
	return (node[0] == '>' ? "<" : ">") + node.substr(1);
}

std::unordered_set<std::string> matchNodes(const KmerMatcher& matcher, const std::string& gfaPath)
{
	std::unordered_set<std::string> result;
	std::ifstream file { gfaPath };
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (!file.good()) break;
		if (line[0] != 'S') continue;
		std::string dummy, nodename, nodeseq;
		std::stringstream sstr { line };
		sstr >> dummy >> nodename >> nodeseq;
		if (nodeseq.size() >= 100000)
		{
			std::string prefix = nodeseq.substr(0, 50000);
			size_t matchLength = matcher.getMatchLength(prefix);
			if (matchLength > 2000)
			{
				result.insert("<" + nodename);
			}
			std::string suffix = nodeseq.substr(nodeseq.size()-50000);
			matchLength = matcher.getMatchLength(suffix);
			if (matchLength > 2000)
			{
				result.insert(">" + nodename);
			}
			continue;
		}
		size_t matchLength = matcher.getMatchLength(nodeseq);
		if (matchLength < 2000) continue;
		result.insert(nodename);
	}
	return result;
}

std::string find(std::unordered_map<std::string, std::string>& parent, std::string key)
{
	while (parent.at(key) != parent.at(parent.at(key))) parent[key] = parent[parent[key]];
	return parent[key];
}

void merge(std::unordered_map<std::string, std::string>& parent, std::string left, std::string right)
{
	left = find(parent, left);
	right = find(parent, right);
	assert(parent.at(left) == left);
	assert(parent.at(right) == right);
	parent[right] = left;
}

std::vector<std::vector<std::string>> extendTangles(const std::unordered_set<std::string>& validNodes, const std::string& gfaPath)
{
	std::unordered_map<std::string, std::string> parent;
	for (auto node : validNodes)
	{
		parent[node] = node;
	}
	phmap::flat_hash_map<std::string, std::vector<std::string>> edgesInvolvingNodeEnds;
	{
		std::ifstream file { gfaPath };
		while (file.good())
		{
			std::string line;
			getline(file, line);
			if (!file.good()) break;
			if (line[0] != 'L') continue;
			std::stringstream sstr { line };
			std::string dummy, fromnode, tonode;
			std::string fromorient, toorient;
			sstr >> dummy >> fromnode >> fromorient >> tonode >> toorient;
			if (validNodes.count(fromnode) == 1 && validNodes.count(tonode) == 1)
			{
				merge(parent, fromnode, tonode);
			}
			if ((fromorient == "+" && validNodes.count(">" + fromnode) == 1) || (fromorient == "-" && validNodes.count("<" + fromnode) == 1))
			{
				if (validNodes.count(tonode) == 1)
				{
					edgesInvolvingNodeEnds[((fromorient == "+") ? ">" : "<") + fromnode].emplace_back(tonode);
				}
			}
			if ((toorient == "+" && validNodes.count("<" + tonode) == 1) || (toorient == "-" && validNodes.count(">" + tonode) == 1))
			{
				if (validNodes.count(fromnode) == 1)
				{
					edgesInvolvingNodeEnds[((toorient == "-") ? ">" : "<") + tonode].emplace_back(fromnode);
				}
			}
		}
	}
	for (auto node : validNodes)
	{
		find(parent, node);
	}
	std::unordered_map<std::string, size_t> tangleSize;
	for (auto node : validNodes)
	{
		tangleSize[find(parent, node)] += 1;
	}
	std::unordered_map<std::string, size_t> tangleNumber;
	size_t numTangles = 0;
	for (auto pair : tangleSize)
	{
		if (pair.second < 10) continue;
		tangleNumber[pair.first] = numTangles;
		numTangles += 1;
	}
	std::vector<std::vector<std::string>> result;
	result.resize(numTangles);
	for (auto node : validNodes)
	{
		auto tangle = find(parent, node);
		if (tangleNumber.count(tangle) == 0) continue;
		result[tangleNumber.at(tangle)].push_back(node);
	}
	for (const auto& pair : edgesInvolvingNodeEnds)
	{
		phmap::flat_hash_set<size_t> possibleTangles;
		for (const auto& node : pair.second)
		{
			assert(parent.count(node) == 1);
			if (tangleNumber.count(find(parent, node)) == 1)
			{
				possibleTangles.insert(tangleNumber.at(find(parent, node)));
			}
		}
		if (possibleTangles.size() != 1) continue;
		result[*possibleTangles.begin()].push_back(pair.first);
	}
	return result;
}

bool hasCycle(const std::vector<std::string>& tangleNodes, const std::unordered_map<std::string, std::vector<std::string>>& edges)
{
	std::unordered_map<std::string, std::unordered_set<std::string>> reachable;
	for (auto node : tangleNodes)
	{
		if (edges.count(">" + node) == 1)
		{
			for (auto edge : edges.at(">" + node))
			{
				reachable[">" + node].insert(edge);
				if (edge == ">" + node) return true;
			}
		}
		if (edges.count("<" + node) == 1)
		{
			for (auto edge : edges.at("<" + node))
			{
				reachable["<" + node].insert(edge);
				if (edge == "<" + node) return true;
			}
		}
	}
	while (true)
	{
		bool changed = false;
		for (auto node : tangleNodes)
		{
			auto oldReachable = reachable[">" + node];
			for (auto node2 : oldReachable)
			{
				if (edges.count(node2) == 0) continue;
				for (auto edge : edges.at(node2))
				{
					if (reachable[">" + node].count(edge) == 1) continue;
					reachable[">" + node].insert(edge);
					changed = true;
					if (edge == ">" + node) return true;
				}
			}
			oldReachable = reachable["<" + node];
			for (auto node2 : oldReachable)
			{
				if (edges.count(node2) == 0) continue;
				for (auto edge : edges.at(node2))
				{
					if (reachable["<" + node].count(edge) == 1) continue;
					reachable["<" + node].insert(edge);
					changed = true;
					if (edge == "<" + node) return true;
				}
			}
		}
		if (!changed) break;
	}
	return false;
}

void filterOutAcyclic(std::vector<std::vector<std::string>>& tangleNodes, const std::string& gfaPath)
{
	std::unordered_map<std::string, std::vector<std::string>> edges;
	{
		std::ifstream file { gfaPath };
		while (file.good())
		{
			std::string line;
			getline(file, line);
			if (!file.good()) break;
			if (line[0] != 'L') continue;
			std::stringstream sstr { line };
			std::string dummy, fromnode, tonode;
			std::string fromorient, toorient;
			sstr >> dummy >> fromnode >> fromorient >> tonode >> toorient;
			assert(fromorient == "+" || fromorient == "-");
			assert(toorient == "+" || toorient == "-");
			fromnode = (fromorient == "+" ? ">" : "<") + fromnode;
			tonode = (toorient == "+" ? ">" : "<") + tonode;
			edges[fromnode].push_back(tonode);
			edges[revnode(tonode)].push_back(revnode(fromnode));
		}
	}
	for (size_t i = tangleNodes.size()-1; i < tangleNodes.size(); i--)
	{
		if (hasCycle(tangleNodes[i], edges)) continue;
		std::swap(tangleNodes[i], tangleNodes.back());
		tangleNodes.pop_back();
	}
}

std::vector<std::vector<std::string>> guessTangles(const KmerMatcher& matcher, const std::string& gfaPath)
{
	std::cerr << "guess tangles from " << gfaPath << std::endl;
	auto nodes = matchNodes(matcher, gfaPath);
	auto tangles = extendTangles(nodes, gfaPath);
	filterOutAcyclic(tangles, gfaPath);
	for (size_t i = 0; i < tangles.size(); i++)
	{
		assert(tangles[i].size() >= 1);
		std::sort(tangles[i].begin(), tangles[i].end());
	}
	std::sort(tangles.begin(), tangles.end(), [](const auto& left, const auto& right)
	{
		if (left.size() > right.size()) return true;
		if (left.size() < right.size()) return false;
		for (size_t i = 0; i < left.size(); i++)
		{
			if (left[i] > right[i]) return true;
			if (left[i] < right[i]) return false;
		}
		return false;
	});
	return tangles;
}
