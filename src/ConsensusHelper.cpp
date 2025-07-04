#include "ConsensusHelper.h"
#include "SequenceAligner.h"
#include "WfaHelper.h"
#include "Logger.h"

std::vector<std::tuple<size_t, size_t, std::string, size_t>> getEdits(const std::string_view& refSequence, const std::string_view& querySequence, size_t maxEdits)
{
	if (refSequence == querySequence) return std::vector<std::tuple<size_t, size_t, std::string, size_t>>{};
	// special case for very small differences
	{
		size_t countEqualAtStart = 0;
		size_t countEqualAtEnd = 0;
		while (countEqualAtStart < refSequence.size() && countEqualAtStart < querySequence.size() && refSequence[countEqualAtStart] == querySequence[countEqualAtStart]) countEqualAtStart += 1;
		while (countEqualAtEnd < refSequence.size() && countEqualAtEnd < querySequence.size() && refSequence[refSequence.size()-1-countEqualAtEnd] == querySequence[querySequence.size()-1-countEqualAtEnd]) countEqualAtEnd += 1;
		if (refSequence.size() == querySequence.size() && countEqualAtEnd+countEqualAtStart+1 == refSequence.size())
		{
			// just a single SNP
			std::vector<std::tuple<size_t, size_t, std::string, size_t>> result;
			result.emplace_back(countEqualAtStart, countEqualAtStart+1, querySequence.substr(countEqualAtStart, 1), countEqualAtStart);
			return result;
		}
		if (querySequence.size() > refSequence.size() && countEqualAtEnd+countEqualAtStart == refSequence.size())
		{
			// just a single insertion
			std::vector<std::tuple<size_t, size_t, std::string, size_t>> result;
			result.emplace_back(countEqualAtStart, countEqualAtStart, querySequence.substr(countEqualAtStart, querySequence.size() - countEqualAtStart - countEqualAtEnd), countEqualAtStart);
			return result;
		}
		if (querySequence.size() < refSequence.size() && countEqualAtEnd+countEqualAtStart == querySequence.size())
		{
			// just a single deletion
			std::vector<std::tuple<size_t, size_t, std::string, size_t>> result;
			result.emplace_back(countEqualAtStart, refSequence.size()-countEqualAtEnd, "", countEqualAtStart);
			return result;
		}
	}
	// reference is horizontal, first bp left, last bp right (one bp one column)
	// query is vertical, first bp top, last bp bottom (one bp one row)
	// edits moves diagonally down-right
	// j moves horizontally, plus right minus left
	// insertion is extra sequence on query (vertical backtrace)
	// deletion is missing sequence on query (horizontal backtrace)
	// wfa matrix pair coord first is edits (row), second is j (column)
	// wfa matrix is triangular, missing parts are implied: edit 0 j 0 is above edit 1 j 1 and up-right from edit 1 j 0
	// real matrix pair coord first is ref pos (column), second is query pos (row)
	const std::pair<size_t, size_t> uninitialized { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
	const size_t mismatchCost = 2;
	const size_t gapOpenCost = 2; // only open cost, and 1bp indel costs gapOpenCost+gapExtendCost
	const size_t gapExtendCost = 1;
	// match cost 0
	// gap close cost 0
	Wfa::WfaMatrix matchMatrix;
	Wfa::WfaMatrix insertionMatrix;
	Wfa::WfaMatrix deletionMatrix;
	matchMatrix.resize(1);
	insertionMatrix.resize(1);
	deletionMatrix.resize(1);
	insertionMatrix[0].resize(1, uninitialized);
	deletionMatrix[0].resize(1, uninitialized);
	matchMatrix[0].resize(1);
	matchMatrix[0][0].first = 0;
	matchMatrix[0][0].second = 0;
	while (matchMatrix[0][0].second < refSequence.size() && matchMatrix[0][0].second < querySequence.size() && refSequence[matchMatrix[0][0].second] == querySequence[matchMatrix[0][0].second]) matchMatrix[0][0].second += 1;
	size_t edits = 0;
	std::pair<size_t, size_t> backtraceWfaPosition = uninitialized;
	for (edits = 1; edits < maxEdits; edits++)
	{
		matchMatrix.emplace_back();
		insertionMatrix.emplace_back();
		deletionMatrix.emplace_back();
		matchMatrix.back().resize(2*edits+1, uninitialized);
		insertionMatrix.back().resize(2*edits+1, uninitialized);
		deletionMatrix.back().resize(2*edits+1, uninitialized);
		for (size_t j = 0; j < 2*edits+1; j++)
		{
			//mismatch
			Wfa::updateMatrix(std::make_pair(edits, j), mismatchCost, 1, 1, matchMatrix, matchMatrix, refSequence.size(), querySequence.size());
			// gap open
			Wfa::updateMatrix(std::make_pair(edits, j), gapOpenCost, 0, 0, matchMatrix, insertionMatrix, refSequence.size(), querySequence.size());
			Wfa::updateMatrix(std::make_pair(edits, j), gapOpenCost, 0, 0, matchMatrix, deletionMatrix, refSequence.size(), querySequence.size());
			// gap extend
			Wfa::updateMatrix(std::make_pair(edits, j), gapExtendCost, 1, 0, insertionMatrix, insertionMatrix, refSequence.size(), querySequence.size());
			Wfa::updateMatrix(std::make_pair(edits, j), gapExtendCost, 0, 1, deletionMatrix, deletionMatrix, refSequence.size(), querySequence.size());
			insertionMatrix[edits][j].second = insertionMatrix[edits][j].first;
			deletionMatrix[edits][j].second = deletionMatrix[edits][j].first;
			// gap close
			Wfa::updateMatrix(std::make_pair(edits, j), 0, 0, 0, insertionMatrix, matchMatrix, refSequence.size(), querySequence.size());
			Wfa::updateMatrix(std::make_pair(edits, j), 0, 0, 0, deletionMatrix, matchMatrix, refSequence.size(), querySequence.size());
			//matches
			matchMatrix[edits][j].second = matchMatrix[edits][j].first;
			size_t offset = 0;
			if (matchMatrix[edits][j] == uninitialized) continue;
			std::pair<size_t, size_t> realPos = Wfa::getRealMatrixPosition(std::make_pair(edits, j), matchMatrix[edits][j].first);
			assert(realPos.first <= refSequence.size() && realPos.second <= querySequence.size());
			while (realPos.first+offset < refSequence.size() && realPos.second+offset < querySequence.size() && refSequence[realPos.first+offset] == querySequence[realPos.second+offset]) offset += 1;
			matchMatrix[edits][j].second += offset;
			if (realPos.first+offset == refSequence.size() && realPos.second+offset == querySequence.size())
			{
				backtraceWfaPosition.first = edits;
				backtraceWfaPosition.second = j;
				break;
			}
		}
		if (backtraceWfaPosition != uninitialized) break;
	}
	if (backtraceWfaPosition == uninitialized) return std::vector<std::tuple<size_t, size_t, std::string, size_t>> {};
	assert(backtraceWfaPosition != uninitialized);
	std::vector<std::tuple<size_t, size_t, std::string, size_t>> result;
	// 0 match 1 insertion 2 deletion
	size_t currentMatrix = 0;
	std::pair<size_t, size_t> lastMatchRealPosition = uninitialized;
	while (backtraceWfaPosition.first != 0)
	{
		switch(currentMatrix)
		{
			case 0:
			{
				assert(backtraceWfaPosition.first < matchMatrix.size());
				assert(backtraceWfaPosition.second < matchMatrix[backtraceWfaPosition.first].size());
				assert(matchMatrix[backtraceWfaPosition.first][backtraceWfaPosition.second] != uninitialized);
				std::pair<size_t, size_t> newMatchRealPosition = Wfa::getRealMatrixPosition(backtraceWfaPosition, matchMatrix[backtraceWfaPosition.first][backtraceWfaPosition.second].second);
				if (lastMatchRealPosition != uninitialized)
				{
					assert(newMatchRealPosition.first <= lastMatchRealPosition.first);
					assert(newMatchRealPosition.second <= lastMatchRealPosition.second);
					result.emplace_back(newMatchRealPosition.first, lastMatchRealPosition.first, querySequence.substr(newMatchRealPosition.second, lastMatchRealPosition.second - newMatchRealPosition.second), newMatchRealPosition.second);
				}
				lastMatchRealPosition = Wfa::getRealMatrixPosition(backtraceWfaPosition, matchMatrix[backtraceWfaPosition.first][backtraceWfaPosition.second].first);
				if (Wfa::canBacktrace(backtraceWfaPosition, mismatchCost, 1, 1, matchMatrix, matchMatrix, refSequence.size(), querySequence.size()))
				{
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, mismatchCost, 1, 1, matchMatrix, matchMatrix);
					break;
				}
				if (Wfa::canBacktrace(backtraceWfaPosition, 0, 0, 0, insertionMatrix, matchMatrix, refSequence.size(), querySequence.size()))
				{
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, 0, 0, 0, matchMatrix, insertionMatrix);
					currentMatrix = 1;
					break;
				}
				else
				{
					assert(Wfa::canBacktrace(backtraceWfaPosition, 0, 0, 0, deletionMatrix, matchMatrix, refSequence.size(), querySequence.size()));
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, 0, 0, 0, matchMatrix, deletionMatrix);
					currentMatrix = 2;
					break;
				}
			}
			case 1:
			{
				assert(backtraceWfaPosition.first < insertionMatrix.size());
				assert(backtraceWfaPosition.second < insertionMatrix[backtraceWfaPosition.first].size());
				assert(insertionMatrix[backtraceWfaPosition.first][backtraceWfaPosition.second] != uninitialized);
				if (Wfa::canBacktrace(backtraceWfaPosition, gapExtendCost, 1, 0, insertionMatrix, insertionMatrix, refSequence.size(), querySequence.size()))
				{
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, gapExtendCost, 1, 0, insertionMatrix, insertionMatrix);
					break;
				}
				else
				{
					assert(Wfa::canBacktrace(backtraceWfaPosition, gapOpenCost, 0, 0, matchMatrix, insertionMatrix, refSequence.size(), querySequence.size()));
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, gapOpenCost, 0, 0, insertionMatrix, matchMatrix);
					currentMatrix = 0;
					break;
				}
			}
			case 2:
			{
				assert(backtraceWfaPosition.first < deletionMatrix.size());
				assert(backtraceWfaPosition.second < deletionMatrix[backtraceWfaPosition.first].size());
				assert(deletionMatrix[backtraceWfaPosition.first][backtraceWfaPosition.second] != uninitialized);
				if (Wfa::canBacktrace(backtraceWfaPosition, gapExtendCost, 0, 1, deletionMatrix, deletionMatrix, refSequence.size(), querySequence.size()))
				{
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, gapExtendCost, 0, 1, deletionMatrix, deletionMatrix);
					break;
				}
				else
				{
					assert(Wfa::canBacktrace(backtraceWfaPosition, gapOpenCost, 0, 0, matchMatrix, deletionMatrix, refSequence.size(), querySequence.size()));
					backtraceWfaPosition = Wfa::getBacktracePosition(backtraceWfaPosition, gapOpenCost, 0, 0, deletionMatrix, matchMatrix);
					currentMatrix = 0;
					break;
				}
			}
		}
	}
	std::pair<size_t, size_t> newMatchRealPosition = Wfa::getRealMatrixPosition(backtraceWfaPosition, matchMatrix[backtraceWfaPosition.first][backtraceWfaPosition.second].second);
	if (lastMatchRealPosition != uninitialized)
	{
		assert(newMatchRealPosition.first <= lastMatchRealPosition.first);
		assert(newMatchRealPosition.second <= lastMatchRealPosition.second);
		result.emplace_back(newMatchRealPosition.first, lastMatchRealPosition.first, querySequence.substr(newMatchRealPosition.second, lastMatchRealPosition.second - newMatchRealPosition.second), newMatchRealPosition.second);
	}
	assert(result.size() >= 1);
	std::reverse(result.begin(), result.end());
	return result;
}

void removeRepeatKmers(std::vector<std::pair<uint64_t, uint64_t>>& matches)
{
	phmap::flat_hash_set<size_t> foundOnce;
	phmap::flat_hash_set<size_t> foundTwice;
	for (const auto& pair : matches)
	{
		if (foundOnce.count(pair.first) == 1) foundTwice.emplace(pair.first);
		foundOnce.emplace(pair.first);
	}
	if (foundTwice.size() == 0) return;
	for (size_t i = matches.size()-1; i < matches.size(); i--)
	{
		if (foundTwice.count(matches[i].first) == 0) continue;
		std::swap(matches[i], matches.back());
		matches.pop_back();
	}
	std::sort(matches.begin(), matches.end());
}

void removeRepeatKmersEitherRepeat(std::vector<std::pair<uint64_t, uint64_t>>& matches)
{
	phmap::flat_hash_set<size_t> foundOnceFirst;
	phmap::flat_hash_set<size_t> foundTwiceFirst;
	phmap::flat_hash_set<size_t> foundOnceSecond;
	phmap::flat_hash_set<size_t> foundTwiceSecond;
	for (const auto& pair : matches)
	{
		if (foundOnceFirst.count(pair.first) == 1) foundTwiceFirst.emplace(pair.first);
		foundOnceFirst.emplace(pair.first);
		if (foundOnceSecond.count(pair.second) == 1) foundTwiceSecond.emplace(pair.second);
		foundOnceSecond.emplace(pair.second);
	}
	if (foundTwiceFirst.size() == 0 && foundTwiceSecond.size() == 0) return;
	for (size_t i = matches.size()-1; i < matches.size(); i--)
	{
		if (foundTwiceFirst.count(matches[i].first) == 0 && foundTwiceSecond.count(matches[i].second) == 0) continue;
		std::swap(matches[i], matches.back());
		matches.pop_back();
	}
	std::sort(matches.begin(), matches.end());
}

void removeNonDiagonalKmers(std::vector<std::tuple<size_t, size_t, size_t>>& kmerMatches)
{
	std::sort(kmerMatches.begin(), kmerMatches.end(), [](auto left, auto right){ return (int)std::get<0>(left) - (int)std::get<1>(left) < (int)std::get<0>(right) - (int)std::get<1>(right); });
	size_t bestClusterSize = 0;
	int bestClusterStart = 0;
	int bestClusterEnd = 0;
	size_t currentClusterStart = 0;
	for (size_t i = 1; i <= kmerMatches.size(); i++)
	{
		if (i == kmerMatches.size() || (int)std::get<0>(kmerMatches[i]) - (int)std::get<1>(kmerMatches[i]) > (int)std::get<0>(kmerMatches[i-1]) - (int)std::get<1>(kmerMatches[i-1]) + 100)
		{
			size_t currentClusterSize = i - currentClusterStart;
			if (currentClusterSize > bestClusterSize)
			{
				bestClusterSize = currentClusterSize;
				bestClusterStart = (int)std::get<0>(kmerMatches[currentClusterStart]) - (int)std::get<1>(kmerMatches[currentClusterStart]);
				bestClusterEnd = (int)std::get<0>(kmerMatches[i-1]) - (int)std::get<1>(kmerMatches[i-1]);
			}
			currentClusterSize = 0;
			currentClusterStart = i;
		}
	}
	for (size_t i = kmerMatches.size()-1; i < kmerMatches.size(); i--)
	{
		int diagonal = std::get<0>(kmerMatches[i]) - std::get<1>(kmerMatches[i]);
		if (diagonal >= bestClusterStart && diagonal <= bestClusterEnd) continue;
		std::swap(kmerMatches[i], kmerMatches.back());
		kmerMatches.pop_back();
	}
	std::sort(kmerMatches.begin(), kmerMatches.end());
}

void removeNonDiagonalKmers(std::vector<std::pair<size_t, size_t>>& kmerMatches)
{
	std::sort(kmerMatches.begin(), kmerMatches.end(), [](auto left, auto right){ return (int)left.first - (int)left.second < (int)right.first - (int)right.second; });
	size_t bestClusterSize = 0;
	int bestClusterStart = 0;
	int bestClusterEnd = 0;
	size_t currentClusterStart = 0;
	for (size_t i = 1; i <= kmerMatches.size(); i++)
	{
		if (i == kmerMatches.size() || (int)kmerMatches[i].first - (int)kmerMatches[i].second > (int)kmerMatches[i-1].first - (int)kmerMatches[i-1].second + 100)
		{
			size_t currentClusterSize = i - currentClusterStart;
			if (currentClusterSize > bestClusterSize)
			{
				bestClusterSize = currentClusterSize;
				bestClusterStart = (int)kmerMatches[currentClusterStart].first - (int)kmerMatches[currentClusterStart].second;
				bestClusterEnd = (int)kmerMatches[i-1].first - (int)kmerMatches[i-1].second;
			}
			currentClusterSize = 0;
			currentClusterStart = i;
		}
	}
	for (size_t i = kmerMatches.size()-1; i < kmerMatches.size(); i--)
	{
		int diagonal = (int)kmerMatches[i].first - (int)kmerMatches[i].second;
		if (diagonal >= bestClusterStart && diagonal <= bestClusterEnd) continue;
		std::swap(kmerMatches[i], kmerMatches.back());
		kmerMatches.pop_back();
	}
	std::sort(kmerMatches.begin(), kmerMatches.end());
}

void removeNonDiagonalKmers(std::vector<std::pair<size_t, size_t>>& kmerMatches, const size_t kmerSize)
{
	if (kmerMatches.size() == 0) return;
	std::sort(kmerMatches.begin(), kmerMatches.end(), [](auto left, auto right){ return left.first < right.first; });
	std::vector<size_t> parent;
	parent.reserve(kmerMatches.size());
	for (size_t i = 0; i < kmerMatches.size(); i++)
	{
		parent.emplace_back(i);
	}
	for (size_t i = 1; i < kmerMatches.size(); i++)
	{
		for (size_t j = i-1; j < kmerMatches.size(); j--)
		{
			if (kmerMatches[j].first+200 < kmerMatches[i].first) break;
			if ((int)kmerMatches[i].first - (int)kmerMatches[i].second > (int)kmerMatches[j].first - (int)kmerMatches[j].second + 10) continue;
			if ((int)kmerMatches[i].first - (int)kmerMatches[i].second < (int)kmerMatches[j].first - (int)kmerMatches[j].second - 10) continue;
			while (parent[i] != parent[parent[i]]) parent[i] = parent[parent[i]];
			while (parent[j] != parent[parent[j]]) parent[j] = parent[parent[j]];
			parent[parent[j]] = parent[i];
		}
	}
	phmap::flat_hash_map<size_t, size_t> clusterKmerCount;
	for (size_t i = 0; i < parent.size(); i++)
	{
		while (parent[i] != parent[parent[i]]) parent[i] = parent[parent[i]];
		clusterKmerCount[parent[i]] += 1;
	}
	std::vector<bool> removeMatch;
	removeMatch.resize(kmerMatches.size(), false);
	for (size_t i = 0; i < parent.size(); i++)
	{
		assert(parent[i] == parent[parent[i]]);
		if (clusterKmerCount.at(parent[i]) < 10) removeMatch[i] = true;
	}
	for (size_t i = kmerMatches.size()-1; i < kmerMatches.size(); i--)
	{
		if (!removeMatch[i]) continue;
		std::swap(kmerMatches[i], kmerMatches.back());
		kmerMatches.pop_back();
	}
	std::sort(kmerMatches.begin(), kmerMatches.end());
}
/*
void removeNonDiagonalKmers(std::vector<std::pair<size_t, size_t>>& kmerMatches, const size_t kmerSize)
{
	if (kmerMatches.size() == 0) return;
	std::sort(kmerMatches.begin(), kmerMatches.end(), [](auto left, auto right){ return (int)left.first - (int)left.second < (int)right.first - (int)right.second; });
	std::vector<bool> removeMatch;
	removeMatch.resize(kmerMatches.size(), false);
	size_t bestClusterSize = 0;
	int bestClusterStart = 0;
	int bestClusterEnd = 0;
	size_t currentClusterStart = 0;
	for (size_t i = 1; i <= kmerMatches.size(); i++)
	{
		if (i < kmerMatches.size() && (int)kmerMatches[i].first - (int)kmerMatches[i].second < (int)kmerMatches[i-1].first - (int)kmerMatches[i-1].second + 10) continue;
		std::sort(kmerMatches.begin()+currentClusterStart, kmerMatches.begin()+i);
		size_t currentSubclusterStart = currentClusterStart;
		size_t currentMatchLength = 0;
		size_t currentMatchEnd = 0;
		for (size_t j = currentClusterStart; j < i; j++)
		{
			if (kmerMatches[j].first > currentMatchEnd+500)
			{
				if (currentMatchLength < 100)
				{
					for (size_t k = currentSubclusterStart; k < j; k++)
					{
						removeMatch[k] = true;
					}
				}
				currentMatchLength = 0;
				currentSubclusterStart = j;
			}
			if (kmerMatches[j].first > currentMatchEnd)
			{
				currentMatchLength += kmerSize;
				currentMatchEnd = kmerMatches[j].first+kmerSize;
			}
			else if (kmerMatches[j].first+kmerSize > currentMatchEnd)
			{
				currentMatchLength += kmerMatches[j].first+kmerSize - currentMatchEnd;
				currentMatchEnd = kmerMatches[j].first+kmerSize;
			}
		}
		if (currentMatchLength < 100)
		{
			for (size_t k = currentSubclusterStart; k < i; k++)
			{
				removeMatch[k] = true;
			}
		}
		currentClusterStart = i;
	}
	for (size_t i = kmerMatches.size()-1; i < kmerMatches.size(); i--)
	{
		if (!removeMatch[i]) continue;
		std::swap(kmerMatches[i], kmerMatches.back());
		kmerMatches.pop_back();
	}
	std::sort(kmerMatches.begin(), kmerMatches.end());
}
*/
std::unordered_set<size_t> getMajorityNodes(const std::vector<OntLoop>& paths)
{
	std::unordered_map<size_t, size_t> nodeCoverage;
	std::unordered_set<size_t> notCore;
	std::unordered_set<size_t> firstNodes;
	std::unordered_set<size_t> lastNodes;
	for (const auto& path : paths)
	{
		firstNodes.insert(path.path[0].id());
		lastNodes.insert(path.path.back().id());
		std::unordered_map<size_t, size_t> coverageHere;
		for (auto node : path.path)
		{
			coverageHere[node.id()] += 1;
		}
		for (auto pair : coverageHere)
		{
			if (pair.second != 1) notCore.insert(pair.first);
			nodeCoverage[pair.first] += 1;
		}
	}
	std::unordered_set<size_t> result;
	for (auto pair : nodeCoverage)
	{
		if (notCore.count(pair.first) == 1) continue;
		if (firstNodes.count(pair.first) == 1 && lastNodes.count(pair.first) == 1) continue;
		if (pair.second > paths.size()/2) result.insert(pair.first);
	}
	return result;
}

std::vector<size_t> getMajorityPath(const std::unordered_set<size_t>& coreNodes, const std::vector<OntLoop>& rawPaths)
{
	std::unordered_map<size_t, size_t> firstCount;
	std::unordered_map<size_t, size_t> coverage;
	std::unordered_map<size_t, std::unordered_map<size_t, size_t>> edges;
	for (const auto& path : rawPaths)
	{
		size_t lastCore = std::numeric_limits<size_t>::max();
		for (auto node : path.path)
		{
			if (coreNodes.count(node.id()) == 0) continue;
			coverage[node.id()] += 1;
			if (lastCore == std::numeric_limits<size_t>::max())
			{
				firstCount[node.id()] += 1;
			}
			else
			{
				edges[lastCore][node.id()] += 1;
			}
			lastCore = node.id();
		}
	}
	std::vector<size_t> result;
	size_t startNode = std::numeric_limits<size_t>::max();
	for (auto pair : firstCount)
	{
		if (startNode == std::numeric_limits<size_t>::max() || pair.second > firstCount.at(startNode))
		{
			startNode = pair.first;
		}
	}
	result.emplace_back(startNode);
	size_t pos = startNode;
	while (true)
	{
		if (edges.count(pos) == 0) break;
		size_t maxNode = std::numeric_limits<size_t>::max();
		for (auto pair : edges.at(pos))
		{
			if (maxNode == std::numeric_limits<size_t>::max() || pair.second > edges.at(pos).at(maxNode))
			{
				maxNode = pair.first;
			}
		}
		assert(maxNode != std::numeric_limits<size_t>::max());
		pos = maxNode;
		result.emplace_back(pos);
	}
	return result;
}

std::vector<std::pair<size_t, size_t>> getIncreasingChain(const std::vector<std::pair<size_t, size_t>>& kmerMatches)
{
	if (kmerMatches.size() == 0) return kmerMatches;
	std::vector<std::vector<size_t>> paretoFrontier;
	paretoFrontier.emplace_back();
	paretoFrontier.emplace_back();
	paretoFrontier[1].emplace_back(0);
	for (size_t i = 1; i < kmerMatches.size(); i++)
	{
		bool found = false;
		for (size_t length = paretoFrontier.size()-1; length < paretoFrontier.size(); length++)
		{
			for (size_t j : paretoFrontier[length])
			{
				assert(kmerMatches[i].first > kmerMatches[j].first);
				assert(kmerMatches[i].second != kmerMatches[j].second);
				if (kmerMatches[i].second < kmerMatches[j].second) continue;
				found = true;
				break;
			}
			if (found)
			{
				assert(length+1 <= paretoFrontier.size());
				if (length+1 == paretoFrontier.size()) paretoFrontier.emplace_back();
				bool thisIsInParetoFrontier = true;
				for (size_t j : paretoFrontier[length+1])
				{
					assert(kmerMatches[i].first > kmerMatches[j].first);
					assert(kmerMatches[i].second != kmerMatches[j].second);
					if (kmerMatches[i].second < kmerMatches[j].second) continue;
					thisIsInParetoFrontier = false;
					break;
				}
				if (thisIsInParetoFrontier) paretoFrontier[length+1].emplace_back(i);
				break;
			}
		}
		if (!found)
		{
			bool thisIsInParetoFrontier = true;
			for (size_t j : paretoFrontier[1])
			{
				assert(kmerMatches[i].first > kmerMatches[j].first);
				assert(kmerMatches[i].second != kmerMatches[j].second);
				if (kmerMatches[i].second < kmerMatches[j].second) continue;
				thisIsInParetoFrontier = false;
				break;
			}
			if (thisIsInParetoFrontier) paretoFrontier[1].emplace_back(i);
		}
	}
	assert(paretoFrontier.back().size() >= 1);
	size_t maxIndex = paretoFrontier.back()[0];
	std::vector<std::pair<size_t, size_t>> result;
	result.emplace_back(kmerMatches[maxIndex]);
	for (size_t length = paretoFrontier.size()-2; length > 0; length--)
	{
		bool found = false;
		for (size_t j : paretoFrontier[length])
		{
			if (kmerMatches[j].first > kmerMatches[maxIndex].first) continue;
			if (kmerMatches[j].second > kmerMatches[maxIndex].second) continue;
			found = true;
			result.emplace_back(kmerMatches[j]);
			maxIndex = j;
			break;
		}
		assert(found);
	}
	std::reverse(result.begin(), result.end());
	return result;
}

std::string getSequence(const std::vector<Node>& nodes, const std::vector<std::string>& nodeSeqs, const std::vector<std::string>& revCompNodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges)
{
	std::string result;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		std::string add;
		if (nodes[i].forward())
		{
			add = nodeSeqs.at(nodes[i].id());
		}
		else
		{
			add = revCompNodeSeqs.at(nodes[i].id());
		}
		if (i > 0)
		{
			size_t overlap = std::numeric_limits<size_t>::max();
			for (auto edge : edges.at(nodes[i-1]))
			{
				if (std::get<0>(edge) != nodes[i]) continue;
				assert(overlap == std::numeric_limits<size_t>::max());
				overlap = std::get<1>(edge);
			}
			assert(overlap != std::numeric_limits<size_t>::max());
			add = add.substr(overlap);
		}
		result += add;
	}
	return result;
}

template <typename F>
void iterateSuccessorKmers(const __uint128_t kmer, const size_t k, F callback)
{
	const __uint128_t mask = ((__uint128_t)1 << (__uint128_t)(2ull*k))-(__uint128_t)1;
	for (size_t i = 0; i < 4; i++)
	{
		callback(((kmer << 2) & mask) + i);
	}
}

template <typename F>
void iteratePredecessorKmers(const __uint128_t kmer, const size_t k, F callback)
{
	for (size_t i = 0; i < 4; i++)
	{
		callback((kmer >> 2) + ((__uint128_t)i << (__uint128_t)(2*(k-1))));
	}
}
bool isTipFw(const __uint128_t kmer, const phmap::flat_hash_set<__uint128_t>& keptKmers, const phmap::flat_hash_map<__uint128_t, size_t>& kmerCounts, const size_t k)
{
	size_t countNeighbors = 0;
	iterateSuccessorKmers(kmer, k, [&keptKmers, &countNeighbors](const __uint128_t neighbor)
	{
		if (keptKmers.count(neighbor) == 1) countNeighbors += 1;
	});
	return countNeighbors == 0;
}

bool isTipBw(const __uint128_t kmer, const phmap::flat_hash_set<__uint128_t>& keptKmers, const phmap::flat_hash_map<__uint128_t, size_t>& kmerCounts, const size_t k)
{
	size_t countNeighbors = 0;
	iteratePredecessorKmers(kmer, k, [&keptKmers, &countNeighbors](const __uint128_t neighbor)
	{
		if (keptKmers.count(neighbor) == 1) countNeighbors += 1;
	});
	return countNeighbors == 0;
}

void checkCanRemoveBw(phmap::flat_hash_set<__uint128_t>& removeThese, const __uint128_t kmer, const phmap::flat_hash_set<__uint128_t>& keptKmers, const phmap::flat_hash_map<__uint128_t, size_t>& kmerCounts, const size_t k)
{
	phmap::flat_hash_set<__uint128_t> visited;
	visited.insert(kmer);
	__uint128_t wanderKmer = kmer;
	for (size_t iteration = 0; iteration < k; iteration++)
	{
		std::vector<__uint128_t> allPredecessors;
		iteratePredecessorKmers(wanderKmer, k, [&keptKmers, &allPredecessors](const __uint128_t testKmer)
		{
			if (keptKmers.count(testKmer) == 0) return;
			allPredecessors.emplace_back(testKmer);
		});
		if (allPredecessors.size() == 0) break;
		if (allPredecessors.size() >= 2)
		{
			bool allPredecessorsHaveGoodAlternate = true;
			for (__uint128_t predecessor : allPredecessors)
			{
				bool thisPredecessorHasGoodAlternate = false;
				iterateSuccessorKmers(predecessor, k, [&keptKmers, &kmerCounts, &thisPredecessorHasGoodAlternate](const __uint128_t neighbor)
				{
					if (keptKmers.count(neighbor) == 0) return;
					if (kmerCounts.at(neighbor) >= 3) thisPredecessorHasGoodAlternate = true;
				});
				if (!thisPredecessorHasGoodAlternate) allPredecessorsHaveGoodAlternate = false;
			}
			if (!allPredecessorsHaveGoodAlternate) return;
			break;
		}
		assert(allPredecessors.size() == 1);
		bool predecessorHasOtherNeighbor = false;
		bool predecessorHasBetterNeigbor = false;
		iterateSuccessorKmers(allPredecessors[0], k, [&keptKmers, &kmerCounts, &predecessorHasOtherNeighbor, &predecessorHasBetterNeigbor, wanderKmer](const __uint128_t neighbor)
		{
			if (keptKmers.count(neighbor) == 0) return;
			if (neighbor == wanderKmer) return;
			predecessorHasOtherNeighbor = true;
			if (kmerCounts.at(neighbor) >= 3) predecessorHasBetterNeigbor = true;
		});
		if (predecessorHasOtherNeighbor && predecessorHasBetterNeigbor) break;
		if (predecessorHasOtherNeighbor && !predecessorHasBetterNeigbor) return;
		if (kmerCounts.at(allPredecessors[0]) >= 3) return;
		wanderKmer = allPredecessors[0];
		visited.insert(wanderKmer);
	}
	assert(visited.size() >= 1);
	for (__uint128_t kmer : visited)
	{
		assert(keptKmers.count(kmer) == 1);
		assert(kmerCounts.at(kmer) == 2);
		removeThese.insert(kmer);
	}
}

void checkCanRemoveFw(phmap::flat_hash_set<__uint128_t>& removeThese, const __uint128_t kmer, const phmap::flat_hash_set<__uint128_t>& keptKmers, const phmap::flat_hash_map<__uint128_t, size_t>& kmerCounts, const size_t k)
{
	phmap::flat_hash_set<__uint128_t> visited;
	visited.insert(kmer);
	__uint128_t wanderKmer = kmer;
	for (size_t iteration = 0; iteration < k; iteration++)
	{
		std::vector<__uint128_t> allSuccessors;
		iterateSuccessorKmers(wanderKmer, k, [&keptKmers, &allSuccessors](const __uint128_t testKmer)
		{
			if (keptKmers.count(testKmer) == 0) return;
			allSuccessors.emplace_back(testKmer);
		});
		if (allSuccessors.size() == 0) break;
		if (allSuccessors.size() >= 2)
		{
			bool allSuccessorsHaveGoodAlternate = true;
			for (__uint128_t successor : allSuccessors)
			{
				bool thisSuccessorHasGoodAlternate = false;
				iteratePredecessorKmers(successor, k, [&keptKmers, &kmerCounts, &thisSuccessorHasGoodAlternate](const __uint128_t neighbor)
				{
					if (keptKmers.count(neighbor) == 0) return;
					if (kmerCounts.at(neighbor) >= 3) thisSuccessorHasGoodAlternate = true;
				});
				if (!thisSuccessorHasGoodAlternate) allSuccessorsHaveGoodAlternate = false;
			}
			if (!allSuccessorsHaveGoodAlternate) return;
			break;
		}
		assert(allSuccessors.size() == 1);
		bool successorHasOtherNeighbor = false;
		bool successorHasBetterNeigbor = false;
		iteratePredecessorKmers(allSuccessors[0], k, [&keptKmers, &kmerCounts, &successorHasOtherNeighbor, &successorHasBetterNeigbor, wanderKmer](const __uint128_t neighbor)
		{
			if (keptKmers.count(neighbor) == 0) return;
			if (neighbor == wanderKmer) return;
			successorHasOtherNeighbor = true;
			if (kmerCounts.at(neighbor) >= 3) successorHasBetterNeigbor = true;
		});
		if (successorHasOtherNeighbor && successorHasBetterNeigbor) break;
		if (successorHasOtherNeighbor && !successorHasBetterNeigbor) return;
		if (kmerCounts.at(allSuccessors[0]) >= 3) return;
		wanderKmer = allSuccessors[0];
		visited.insert(wanderKmer);
	}
	assert(visited.size() >= 1);
	for (__uint128_t kmer : visited)
	{
		assert(keptKmers.count(kmer) == 1);
		assert(kmerCounts.at(kmer) == 2);
		removeThese.insert(kmer);
	}
}

void addPath(std::string& result, const phmap::flat_hash_set<__uint128_t>& keptKmers, const size_t lastMatchPos, const size_t thisMatchPos, const std::string& rawSequence, const size_t k)
{
	const __uint128_t mask = ((__uint128_t)1 << (__uint128_t)(2ull*k))-(__uint128_t)1;
	__uint128_t lastMatchKmer = 0;
	__uint128_t thisMatchKmer = 0;
	for (size_t i = 0; i < k; i++)
	{
		lastMatchKmer <<= 2;
		thisMatchKmer <<= 2;
		switch(rawSequence[lastMatchPos+i])
		{
			case 'A':
				lastMatchKmer += 0;
				break;
			case 'C':
				lastMatchKmer += 1;
				break;
			case 'G':
				lastMatchKmer += 2;
				break;
			case 'T':
				lastMatchKmer += 3;
				break;
			default:
				assert(false);
		}
		switch(rawSequence[thisMatchPos+i])
		{
			case 'A':
				thisMatchKmer += 0;
				break;
			case 'C':
				thisMatchKmer += 1;
				break;
			case 'G':
				thisMatchKmer += 2;
				break;
			case 'T':
				thisMatchKmer += 3;
				break;
			default:
				assert(false);
		}
	}
	if (lastMatchKmer == thisMatchKmer)
	{
		result += rawSequence.substr(lastMatchPos + k, thisMatchPos - lastMatchPos);
		return;
	}
	assert(keptKmers.count(lastMatchKmer) == 1);
	assert(keptKmers.count(thisMatchKmer) == 1);
	std::string unambiguousExtensionFromLast;
	__uint128_t wanderKmer = lastMatchKmer;
	phmap::flat_hash_set<__uint128_t> visited;
	while (wanderKmer != thisMatchKmer)
	{
		size_t uniqueSuccessor = 4;
		for (size_t i = 0; i < 4; i++)
		{
			__uint128_t testKmer = ((wanderKmer << (__uint128_t)2) & mask) + i;
			if (keptKmers.count(testKmer) == 0) continue;
			if (uniqueSuccessor == 4)
			{
				uniqueSuccessor = i;
			}
			else
			{
				uniqueSuccessor = 5;
			}
		}
		if (uniqueSuccessor >= 4)
		{
			break;
		}
		unambiguousExtensionFromLast += "ACGT"[uniqueSuccessor];
		size_t nextKmer = ((wanderKmer << 2) & mask) + uniqueSuccessor;
		wanderKmer = nextKmer;
		if (visited.count(wanderKmer) == 1) break;
		visited.insert(wanderKmer);
	}
	if (wanderKmer == thisMatchKmer)
	{
		assert(unambiguousExtensionFromLast.size() >= 1);
		size_t pathLength = visited.size();
		if (pathLength + 50 > thisMatchPos-lastMatchPos)
		{
			if (thisMatchPos-lastMatchPos+50 > pathLength)
			{
				result += unambiguousExtensionFromLast;
				return;
			}
		}
	}
	visited.clear();
	wanderKmer = thisMatchKmer;
	unambiguousExtensionFromLast.clear();
	while (wanderKmer != lastMatchKmer)
	{
		size_t uniquePredecessor = 4;
		for (size_t i = 0; i < 4; i++)
		{
			__uint128_t testKmer = ((wanderKmer >> 2)) + (i << (2ull*(k-1ull)));
			if (keptKmers.count(testKmer) == 0) continue;
			if (uniquePredecessor == 4)
			{
				uniquePredecessor = i;
			}
			else
			{
				uniquePredecessor = 5;
			}
		}
		if (uniquePredecessor >= 4)
		{
			break;
		}
		unambiguousExtensionFromLast += "ACGT"[wanderKmer & 3];
		size_t nextKmer = ((wanderKmer >> 2)) + (uniquePredecessor << (2ull*(k-1ull)));
		wanderKmer = nextKmer;
		if (visited.count(wanderKmer) == 1) break;
		visited.insert(wanderKmer);
	}
	if (wanderKmer == lastMatchKmer)
	{
		assert(unambiguousExtensionFromLast.size() >= 1);
		size_t pathLength = visited.size();
		if (pathLength + 20 > thisMatchPos-lastMatchPos)
		{
			if (thisMatchPos-lastMatchPos+20 > pathLength)
			{
				std::reverse(unambiguousExtensionFromLast.begin(), unambiguousExtensionFromLast.end());
				result += unambiguousExtensionFromLast;
				return;
			}
		}
	}
	result += rawSequence.substr(lastMatchPos + k, thisMatchPos - lastMatchPos);
	return;
}

void removeTips(phmap::flat_hash_set<__uint128_t>& keptKmers, const phmap::flat_hash_map<__uint128_t, size_t>& kmerCounts, const size_t k)
{
	phmap::flat_hash_set<__uint128_t> removeThese;
	for (__uint128_t kmer : keptKmers)
	{
		assert(kmerCounts.count(kmer) == 1);
		assert(kmerCounts.at(kmer) >= 2);
		if (kmerCounts.at(kmer) != 2) continue;
		if (isTipFw(kmer, keptKmers, kmerCounts, k))
		{
			checkCanRemoveBw(removeThese, kmer, keptKmers, kmerCounts, k);
		}
		if (isTipBw(kmer, keptKmers, kmerCounts, k))
		{
			checkCanRemoveFw(removeThese, kmer, keptKmers, kmerCounts, k);
		}
	}
	for (__uint128_t kmer : removeThese)
	{
		assert(keptKmers.count(kmer) == 1);
		keptKmers.erase(kmer);
	}
}

std::vector<std::string> getSelfCorrectedLoopsOneIteration(const std::vector<std::string>& cluster, const size_t k)
{
	phmap::flat_hash_map<__uint128_t, size_t> kmerCounts;
	for (size_t i = 0; i < cluster.size(); i++)
	{
		phmap::flat_hash_set<__uint128_t> presentKmersThisRead;
		iterateKmers(cluster[i], k, [&presentKmersThisRead](const __uint128_t kmer, const size_t position) { presentKmersThisRead.insert(kmer); });
		for (auto kmer : presentKmersThisRead)
		{
			kmerCounts[kmer] += 1;
		}
	}
	phmap::flat_hash_set<__uint128_t> keptKmers;
	for (auto pair : kmerCounts)
	{
		if (pair.second == 1) continue;
		keptKmers.insert(pair.first);
	}
	removeTips(keptKmers, kmerCounts, k);
	std::vector<std::string> result;
	for (size_t i = 0; i < cluster.size(); i++)
	{
		std::string resultHere;
		size_t lastKmerMatch = std::numeric_limits<size_t>::max();
		iterateKmers(cluster[i], k, [&keptKmers, &cluster, &resultHere, &lastKmerMatch, i, k](const __uint128_t kmer, const size_t pos)
		{
			if (keptKmers.count(kmer) == 0) return;
			if (lastKmerMatch == std::numeric_limits<size_t>::max())
			{
				resultHere += cluster[i].substr(0, pos+k);
				lastKmerMatch = pos;
				return;
			}
			if (pos == lastKmerMatch+1)
			{
				resultHere += cluster[i][lastKmerMatch+k];
				lastKmerMatch = pos;
				return;
			}
			addPath(resultHere, keptKmers, lastKmerMatch, pos, cluster[i], k);
			lastKmerMatch = pos;
		});
		if (lastKmerMatch == std::numeric_limits<size_t>::max())
		{
			resultHere = cluster[i];
		}
		else
		{
			resultHere += cluster[i].substr(lastKmerMatch+k);
		}
		result.emplace_back(resultHere);
	}
	return result;
}

std::vector<std::string> getSelfCorrectedLoops(const std::vector<OntLoop>& cluster)
{
	std::vector<std::string> raws;
	for (size_t i = 0; i < cluster.size(); i++)
	{
		raws.emplace_back(cluster[i].rawSequence);
	}
	std::vector<std::string> partwiseCorrected;
	partwiseCorrected.resize(cluster.size());
	for (size_t i = 0; i < 10; i++)
	{
		std::vector<std::string> parts;
		for (size_t j = 0; j < cluster.size(); j++)
		{
			parts.emplace_back(raws[j].substr(raws[j].size()*i/10, raws[j].size()*(i+1)/10 - raws[j].size()*i/10));
		}
		std::vector<std::string> correctedParts = getSelfCorrectedLoopsOneIteration(parts, 15);
		correctedParts = getSelfCorrectedLoopsOneIteration(correctedParts, 31);
		correctedParts = getSelfCorrectedLoopsOneIteration(correctedParts, 63);
		for (size_t j = 0; j < cluster.size(); j++)
		{
			partwiseCorrected[j] += correctedParts[j];
		}
	}
	phmap::flat_hash_map<uint64_t, size_t> kmerCounts;
	std::vector<std::string> fullCorrected = getSelfCorrectedLoopsOneIteration(partwiseCorrected, 31);
	fullCorrected = getSelfCorrectedLoopsOneIteration(partwiseCorrected, 63);
	return fullCorrected;
}

std::string getPolishedSequence(const std::string& rawConsensus, const std::vector<std::string>& loopSequences, const size_t numThreads)
{
	std::vector<std::map<std::tuple<size_t, size_t, std::string>, size_t>> editCountsPerThread;
	editCountsPerThread.resize(numThreads);
	iterateEdits(rawConsensus, loopSequences, numThreads, [&editCountsPerThread](const size_t threadIndex, const size_t readIndex, const size_t firstMatchRef, const size_t lastMatchRef, const size_t firstMatchRead, const size_t lastMatchRead, const std::tuple<size_t, size_t, std::string, size_t>& edit)
	{
		assert(threadIndex < editCountsPerThread.size());
		auto filterEdit = std::make_tuple(std::get<0>(edit), std::get<1>(edit), std::get<2>(edit));
		editCountsPerThread[threadIndex][filterEdit] += 1;
	});
	std::map<std::tuple<size_t, size_t, std::string>, size_t> editCounts;
	for (size_t i = 0; i < editCountsPerThread.size(); i++)
	{
		for (const auto& pair : editCountsPerThread[i])
		{
			editCounts[pair.first] += pair.second;
		}
	}
	std::vector<std::tuple<size_t, size_t, std::string>> pickedEdits;
	for (const auto& pair : editCounts)
	{
		if (pair.second*1.5 <= loopSequences.size()) continue; // only pick edits supported by strictly >2/3 reads
		pickedEdits.emplace_back(pair.first);
	}
	std::sort(pickedEdits.begin(), pickedEdits.end());
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "picked " << pickedEdits.size() << " edits" << std::endl;
	std::string result;
	size_t lastMatch = 0;
	for (size_t i = 0; i < pickedEdits.size(); i++)
	{
		assert(i == 0 || std::get<0>(pickedEdits[i]) >= std::get<1>(pickedEdits[i-1]));
		assert(std::get<0>(pickedEdits[i]) >= lastMatch);
		result += rawConsensus.substr(lastMatch, std::get<0>(pickedEdits[i]) - lastMatch);
		result += std::get<2>(pickedEdits[i]);
		lastMatch = std::get<1>(pickedEdits[i]);
	}
	result += rawConsensus.substr(lastMatch);
	return result;
}

std::string getCentralModalAllele(const std::vector<std::string>& sequences, const std::vector<size_t>& startIndices, const std::vector<size_t>& endIndices)
{
	assert(startIndices.size() == sequences.size());
	assert(endIndices.size() == sequences.size());
	phmap::flat_hash_map<std::string, size_t> alleleCounts;
	for (size_t i = 0; i < sequences.size(); i++)
	{
		assert(endIndices[i] > startIndices[i]);
		assert(endIndices[i] + 1 <= sequences[i].size());
		std::string alleleHere = sequences[i].substr(startIndices[i], endIndices[i]+1-startIndices[i]);
		if (!(alleleHere.size() >= 2))
		{
			std::cerr << sequences[i].size() << " " << startIndices[i] << " " << endIndices[i] << " " << alleleHere.size() << std::endl;
		}
		assert(alleleHere.size() >= 2);
		alleleCounts[alleleHere] += 1;
	}
	std::vector<std::string> modalAlleles;
	size_t modalAlleleCoverage = 0;
	for (const auto& pair : alleleCounts)
	{
		assert(pair.first.size() >= 2);
		if (pair.second < modalAlleleCoverage) continue;
		if (pair.second > modalAlleleCoverage)
		{
			modalAlleles.clear();
			modalAlleleCoverage = pair.second;
		}
		modalAlleles.emplace_back(pair.first);
	}
	assert(modalAlleles.size() >= 1);
	if (modalAlleles.size() == 1) return modalAlleles[0];
	if (modalAlleles.size() == sequences.size()) return modalAlleles[0];
	std::vector<size_t> weightedEdits;
	weightedEdits.resize(modalAlleles.size(), 0);
	for (size_t i = 0; i < modalAlleles.size(); i++)
	{
		for (const auto& pair : alleleCounts)
		{
			if (pair.second == modalAlleleCoverage) continue;
			if (pair.first == modalAlleles[i]) continue;
			size_t edits = getEditDistanceWfa(modalAlleles[i], pair.first);
			weightedEdits[i] += edits * pair.second;
		}
		for (size_t j = 0; j < i; j++)
		{
			size_t edits = getEditDistanceWfa(modalAlleles[i], modalAlleles[j]);
			weightedEdits[i] += edits * modalAlleleCoverage;
			weightedEdits[j] += edits * modalAlleleCoverage;
		}
	}
	size_t bestIndex = 0;
	for (size_t i = 1; i < modalAlleles.size(); i++)
	{
		if (weightedEdits[i] < weightedEdits[bestIndex]) bestIndex = i;
	}
	assert(bestIndex < modalAlleles.size());
	assert(modalAlleles[bestIndex].size() >= 2);
	return modalAlleles[bestIndex];
}

std::string polishByBubbles(const std::string& refSequence, const std::vector<std::string>& ontSequences)
{
	std::vector<std::string> sequencesWithRef;
	sequencesWithRef.emplace_back(refSequence);
	sequencesWithRef.insert(sequencesWithRef.end(), ontSequences.begin(), ontSequences.end());
	std::vector<std::vector<size_t>> kmerChain = getMatchBases(refSequence, sequencesWithRef);
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "bubble polishing has " << kmerChain.size() << " shared kmers" << std::endl;
	if (kmerChain.size() < 2) return refSequence;
	assert(kmerChain[0].size() == sequencesWithRef.size());
	std::string result;
	result = refSequence.substr(0, kmerChain[0][0]);
	for (size_t i = 1; i < kmerChain.size(); i++)
	{
		assert(kmerChain[i][0] > kmerChain[i-1][0]);
		if (kmerChain[i][0] > kmerChain[i-1][0] + 5000)
		{
			result += refSequence.substr(kmerChain[i-1][0], kmerChain[i][0]-kmerChain[i-1][0]);
			continue;
		}
		std::string alleleHere = getCentralModalAllele(sequencesWithRef, kmerChain[i-1], kmerChain[i]);
		assert(alleleHere.size() >= 2);
		result += alleleHere.substr(0, alleleHere.size()-1);
	}
	result += refSequence.substr(kmerChain.back()[0]);
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "got corrected sequence" << std::endl;
	return result;
}

std::string polishConsensus(const std::string rawConsensus, const std::vector<std::string>& seqs, const size_t numThreads)
{
	std::string result = polishByBubbles(rawConsensus, seqs);
	for (size_t j = 0; j < 5; j++)
	{
		std::string newCorrection = getPolishedSequence(result, seqs, numThreads);
		if (newCorrection == result) break;
		result = newCorrection;
	}
	return result;
}

phmap::flat_hash_map<uint64_t, size_t> getRefKmers(const std::string_view& ref, const size_t k)
{
	assert(k <= 31);
	const uint64_t mask = (1ull << (2ull * k)) - 1;
	phmap::flat_hash_map<uint64_t, size_t> result;
	assert(ref.size() >= k);
	uint64_t kmer = 0;
	for (size_t i = 0; i < k; i++)
	{
		kmer <<= 2;
		switch(ref[i])
		{
			case 'A':
			case 'a':
				kmer += 0;
				break;
			case 'C':
			case 'c':
				kmer += 1;
				break;
			case 'G':
			case 'g':
				kmer += 2;
				break;
			case 'T':
			case 't':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
	}
	result[kmer] = 0;
	phmap::flat_hash_set<uint64_t> eraseThese;
	for (size_t i = k; i < ref.size(); i++)
	{
		kmer <<= 2;
		kmer &= mask;
		switch(ref[i])
		{
			case 'A':
			case 'a':
				kmer += 0;
				break;
			case 'C':
			case 'c':
				kmer += 1;
				break;
			case 'G':
			case 'g':
				kmer += 2;
				break;
			case 'T':
			case 't':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
		if (result.count(kmer) == 1) eraseThese.emplace(kmer);
		result[kmer] = i-k+1;
	}
	for (uint64_t kmer : eraseThese)
	{
		assert(result.count(kmer) == 1);
		result.erase(kmer);
	}
	return result;
}

std::vector<std::pair<size_t, size_t>> getKmerAnchors(const std::string_view& ref, const phmap::flat_hash_map<uint64_t, size_t>& refKmers, const std::string_view& query, const size_t k)
{
	assert(k <= 31);
	const uint64_t mask = (1ull << (2ull * k)) - 1;
	assert(ref.size() >= k);
	uint64_t kmer = 0;
	std::vector<std::pair<size_t, size_t>> kmerMatches;
	for (size_t i = 0; i < k; i++)
	{
		kmer <<= 2;
		switch(query[i])
		{
			case 'A':
			case 'a':
				kmer += 0;
				break;
			case 'C':
			case 'c':
				kmer += 1;
				break;
			case 'G':
			case 'g':
				kmer += 2;
				break;
			case 'T':
			case 't':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
	}
	if (refKmers.count(kmer) == 1) kmerMatches.emplace_back(refKmers.at(kmer), 0);
	for (size_t i = k; i < query.size(); i++)
	{
		kmer <<= 2;
		kmer &= mask;
		switch(query[i])
		{
			case 'A':
			case 'a':
				kmer += 0;
				break;
			case 'C':
			case 'c':
				kmer += 1;
				break;
			case 'G':
			case 'g':
				kmer += 2;
				break;
			case 'T':
			case 't':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
		if (refKmers.count(kmer) == 1) kmerMatches.emplace_back(refKmers.at(kmer), i-k+1);
	}
	for (auto pair : kmerMatches)
	{
		assert(ref.substr(pair.first, k) == query.substr(pair.second, k));
	}
	std::sort(kmerMatches.begin(), kmerMatches.end());
	removeNonDiagonalKmers(kmerMatches, k);
	removeNonDiagonalKmers(kmerMatches);
	removeRepeatKmers(kmerMatches);
	std::vector<std::pair<size_t, size_t>> chain = getIncreasingChain(kmerMatches);
	return chain;
}

std::vector<Node> getConsensusPath(const std::vector<OntLoop>& rawPaths, const GfaGraph& graph)
{
	assert(rawPaths.size() >= 1);
	auto coreNodes = getCoreNodes(rawPaths);
	if (coreNodes.size() == 0)
	{
		coreNodes = getMajorityNodes(rawPaths);
	}
	assert(coreNodes.size() >= 1);
	std::unordered_map<std::vector<Node>, size_t> alleleCounts;
	Node fakeStart { graph.numNodes()+1, true };
	Node fakeEnd { graph.numNodes()+2, true };
	for (const auto& read : rawPaths)
	{
		const auto& path = read.path;
		size_t lastCore = std::numeric_limits<size_t>::max();
		for (size_t i = 0; i < path.size(); i++)
		{
			assert(path[i].id() != fakeStart.id());
			assert(path[i].id() != fakeEnd.id());
			if (coreNodes.count(path[i].id()) == 0) continue;
			std::vector<Node> subpath { path.begin() + (lastCore == std::numeric_limits<size_t>::max() ? 0 : lastCore), path.begin() + i + 1 };
			assert(subpath.size() >= 2 || (i == 0 && subpath.size() == 1));
			if (lastCore == std::numeric_limits<size_t>::max()) subpath.insert(subpath.begin(), fakeStart);
			alleleCounts[subpath] += 1;
			lastCore = i;
		}
		std::vector<Node> subpath { path.begin() + (lastCore == std::numeric_limits<size_t>::max() ? 0 : lastCore), path.end()};
		assert(subpath.size() >= 2 || (lastCore == path.size()-1 && subpath.size() == 1));
		if (lastCore == std::numeric_limits<size_t>::max()) subpath.insert(subpath.begin(), fakeStart);
		subpath.insert(subpath.end(), fakeEnd);
		alleleCounts[subpath] += 1;
	}
	auto corePath = getMajorityPath(coreNodes, rawPaths);
	std::unordered_map<size_t, size_t> coreNodePositionInPath;
	for (size_t i = 0; i < corePath.size(); i++)
	{
		coreNodePositionInPath[corePath[i]] = i+1;
	}
	coreNodePositionInPath[fakeStart.id()] = 0;
	coreNodePositionInPath[fakeEnd.id()] = corePath.size()+1;
	std::vector<std::pair<std::vector<Node>, size_t>> bestAlleles;
	bestAlleles.resize(corePath.size()+1, std::make_pair(std::vector<Node>{}, 0));
	for (auto pair : alleleCounts)
	{
		if (coreNodePositionInPath.count(pair.first[0].id()) == 0) continue;
		if (coreNodePositionInPath.count(pair.first.back().id()) == 0) continue;
		if (coreNodePositionInPath.at(pair.first.back().id()) != coreNodePositionInPath.at(pair.first[0].id()) + 1) continue;
		size_t index = coreNodePositionInPath.at(pair.first[0].id());
		if (pair.second > bestAlleles[index].second) bestAlleles[index] = pair;
	}
	std::vector<Node> consensusPath;
	for (size_t i = 0; i < bestAlleles.size(); i++)
	{
		assert(bestAlleles[i].second >= 1);
		assert(bestAlleles[i].first.size() >= 1);
		if (i == 0)
		{
			consensusPath.insert(consensusPath.end(), bestAlleles[i].first.begin(), bestAlleles[i].first.end());
		}
		else
		{
			assert(consensusPath.size() >= 1);
			assert(consensusPath.back() == bestAlleles[i].first[0]);
			consensusPath.insert(consensusPath.end(), bestAlleles[i].first.begin()+1, bestAlleles[i].first.end());
		}
	}
	assert(consensusPath.size() >= 3);
	assert(consensusPath[0] == fakeStart);
	assert(consensusPath.back() == fakeEnd);
	consensusPath.erase(consensusPath.begin());
	consensusPath.erase(consensusPath.begin()+consensusPath.size()-1);
	assert(consensusPath.size() >= 1);
	return consensusPath;
}

std::string getConsensusSequence(const std::vector<Node>& consensusPath, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
{
	std::string consensusSeq = getSequence(consensusPath, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
	consensusSeq = consensusSeq.substr(pathStartClip.at(consensusPath[0]), consensusSeq.size() - pathStartClip.at(consensusPath[0]) - pathEndClip.at(consensusPath.back()));
	return consensusSeq;
}

std::vector<std::tuple<size_t, size_t, size_t>> getMatchRegionsRec(const std::string_view& refSequence, const std::string_view& querySequence, const size_t k)
{
	if (refSequence.size() < k) return std::vector<std::tuple<size_t, size_t, size_t>> {};
	if (querySequence.size() < k) return std::vector<std::tuple<size_t, size_t, size_t>> {};
	if (refSequence.size() == querySequence.size())
	{
		if (refSequence == querySequence)
		{
			std::vector<std::tuple<size_t, size_t, size_t>> result;
			result.emplace_back(0, 0, refSequence.size());
			return result;
		}
	}
	auto refKmers = getRefKmers(refSequence, k);
	auto kmerMatches = getKmerAnchors(refSequence, refKmers, querySequence, k);
	for (size_t i = 0; i < kmerMatches.size(); i++)
	{
		assert(refSequence.substr(kmerMatches[i].first, k) == querySequence.substr(kmerMatches[i].second, k));
	}
	for (size_t i = 1; i < kmerMatches.size(); i++)
	{
		assert(kmerMatches[i].first > kmerMatches[i-1].first);
		assert(kmerMatches[i].second > kmerMatches[i-1].second);
	}
	std::vector<std::tuple<size_t, size_t, size_t>> result;
	size_t lastRefPos = 0;
	size_t lastQueryPos = 0;
	for (size_t i = 0; i < kmerMatches.size(); i++)
	{
		if (i == 0 || (int)kmerMatches[i].first-(int)kmerMatches[i].second != (int)kmerMatches[i-1].first-(int)kmerMatches[i-1].second || kmerMatches[i].first > kmerMatches[i-1].first+k)
		{
			size_t refOffset = lastRefPos+1;
			size_t queryOffset = lastQueryPos+1;
			auto recResult = getMatchRegionsRec(std::string_view { refSequence.begin()+lastRefPos+1, kmerMatches[i].first+k-2 - lastRefPos }, std::string_view { querySequence.begin()+lastQueryPos+1, kmerMatches[i].second+k-2 - lastQueryPos }, k);
			if (recResult.size() == 0 && k > 11)
			{
				refOffset = lastRefPos;
				queryOffset = lastQueryPos;
				recResult = getMatchRegionsRec(std::string_view { refSequence.begin()+lastRefPos, kmerMatches[i].first+k - lastRefPos }, std::string_view { querySequence.begin()+lastQueryPos, kmerMatches[i].second+k - lastQueryPos }, k-10);
			}
			for (size_t j = 1; j < recResult.size(); j++)
			{
				if (!(std::get<0>(recResult[j]) > std::get<0>(recResult[j-1]) || (std::get<0>(recResult[j]) == std::get<0>(recResult[j-1]) && std::get<1>(recResult[j]) == std::get<1>(recResult[j-1]) && std::get<2>(recResult[j]) != std::get<2>(recResult[j-1]))))
				{
					for (size_t m = 0; m < recResult.size(); m++)
					{
						std::cerr << m << " " << std::get<0>(recResult[m]) << " " << std::get<1>(recResult[m]) << " " << std::get<2>(recResult[m]) << std::endl;
					}
					std::cerr << j << std::endl;
				}
				assert(std::get<0>(recResult[j]) > std::get<0>(recResult[j-1]) || (std::get<0>(recResult[j]) == std::get<0>(recResult[j-1]) && std::get<1>(recResult[j]) == std::get<1>(recResult[j-1]) && std::get<2>(recResult[j]) != std::get<2>(recResult[j-1])));
				assert(std::get<1>(recResult[j]) > std::get<1>(recResult[j-1]) || (std::get<0>(recResult[j]) == std::get<0>(recResult[j-1]) && std::get<1>(recResult[j]) == std::get<1>(recResult[j-1]) && std::get<2>(recResult[j]) != std::get<2>(recResult[j-1])));
			}
			for (auto t : recResult)
			{
				if (std::get<0>(t)+refOffset >= kmerMatches[i].first) continue;
				if (std::get<1>(t)+queryOffset >= kmerMatches[i].second) continue;
				if (i > 0 && std::get<0>(t)+std::get<2>(t)+refOffset <= kmerMatches[i-1].first+k) continue;
				if (i > 0 && std::get<1>(t)+std::get<2>(t)+queryOffset <= kmerMatches[i-1].second+k) continue;
				result.emplace_back(std::get<0>(t)+refOffset, std::get<1>(t)+queryOffset, std::get<2>(t));
				assert(refSequence.substr(std::get<0>(result.back()), std::get<2>(result.back())) == querySequence.substr(std::get<1>(result.back()), std::get<2>(result.back())));
			}
		}
		lastRefPos = kmerMatches[i].first;
		lastQueryPos = kmerMatches[i].second;
		result.emplace_back(kmerMatches[i].first, kmerMatches[i].second, k);
	}
	std::vector<std::tuple<size_t, size_t, size_t>> recResult;
	size_t refOffset = lastRefPos;
	size_t queryOffset = lastQueryPos;
	if (kmerMatches.size() >= 1)
	{
		refOffset = lastRefPos+1;
		queryOffset = lastQueryPos+1;
		recResult = getMatchRegionsRec(std::string_view { refSequence.begin()+lastRefPos+1, refSequence.size() - lastRefPos - 1 }, std::string_view { querySequence.begin()+lastQueryPos + 1, querySequence.size() - lastQueryPos-1 }, k);
	}
	if (recResult.size() == 0 && k > 11)
	{
		refOffset = lastRefPos;
		queryOffset = lastQueryPos;
		recResult = getMatchRegionsRec(std::string_view { refSequence.begin()+lastRefPos, refSequence.size() - lastRefPos }, std::string_view { querySequence.begin()+lastQueryPos, querySequence.size() - lastQueryPos }, k-10);
	}
	for (auto t : recResult)
	{
		if (kmerMatches.size() > 0 && std::get<0>(t)+std::get<2>(t)+refOffset <= kmerMatches.back().first+k) continue;
		if (kmerMatches.size() > 0 && std::get<1>(t)+std::get<2>(t)+queryOffset <= kmerMatches.back().second+k) continue;
		result.emplace_back(std::get<0>(t)+refOffset, std::get<1>(t)+queryOffset, std::get<2>(t));
		assert(refSequence.substr(std::get<0>(result.back()), std::get<2>(result.back())) == querySequence.substr(std::get<1>(result.back()), std::get<2>(result.back())));
	}
	return result;
}

std::vector<std::tuple<size_t, size_t, size_t>> getMatchRegions(const std::string& refSequence, const std::string& querySequence)
{
	auto kmerMatches = getMatchRegionsRec(refSequence, querySequence, 31);
	removeNonDiagonalKmers(kmerMatches);
	if (kmerMatches.size() == 0) return kmerMatches;
	std::sort(kmerMatches.begin(), kmerMatches.end());
	for (size_t j = 1; j < kmerMatches.size(); j++)
	{
		assert(std::get<0>(kmerMatches[j]) > std::get<0>(kmerMatches[j-1]));
		assert(std::get<1>(kmerMatches[j]) > std::get<1>(kmerMatches[j-1]));
	}
	for (size_t j = 0; j < kmerMatches.size(); j++)
	{
		assert(refSequence.substr(std::get<0>(kmerMatches[j]), std::get<2>(kmerMatches[j])) == querySequence.substr(std::get<1>(kmerMatches[j]), std::get<2>(kmerMatches[j])));
	}
	std::sort(kmerMatches.begin(), kmerMatches.end(), [](auto left, auto right) {
		int leftDiagonal = (int)std::get<0>(left)-(int)std::get<1>(left);
		int rightDiagonal = (int)std::get<0>(right)-(int)std::get<1>(right);
		if (leftDiagonal < rightDiagonal) return true;
		if (leftDiagonal > rightDiagonal) return false;
		assert(leftDiagonal == rightDiagonal);
		return std::get<0>(left) < std::get<0>(right);
	});
	std::vector<std::tuple<size_t, size_t, size_t>> diagonals;
	for (auto t : kmerMatches)
	{
		if (diagonals.size() == 0)
		{
			diagonals.emplace_back(t);
			continue;
		}
		if ((int)std::get<0>(t)-(int)std::get<1>(t) != (int)std::get<0>(diagonals.back())-(int)std::get<1>(diagonals.back()))
		{
			diagonals.emplace_back(t);
			continue;
		}
		if (std::get<0>(t) > std::get<0>(diagonals.back()) + std::get<2>(diagonals.back()))
		{
			diagonals.emplace_back(t);
			continue;
		}
		assert(std::get<0>(t) >= std::get<0>(diagonals.back()));
		std::get<2>(diagonals.back()) = std::max(std::get<2>(diagonals.back()), std::get<0>(t) + std::get<2>(t) - std::get<0>(diagonals.back()));
	}
	assert(diagonals.size() >= 1);
	std::vector<bool> contained;
	contained.resize(diagonals.size(), false);
	std::sort(diagonals.begin(), diagonals.end(), [](auto left, auto right)
	{
		if (std::get<1>(left) < std::get<1>(right)) return true;
		if (std::get<1>(left) > std::get<1>(right)) return false;
		if (std::get<2>(left) > std::get<2>(right)) return true;
		if (std::get<2>(left) < std::get<2>(right)) return false;
		return std::get<0>(left) < std::get<0>(right);
	});
	size_t lastMatchEnd = std::get<1>(diagonals[0])+std::get<2>(diagonals[0]);
	for (size_t i = 1; i < diagonals.size(); i++)
	{
		size_t endHere = std::get<1>(diagonals[i]) + std::get<2>(diagonals[i]);
		if (endHere <= lastMatchEnd) contained[i] = true;
		lastMatchEnd = std::max(lastMatchEnd, endHere);
	}
	for (size_t i = diagonals.size()-1; i < diagonals.size(); i--)
	{
		if (!contained[i]) continue;
		std::swap(diagonals[i], diagonals.back());
		diagonals.pop_back();
	}
	assert(diagonals.size() >= 1);
	std::sort(diagonals.begin(), diagonals.end(), [](auto left, auto right)
	{
		if (std::get<0>(left) < std::get<0>(right)) return true;
		if (std::get<0>(left) > std::get<0>(right)) return false;
		if (std::get<2>(left) > std::get<2>(right)) return true;
		if (std::get<2>(left) < std::get<2>(right)) return false;
		return std::get<1>(left) < std::get<1>(right);
	});
	contained.resize(0, false);
	contained.resize(diagonals.size(), false);
	lastMatchEnd = std::get<0>(diagonals[0])+std::get<2>(diagonals[0]);
	for (size_t i = 1; i < diagonals.size(); i++)
	{
		size_t endHere = std::get<0>(diagonals[i]) + std::get<2>(diagonals[i]);
		if (endHere <= lastMatchEnd) contained[i] = true;
		lastMatchEnd = std::max(lastMatchEnd, endHere);
	}
	for (size_t i = diagonals.size()-1; i < diagonals.size(); i--)
	{
		if (!contained[i]) continue;
		std::swap(diagonals[i], diagonals.back());
		diagonals.pop_back();
	}
	std::sort(diagonals.begin(), diagonals.end());
	assert(diagonals.size() >= 1);
	for (size_t j = 0; j < diagonals.size(); j++)
	{
		assert(refSequence.substr(std::get<0>(diagonals[j]), std::get<2>(diagonals[j])) == querySequence.substr(std::get<1>(diagonals[j]), std::get<2>(diagonals[j])));
	}
	return diagonals;
}

std::vector<std::vector<size_t>> getMatchBases(const std::string& consensusSeq, const std::vector<std::string>& loopSequences)
{
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> matchRegions;
	std::vector<std::pair<size_t, size_t>> refAllowedAreas;
	std::tie(refAllowedAreas, matchRegions) = getRegionsConservedInAllSequences(consensusSeq, loopSequences);
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "size " << loopSequences.size() << " ref length " << consensusSeq.size() << " all-present ref regions:";
	for (size_t i = 0; i < refAllowedAreas.size(); i++)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << refAllowedAreas[i].first << "-" << refAllowedAreas[i].second;
	}
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
	if (refAllowedAreas.size() == 0) return std::vector<std::vector<size_t>> {};
	for (size_t i = 0; i < refAllowedAreas.size(); i++)
	{
		assert(refAllowedAreas[i].second > refAllowedAreas[i].first);
		assert(i == 0 || refAllowedAreas[i].first > refAllowedAreas[i-1].second);
	}
	std::vector<std::vector<size_t>> matchesPerRead;
	matchesPerRead.resize(matchRegions.size());
	for (size_t i = 0; i < matchRegions.size(); i++)
	{
		std::sort(matchRegions[i].begin(), matchRegions[i].end(), [](auto left, auto right)
		{
			if (std::get<0>(left) < std::get<0>(right)) return true;
			if (std::get<0>(left) > std::get<0>(right)) return false;
			if (std::get<2>(left) > std::get<2>(right)) return true;
			if (std::get<2>(left) < std::get<2>(right)) return false;
			return std::get<1>(left) < std::get<1>(right);
		});
		size_t readIndex = 0;
		size_t refIndex = 0;
		while (refIndex < refAllowedAreas.size() && readIndex < matchRegions[i].size())
		{
			if (std::get<0>(matchRegions[i][readIndex])+std::get<2>(matchRegions[i][readIndex]) <= refAllowedAreas[refIndex].first)
			{
				readIndex += 1;
				continue;
			}
			if (refAllowedAreas[refIndex].second <= std::get<0>(matchRegions[i][readIndex]))
			{
				refIndex += 1;
				continue;
			}
			assert(refAllowedAreas[refIndex].second > std::get<0>(matchRegions[i][readIndex]));
			assert(std::get<0>(matchRegions[i][readIndex]) + std::get<2>(matchRegions[i][readIndex]) > refAllowedAreas[refIndex].first);
			size_t intersectRefStart = std::max(refAllowedAreas[refIndex].first, std::get<0>(matchRegions[i][readIndex]));
			size_t intersectRefEnd = std::min(refAllowedAreas[refIndex].second, std::get<0>(matchRegions[i][readIndex]) + std::get<2>(matchRegions[i][readIndex]));
			assert(intersectRefEnd > intersectRefStart);
			for (size_t j = 0; j < intersectRefEnd - intersectRefStart; j++)
			{
				size_t matchPos = std::get<1>(matchRegions[i][readIndex]) + (intersectRefStart - std::get<0>(matchRegions[i][readIndex])) + j;
				if (!(matchesPerRead[i].size() == 0 || matchPos > matchesPerRead[i].back()))
				{
					std::cerr << matchPos << " " << matchesPerRead[i].back() << std::endl;
				}
				assert(matchesPerRead[i].size() == 0 || matchPos > matchesPerRead[i].back());
				matchesPerRead[i].emplace_back(matchPos);
			}
			if (std::get<0>(matchRegions[i][readIndex]) + std::get<2>(matchRegions[i][readIndex]) < refAllowedAreas[refIndex].second)
			{
				readIndex += 1;
			}
			else
			{
				refIndex += 1;
			}
		}
	}
	for (size_t i = 1; i < matchesPerRead.size(); i++)
	{
		assert(matchesPerRead[i].size() == matchesPerRead[0].size());
	}
	std::vector<std::vector<size_t>> result;
	result.resize(matchesPerRead[0].size());
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i].resize(loopSequences.size());
		for (size_t j = 0; j < loopSequences.size(); j++)
		{
			result[i][j] = matchesPerRead[j][i];
		}
	}
	return result;
}

std::vector<std::pair<size_t, size_t>> mergeSpans(const std::vector<std::pair<size_t, size_t>>& raws)
{
	std::vector<std::pair<size_t, size_t>> result;
	for (auto pair : raws)
	{
		if (pair.second == pair.first) continue;
		if (result.size() == 0)
		{
			result.emplace_back(pair.first, pair.second);
			continue;
		}
		assert(pair.first >= result.back().first);
		if (pair.first <= result.back().second)
		{
			result.back().second = std::max(result.back().second, pair.second);
			continue;
		}
		result.emplace_back(pair);
	}
	return result;
}

std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>> getRegionsConservedInAllSequences(const std::string& consensusSeq, const std::vector<std::string>& loopSequences)
{
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> matchRegions;
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		matchRegions.emplace_back(getMatchRegions(consensusSeq, loopSequences[i]));
		if (matchRegions.back().size() == 0) return std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>>{};
	}
	std::vector<std::pair<size_t, size_t>> refForbiddenAreas;
	for (size_t i = 0; i < matchRegions.size(); i++)
	{
		std::sort(matchRegions[i].begin(), matchRegions[i].end(), [](auto left, auto right)
		{
			if (std::get<0>(left) < std::get<0>(right)) return true;
			if (std::get<0>(left) > std::get<0>(right)) return false;
			if (std::get<2>(left) > std::get<2>(right)) return true;
			if (std::get<2>(left) < std::get<2>(right)) return false;
			return std::get<1>(left) < std::get<1>(right);
		});
		assert(matchRegions[i].size() >= 1);
		refForbiddenAreas.emplace_back(0, std::get<0>(matchRegions[i][0]));
		refForbiddenAreas.emplace_back(std::get<0>(matchRegions[i].back()) + std::get<2>(matchRegions[i].back()), consensusSeq.size());
		for (size_t j = 1; j < matchRegions[i].size(); j++)
		{
			size_t startHere = std::get<0>(matchRegions[i][j]);
			size_t prevEnd = std::get<0>(matchRegions[i][j-1]) + std::get<2>(matchRegions[i][j-1]); 
			if (startHere > prevEnd)
			{
				refForbiddenAreas.emplace_back(prevEnd, startHere);
			}
			else
			{
				refForbiddenAreas.emplace_back(startHere, prevEnd);
			}
		}
		std::sort(matchRegions[i].begin(), matchRegions[i].end(), [](auto left, auto right)
		{
			if (std::get<1>(left) < std::get<1>(right)) return true;
			if (std::get<1>(left) > std::get<1>(right)) return false;
			if (std::get<2>(left) > std::get<2>(right)) return true;
			if (std::get<2>(left) < std::get<2>(right)) return false;
			return std::get<0>(left) < std::get<0>(right);
		});
		for (size_t j = 1; j < matchRegions[i].size(); j++)
		{
			size_t startHere = std::get<1>(matchRegions[i][j]);
			size_t prevEnd = std::get<1>(matchRegions[i][j-1]) + std::get<2>(matchRegions[i][j-1]); 
			if (startHere >= prevEnd) continue;
			assert(std::get<1>(matchRegions[i][j]) > std::get<1>(matchRegions[i][j-1]));
			size_t prevEndRefOffset = std::get<0>(matchRegions[i][j-1]) + std::get<2>(matchRegions[i][j-1]);
			size_t startHereRefOffset = std::get<1>(matchRegions[i][j]) - std::get<1>(matchRegions[i][j-1]) + std::get<0>(matchRegions[i][j-1]);
			assert(prevEndRefOffset > startHereRefOffset);
			refForbiddenAreas.emplace_back(startHereRefOffset, prevEndRefOffset);
		}
	}
	std::sort(refForbiddenAreas.begin(), refForbiddenAreas.end());
	refForbiddenAreas = mergeSpans(refForbiddenAreas);
	std::vector<std::pair<size_t, size_t>> refAllowedAreas;
	size_t lastForbiddenEnd = 0;
	for (auto pair : refForbiddenAreas)
	{
		if (pair.first > lastForbiddenEnd) refAllowedAreas.emplace_back(lastForbiddenEnd, pair.first);
		assert(pair.second > lastForbiddenEnd);
		lastForbiddenEnd = pair.second;
	}
	if (lastForbiddenEnd < consensusSeq.size()) refAllowedAreas.emplace_back(lastForbiddenEnd, consensusSeq.size());
	return std::make_pair(std::move(refAllowedAreas), std::move(matchRegions));
}
