#include <unordered_set>
#include "SequenceAligner.h"
#include "RibotinUtils.h"

size_t getEditDistancePossiblyMemoized(const std::vector<Node>& left, const std::vector<Node>& right, const size_t leftStartClipBp, const size_t rightStartClipBp, const size_t leftEndClipBp, const size_t rightEndClipBp, const GfaGraph& graph, const size_t maxEdits, std::unordered_map<std::pair<std::vector<Node>, std::vector<Node>>, size_t>& memoizedEditDistances, std::mutex& memoizationMutex)
{
	if (left == right) return 0;
	auto key = canon(left, right);
	size_t add;
	if (leftStartClipBp == 0 && rightStartClipBp == 0 && leftEndClipBp == 0 && rightEndClipBp == 0)
	{
		std::lock_guard<std::mutex> lock { memoizationMutex };
		if (memoizedEditDistances.count(key) == 1) return memoizedEditDistances.at(key);
	}
	auto leftSubseq = getSequenceView(left, graph.nodeSeqsTwobit, graph.revCompNodeSeqsTwobit, graph.edges, leftStartClipBp, leftEndClipBp);
	auto rightSubseq = getSequenceView(right, graph.nodeSeqsTwobit, graph.revCompNodeSeqsTwobit, graph.edges, rightStartClipBp, rightEndClipBp);
	add = getEditDistanceWfa(leftSubseq, rightSubseq, maxEdits);
	assert(add == getEditDistanceWfa(rightSubseq, leftSubseq, maxEdits));
	if (leftStartClipBp == 0 && rightStartClipBp == 0 && leftEndClipBp == 0 && rightEndClipBp == 0)
	{
		std::lock_guard<std::mutex> lock { memoizationMutex };
		if (memoizedEditDistances.count(key) == 0)
		{
			memoizedEditDistances[key] = add;
		}
		else
		{
			assert(memoizedEditDistances.at(key) == add);
		}
	}
	return add;
}

__attribute__((always_inline)) inline size_t getBitvectorMatchLength(uint64_t leftFirst, uint64_t leftSecond, uint64_t rightFirst, uint64_t rightSecond)
{
	uint64_t mismatch = (leftFirst ^ rightFirst) | (leftSecond ^ rightSecond);
	uint64_t prefixZerosPlusOneOne = (mismatch-1) ^ mismatch;
	uint64_t prefixMatchPlusOne = popcount(prefixZerosPlusOneOne);
	assert(prefixMatchPlusOne >= 1);
	return prefixMatchPlusOne-1;
}

__attribute__((always_inline)) inline size_t getMatchLength(const std::string& left, const std::string& right, const size_t leftStart, const size_t rightStart)
{
	if (leftStart == left.size() || rightStart == right.size()) return 0;
	size_t result = 0;
	while (leftStart+result < left.size() && rightStart+result < right.size())
	{
		if (left[leftStart+result] != right[rightStart+result]) break;
		result += 1;
	}
	return result;
}

__attribute__((always_inline)) inline size_t getMatchLength(const TwobitString& left, const TwobitString& right, const size_t leftStart, const size_t rightStart)
{
	size_t result = 0;
	while (true)
	{
		auto leftBits = left.getBitvectorsPossiblyTruncated(leftStart+result);
		auto rightBits = right.getBitvectorsPossiblyTruncated(rightStart+result);
		auto matchHere = getBitvectorMatchLength(leftBits.first, leftBits.second, rightBits.first, rightBits.second);
		size_t maxMatch = std::min(64 - ((leftStart+result) % 64), 64 - ((rightStart+result) % 64));
		matchHere = std::min(matchHere, maxMatch);
		if (matchHere == 0) return result;
		result += matchHere;
		if (leftStart+result >= left.size() || rightStart+result >= right.size())
		{
			result = std::min(result, left.size()-leftStart);
			result = std::min(result, right.size()-rightStart);
			return result;
		}
	}
}

size_t getMatchLength(const PathSequenceView& left, const PathSequenceView& right, const size_t leftStart, const size_t rightStart)
{
	size_t result = 0;
	size_t leftIndex = leftStart + left.leftClip;
	size_t rightIndex = rightStart + right.leftClip;
	size_t leftSize = left.size();
	size_t rightSize = right.size();
	if (leftIndex == leftSize + left.leftClip) return 0;
	if (rightIndex == rightSize + right.leftClip) return 0;
	assert(leftIndex < leftSize + left.leftClip);
	assert(rightIndex < rightSize + right.leftClip);
	size_t leftStartNode = 0;
	size_t rightStartNode = 0;
	while (true)
	{
		while (leftStartNode+1 < left.nodeStartPoses.size() && left.nodeStartPoses[leftStartNode+1] <= leftIndex) leftStartNode += 1;
		while (rightStartNode+1 < right.nodeStartPoses.size() && right.nodeStartPoses[rightStartNode+1] <= rightIndex) rightStartNode += 1;
		size_t leftOffset = leftIndex - left.nodeStartPoses[leftStartNode] + left.nodeStartOverlap[leftStartNode];
		size_t rightOffset = rightIndex - right.nodeStartPoses[rightStartNode] + right.nodeStartOverlap[rightStartNode];
		auto matchHere = getMatchLength(*left.nodePointers[leftStartNode], *right.nodePointers[rightStartNode], leftOffset, rightOffset);
		if (matchHere == 0) return result;
		result += matchHere;
		if (result >= leftSize - leftStart || result >= rightSize - rightStart)
		{
			result = std::min(result, leftSize - leftStart);
			result = std::min(result, rightSize - rightStart);
			return result;
		}
		leftIndex += matchHere;
		rightIndex += matchHere;
	}
}

size_t getEditDistanceWfa(const std::string& left, const std::string& right)
{
	std::vector<size_t> offsets;
	std::vector<size_t> nextOffsets;
	offsets.resize(1, 0);
	const size_t leftSize = left.size();
	const size_t rightSize = right.size();
	offsets[0] = getMatchLength(left, right, 0, 0);
	if (offsets[0] == leftSize) return rightSize - leftSize;
	if (offsets[0] == rightSize) return leftSize - rightSize;
	size_t zeroPos = 1;
	for (size_t score = 0; score < std::max(left.size(), right.size()); score++)
	{
		nextOffsets.resize(offsets.size()+2);
		for (size_t i = 0; i < nextOffsets.size(); i++)
		{
			nextOffsets[i] = 0;
			if (i < zeroPos && zeroPos-i >= left.size()) continue;
			if (i >= zeroPos + right.size()) break;
			if (i >= 1 && i < nextOffsets.size()-1) nextOffsets[i] = offsets[i-1]+1;
			if (i >= 2) nextOffsets[i] = std::max(nextOffsets[i], offsets[i-2]);
			if (i < nextOffsets.size()-2) nextOffsets[i] = std::max(nextOffsets[i], offsets[i]+1);
			assert(nextOffsets[i] + i >= zeroPos);
			nextOffsets[i] = std::min(nextOffsets[i], rightSize - i + zeroPos);
			nextOffsets[i] = std::min(nextOffsets[i], leftSize);
			assert(nextOffsets[i] + i >= zeroPos);
			assert(nextOffsets[i]+i-zeroPos <= rightSize);
			assert(nextOffsets[i] <= leftSize);
			nextOffsets[i] += getMatchLength(left, right, nextOffsets[i], nextOffsets[i]+i-zeroPos);
			assert(nextOffsets[i]+i-zeroPos <= rightSize);
			assert(nextOffsets[i] <= leftSize);
			if (nextOffsets[i] == leftSize && nextOffsets[i]+i-zeroPos == rightSize) return score;
			assert(nextOffsets[i] + i >= zeroPos);
		}
		std::swap(offsets, nextOffsets);
		zeroPos += 1;
		assert(score < leftSize+rightSize);
	}
	assert(false);
	return std::max(left.size(), right.size())+1;
}

size_t getEditDistanceWfa(const PathSequenceView& left, const PathSequenceView& right, const size_t maxEdits)
{
	if (left.size() > right.size()+maxEdits) return maxEdits+1;
	if (right.size() > left.size()+maxEdits) return maxEdits+1;
	std::vector<size_t> offsets;
	std::vector<size_t> nextOffsets;
	offsets.resize(1, 0);
	const size_t leftSize = left.size();
	const size_t rightSize = right.size();
	offsets[0] = getMatchLength(left, right, 0, 0);
	if (offsets[0] == leftSize) return rightSize - leftSize;
	if (offsets[0] == rightSize) return leftSize - rightSize;
	size_t zeroPos = 1;
	for (size_t score = 0; score < maxEdits; score++)
	{
		nextOffsets.resize(offsets.size()+2);
		for (size_t i = 0; i < nextOffsets.size(); i++)
		{
			nextOffsets[i] = 0;
			if (zeroPos > i && zeroPos-i > leftSize) continue;
			if (i >= zeroPos + rightSize) break;
			if (i >= 1 && i < nextOffsets.size()-1) nextOffsets[i] = offsets[i-1]+1;
			if (i >= 2) nextOffsets[i] = std::max(nextOffsets[i], offsets[i-2]);
			if (i < nextOffsets.size()-2) nextOffsets[i] = std::max(nextOffsets[i], offsets[i]+1);
			assert(nextOffsets[i] + i >= zeroPos);
			nextOffsets[i] = std::min(nextOffsets[i], rightSize - i + zeroPos);
			assert(nextOffsets[i]+i-zeroPos <= rightSize);
			nextOffsets[i] = std::min(nextOffsets[i], leftSize);
			assert(nextOffsets[i]+i-zeroPos <= rightSize);
			assert(nextOffsets[i] <= leftSize);
			nextOffsets[i] += getMatchLength(left, right, nextOffsets[i], nextOffsets[i]+i-zeroPos);
			if (nextOffsets[i] == leftSize && nextOffsets[i]+i-zeroPos == rightSize) return score;
		}
		std::swap(offsets, nextOffsets);
		zeroPos += 1;
		assert(score <= std::max(leftSize, rightSize));
	}
	return maxEdits+1;
}

std::vector<std::pair<size_t, size_t>> getNodeMatches(const std::vector<Node>& left, size_t leftIndex, size_t leftStart, size_t leftEnd, const std::vector<Node>& right, size_t rightIndex, size_t rightStart, size_t rightEnd, const std::vector<phmap::flat_hash_map<Node, size_t>>& nodeCountIndex, const std::vector<phmap::flat_hash_map<Node, size_t>>& nodePosIndex)
{
	std::vector<std::pair<size_t, size_t>> unfilteredMatches;
	if (leftEnd == leftStart) return unfilteredMatches;
	if (rightEnd == rightStart) return unfilteredMatches;
	assert(leftEnd > leftStart);
	assert(rightEnd > rightStart);
	assert(leftEnd <= left.size());
	assert(rightEnd <= right.size());
	bool allIncreasing = true;
	for (size_t i = leftStart; i < leftEnd; i++)
	{
		auto found = nodeCountIndex[rightIndex].find(left[i]);
		if (found == nodeCountIndex[rightIndex].end()) continue;
		if (found->second != 1) continue;
		found = nodeCountIndex[leftIndex].find(left[i]);
		if (found == nodeCountIndex[leftIndex].end()) continue;
		if (found->second != 1) continue;
		size_t rightPos = nodePosIndex[rightIndex].at(left[i]);
		if (rightPos < rightStart || rightPos >= rightEnd) continue;
		unfilteredMatches.emplace_back(i, nodePosIndex[rightIndex].at(left[i]));
		if (unfilteredMatches.size() >= 2 && unfilteredMatches.back().second <= unfilteredMatches[unfilteredMatches.size()-2].second) allIncreasing = false;
	}
	if (allIncreasing) return unfilteredMatches;
	std::vector<bool> matchInConflict;
	matchInConflict.resize(unfilteredMatches.size(), false);
	for (size_t i = 1; i < unfilteredMatches.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (unfilteredMatches[i].first <= unfilteredMatches[j].first || unfilteredMatches[i].second <= unfilteredMatches[j].second)
			{
				matchInConflict[i] = true;
				matchInConflict[j] = true;
			}
		}
	}
	std::vector<std::pair<size_t, size_t>> result;
	for (size_t i = 0; i < unfilteredMatches.size(); i++)
	{
		if (matchInConflict[i]) continue;
		result.emplace_back(unfilteredMatches[i]);
	}
	assert(result.size()+2 <= unfilteredMatches.size());
	return result;
}

size_t getEditDistance(const std::vector<Node>& left, const size_t leftIndex, const std::vector<Node>& right, const size_t rightIndex, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t maxEdits, const std::unordered_set<size_t>& coreNodes, const std::vector<phmap::flat_hash_map<Node, size_t>>& nodeCountIndex, const std::vector<phmap::flat_hash_map<Node, size_t>>& nodePosIndex, std::unordered_map<std::pair<std::vector<Node>, std::vector<Node>>, size_t>& memoizedEditDistances, std::mutex& memoizationMutex)
{
	std::vector<size_t> leftCoreMatchPositions;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (coreNodes.count(left[i].id()) == 0) continue;
		leftCoreMatchPositions.push_back(i);
	}
	std::vector<size_t> rightCoreMatchPositions;
	for (size_t i = 0; i < right.size(); i++)
	{
		if (coreNodes.count(right[i].id()) == 0) continue;
		rightCoreMatchPositions.push_back(i);
	}
	assert(leftCoreMatchPositions.size() == coreNodes.size());
	assert(rightCoreMatchPositions.size() == coreNodes.size());
	for (size_t i = 0; i < leftCoreMatchPositions.size(); i++)
	{
		assert(left[leftCoreMatchPositions[i]] == right[rightCoreMatchPositions[i]]);
	}
	std::vector<std::pair<size_t, size_t>> nodeMatches;
	size_t lastLeftStart = 0;
	size_t lastRightStart = 0;
	for (size_t i = 0; i < leftCoreMatchPositions.size(); i++)
	{
		auto addMatches = getNodeMatches(left, leftIndex, lastLeftStart, leftCoreMatchPositions[i], right, rightIndex, lastRightStart, rightCoreMatchPositions[i], nodeCountIndex, nodePosIndex);
		nodeMatches.insert(nodeMatches.end(), addMatches.begin(), addMatches.end());
		nodeMatches.emplace_back(leftCoreMatchPositions[i], rightCoreMatchPositions[i]);
		lastLeftStart = leftCoreMatchPositions[i]+1;
		lastRightStart = rightCoreMatchPositions[i]+1;
	}
	auto addMatches = getNodeMatches(left, leftIndex, lastLeftStart, left.size(), right, rightIndex, lastRightStart, right.size(), nodeCountIndex, nodePosIndex);
	nodeMatches.insert(nodeMatches.end(), addMatches.begin(), addMatches.end());
	assert(nodeMatches.size() >= coreNodes.size());
	if (nodeMatches.size() == 0)
	{
		return getEditDistancePossiblyMemoized(left, right, pathStartClip.at(left[0]), pathStartClip.at(right[0]), pathEndClip.at(left.back()), pathEndClip.at(right.back()), graph, maxEdits, memoizedEditDistances, memoizationMutex);
	}
	assert(nodeMatches.size() >= 1);
	size_t result = 0;
	size_t add = 0;
	for (size_t i = 1; i < nodeMatches.size(); i++)
	{
		assert(nodeMatches[i].first > nodeMatches[i-1].first);
		assert(nodeMatches[i].second > nodeMatches[i-1].second);
		if (nodeMatches[i].first == nodeMatches[i-1].first+1 && nodeMatches[i].second == nodeMatches[i-1].second+1) continue;
		std::vector<Node> leftPath { left.begin() + nodeMatches[i-1].first, left.begin()+nodeMatches[i].first+1 };
		std::vector<Node> rightPath { right.begin() + nodeMatches[i-1].second, right.begin()+nodeMatches[i].second+1 };
		add = getEditDistancePossiblyMemoized(leftPath, rightPath, 0, 0, 0, 0, graph, maxEdits, memoizedEditDistances, memoizationMutex);
		result += add;
		if (result >= maxEdits) return maxEdits+1;
	}
	std::vector<Node> leftPath { left.begin(), left.begin()+nodeMatches[0].first+1 };
	std::vector<Node> rightPath { right.begin(), right.begin()+nodeMatches[0].second+1 };
	add = getEditDistancePossiblyMemoized(leftPath, rightPath, pathStartClip.at(leftPath[0]), pathStartClip.at(rightPath[0]), 0, 0, graph, maxEdits, memoizedEditDistances, memoizationMutex);
	result += add;
	if (result >= maxEdits) return maxEdits+1;
	leftPath = std::vector<Node> { left.begin()+nodeMatches.back().first, left.end() };
	rightPath = std::vector<Node> { right.begin()+nodeMatches.back().second, right.end() };
	add = getEditDistancePossiblyMemoized(leftPath, rightPath, 0, 0, pathEndClip.at(leftPath.back()), pathEndClip.at(rightPath.back()), graph, maxEdits, memoizedEditDistances, memoizationMutex);
	result += add;
	if (result >= maxEdits) return maxEdits+1;
	return result;
}

size_t getEditDistance(const std::vector<Node>& left, const size_t leftIndex, const std::vector<Node>& right, const size_t rightIndex, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t maxEdits, const std::unordered_set<size_t>& coreNodes, const std::vector<phmap::flat_hash_map<Node, size_t>>& nodeCountIndex, const std::vector<phmap::flat_hash_map<Node, size_t>>& nodePosIndex, std::unordered_map<std::pair<std::vector<Node>, std::vector<Node>>, size_t>& memoizedEditDistances)
{
	std::mutex tmpMutex;
	return getEditDistance(left, leftIndex, right, rightIndex, graph, pathStartClip, pathEndClip, maxEdits, coreNodes, nodeCountIndex, nodePosIndex, memoizedEditDistances, tmpMutex);
}
