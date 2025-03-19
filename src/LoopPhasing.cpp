#include <set>
#include <fstream>
#include "LoopPhasing.h"
#include "RibotinUtils.h"
#include "FileHelper.h"
#include "ConsensusHelper.h"
#include "SequenceAligner.h"
#include "Logger.h"

bool isIndel(const std::tuple<size_t, size_t, std::string>& variant)
{
	if (std::get<1>(variant)-std::get<0>(variant) == std::get<2>(variant).size()) return false;
	return true;
}

bool isSNP(const std::tuple<size_t, size_t, std::string>& variant)
{
	if (std::get<1>(variant)-std::get<0>(variant) == 1 && std::get<2>(variant).size() == 1) return true;
	return false;
}

bool isRepeat(const std::string& seq, const size_t repeatlen)
{
	assert(seq.size() % repeatlen == 0);
	for (size_t i = 1; i*repeatlen < seq.size(); i++)
	{
		if (seq.substr(i*repeatlen, repeatlen) != seq.substr(0, repeatlen)) return false;
	}
	return true;
}

bool isMicrosatelliteOrHomopolymerIndel(const std::tuple<size_t, size_t, std::string>& variant, const std::string& refSequence)
{
	size_t deletionLength = std::get<1>(variant)-std::get<0>(variant);
	size_t insertionLength = std::get<2>(variant).size();
	if (deletionLength > 0 && insertionLength > 0) return false;
	for (size_t i = 1; i < 10; i++)
	{
		if (deletionLength >= i && deletionLength % i == 0)
		{
			if (isRepeat(refSequence.substr(std::get<0>(variant), std::get<1>(variant)-std::get<0>(variant)), i))
			{
				if (refSequence.substr(std::get<1>(variant), i) == refSequence.substr(std::get<0>(variant), i))
				{
					return true;
				}
				if (std::get<0>(variant) >= i && refSequence.substr(std::get<0>(variant)-i, i) == refSequence.substr(std::get<0>(variant), i))
				{
					return true;
				}
			}
		}
		if (insertionLength >= i && insertionLength % i == 0)
		{
			if (isRepeat(std::get<2>(variant), i))
			{
				if (refSequence.substr(std::get<1>(variant), i) == std::get<2>(variant).substr(0, i))
				{
					return true;
				}
				if (std::get<0>(variant) >= i && refSequence.substr(std::get<0>(variant)-i, i) == std::get<2>(variant).substr(0, i))
				{
					return true;
				}
			}
		}
	}
	return false;
}

bool isBigIndel(const std::tuple<size_t, size_t, std::string>& variant)
{
	if (std::get<1>(variant)-std::get<0>(variant) >= 10) return true;
	if (std::get<2>(variant).size() >= 10) return true;
	return false;
}

void clusterBubbleAlleles(std::vector<std::vector<size_t>>& result, const std::vector<std::string>& loopSequences, const std::vector<size_t>& previousMatchIndices, const std::vector<size_t>& currentMatchIndices)
{
	const size_t minCoverage = std::max<size_t>(5, loopSequences.size()*0.05);
	const double maxDivergence = 0.1;
	size_t maxClusterDistance = 20;
	if (currentMatchIndices[0] == previousMatchIndices[0] + 1)
	{
		bool allIdentical = true;
		for (size_t i = 1; i < loopSequences.size(); i++)
		{
			if (currentMatchIndices[i] - previousMatchIndices[i] != currentMatchIndices[0] - previousMatchIndices[0])
			{
				allIdentical = false;
			}
		}
		if (allIdentical) return;
	}
	std::vector<std::string> substrings;
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		assert(loopSequences[i].size() >= currentMatchIndices[i]);
		assert(currentMatchIndices[i] > previousMatchIndices[i]);
		substrings.emplace_back(loopSequences[i].substr(previousMatchIndices[i], currentMatchIndices[i]-previousMatchIndices[i]+1));
		assert(substrings.back().size() >= 2);
		maxClusterDistance = std::max<size_t>(maxClusterDistance, substrings[i].size() * maxDivergence);
	}
	std::vector<std::vector<size_t>> distanceMatrix;
	distanceMatrix.emplace_back();
	size_t maxDistance = 0;
	for (size_t i = 1; i < loopSequences.size(); i++)
	{
		distanceMatrix.emplace_back();
		for (size_t j = 0; j < i; j++)
		{
			distanceMatrix.back().emplace_back(getEditDistanceWfa(substrings[i], substrings[j]));
			maxDistance = std::max(maxDistance, distanceMatrix.back().back());
		}
	}
	if (maxDistance < maxClusterDistance) return;
	std::vector<size_t> parent;
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		parent.emplace_back(i);
	}
	for (size_t i = 1; i < loopSequences.size(); i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (distanceMatrix[i][j] > maxClusterDistance) continue;
			merge(parent, i, j);
		}
	}
	phmap::flat_hash_set<size_t> distinctClusters;
	std::vector<std::vector<size_t>> clusters;
	clusters.resize(loopSequences.size());
	for (size_t i = 0; i < parent.size(); i++)
	{
		clusters[find(parent, i)].emplace_back(i);
		distinctClusters.insert(find(parent, i));
	}
	if (distinctClusters.size() < 2) return;
	bool allClustersAreCliques = true;
	bool hasSmallCluster = false;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		if (clusters[i].size() > 0 && clusters[i].size() < minCoverage)
		{
			hasSmallCluster = true;
			break;
		}
		for (size_t j = 1; j < clusters[i].size(); j++)
		{
			for (size_t k = 0; k < j; k++)
			{
				size_t smaller = std::min(clusters[i][j], clusters[i][k]);
				size_t bigger = std::max(clusters[i][j], clusters[i][k]);
				if (distanceMatrix[bigger][smaller] > maxClusterDistance)
				{
					allClustersAreCliques = false;
					break;
				}
			}
			if (!allClustersAreCliques) break;
		}
		if (!allClustersAreCliques) break;
	}
	if (hasSmallCluster) return;
	if (!allClustersAreCliques) return;
	// distance between clusters is always at least maxClusterDistance+1 because of union-find
	for (size_t i = 0; i < clusters.size(); i++)
	{
		for (size_t read : clusters[i])
		{
			result[read].emplace_back(i);
		}
	}
}

std::vector<std::tuple<size_t, size_t, size_t>> getMatchRegionsRec(const std::string_view& refSequence, const std::string_view& querySequence, const size_t k)
{
	if (refSequence.size() < k) return std::vector<std::tuple<size_t, size_t, size_t>> {};
	if (querySequence.size() < k) return std::vector<std::tuple<size_t, size_t, size_t>> {};
	auto refKmers = getRefKmers(refSequence, k);
	auto kmerMatches = getKmerAnchors(refSequence, refKmers, querySequence, k);
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
			auto recResult = getMatchRegionsRec(std::string_view { refSequence.begin()+lastRefPos+1, kmerMatches[i].first+k-2 - lastRefPos }, std::string_view { querySequence.begin()+lastQueryPos+1, kmerMatches[i].second+k-2 - lastQueryPos }, k);
			if (recResult.size() == 0 && k > 11)
			{
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
				if (std::get<0>(t)+lastRefPos >= kmerMatches[i].first) continue;
				if (std::get<1>(t)+lastQueryPos >= kmerMatches[i].second) continue;
				if (i > 0 && std::get<0>(t)+std::get<2>(t)+lastRefPos <= kmerMatches[i-1].first+k) continue;
				if (i > 0 && std::get<1>(t)+std::get<2>(t)+lastQueryPos <= kmerMatches[i-1].second+k) continue;
				result.emplace_back(std::get<0>(t)+lastRefPos, std::get<1>(t)+lastQueryPos, std::get<2>(t));
			}
		}
		lastRefPos = kmerMatches[i].first;
		lastQueryPos = kmerMatches[i].second;
		result.emplace_back(kmerMatches[i].first, kmerMatches[i].second, k);
	}
	std::vector<std::tuple<size_t, size_t, size_t>> recResult;
	if (kmerMatches.size() >= 1) recResult = getMatchRegionsRec(std::string_view { refSequence.begin()+lastRefPos+1, refSequence.size() - lastRefPos - 1 }, std::string_view { querySequence.begin()+lastQueryPos + 1, querySequence.size() - lastQueryPos-1 }, k);
	if (recResult.size() == 0 && k > 11)
	{
		recResult = getMatchRegionsRec(std::string_view { refSequence.begin()+lastRefPos, refSequence.size() - lastRefPos }, std::string_view { querySequence.begin()+lastQueryPos, querySequence.size() - lastQueryPos }, k-10);
	}
	for (auto t : recResult)
	{
		if (kmerMatches.size() > 0 && std::get<0>(t)+std::get<2>(t)+lastRefPos <= kmerMatches.back().first+k) continue;
		if (kmerMatches.size() > 0 && std::get<1>(t)+std::get<2>(t)+lastQueryPos <= kmerMatches.back().second+k) continue;
		result.emplace_back(std::get<0>(t)+lastRefPos, std::get<1>(t)+lastQueryPos, std::get<2>(t));
	}
	return result;
}

std::vector<std::tuple<size_t, size_t, size_t>> getMatchRegions(const std::string& refSequence, const std::string& querySequence)
{
	auto kmerMatches = getMatchRegionsRec(refSequence, querySequence, 31);
	if (kmerMatches.size() == 0) return kmerMatches;
	std::sort(kmerMatches.begin(), kmerMatches.end());
	for (size_t j = 1; j < kmerMatches.size(); j++)
	{
		assert(std::get<0>(kmerMatches[j]) > std::get<0>(kmerMatches[j-1]));
		assert(std::get<1>(kmerMatches[j]) > std::get<1>(kmerMatches[j-1]));
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
	return diagonals;
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

std::vector<std::vector<size_t>> getMatchBases(const std::string& consensusSeq, const std::vector<std::string>& loopSequences)
{
	std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> matchRegions;
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		matchRegions.emplace_back(getMatchRegions(consensusSeq, loopSequences[i]));
		if (matchRegions.back().size() == 0) return std::vector<std::vector<size_t>>{};
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
	if (refAllowedAreas.size() == 0) return std::vector<std::vector<size_t>> {};
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "match positions in consensus (consensus len " << consensusSeq.size() << " size " << loopSequences.size() << "):";
	for (size_t i = 0; i < refAllowedAreas.size(); i++)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << refAllowedAreas[i].first << "-" << refAllowedAreas[i].second;
	}
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
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

std::vector<std::vector<size_t>> splitByBubbleAlleles(const std::vector<OntLoop>& loops, const std::vector<std::string>& loopSequences, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t numThreads)
{
	auto path = getConsensusPath(loops, graph);
	std::string consensusSeq = getConsensusSequence(path, graph, pathStartClip, pathEndClip);
	consensusSeq = polishConsensus(consensusSeq, loopSequences, numThreads);
	std::vector<std::vector<size_t>> kmerChain = getMatchBases(consensusSeq, loopSequences);
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "cluster with " << loopSequences.size() << " reads has " << kmerChain.size() << " all-present kmers in chain" << std::endl;
	std::vector<std::vector<size_t>> allelesPerRead;
	allelesPerRead.resize(loopSequences.size());
	for (size_t i = 1; i < kmerChain.size(); i++)
	{
		clusterBubbleAlleles(allelesPerRead, loopSequences, kmerChain[i-1], kmerChain[i]);
	}
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "cluster with " << loopSequences.size() << " reads has " << allelesPerRead[0].size() << " alleles" << std::endl;
	std::map<std::vector<size_t>, size_t> allelesToCluster;
	std::vector<std::vector<size_t>> result;
	for (size_t i = 0; i < allelesPerRead.size(); i++)
	{
		if (allelesToCluster.count(allelesPerRead[i]) == 0)
		{
			allelesToCluster[allelesPerRead[i]] = result.size();
			result.emplace_back();
		}
		result[allelesToCluster.at(allelesPerRead[i])].emplace_back(i);
	}
	return result;
}

std::string getConsensus(const std::vector<OntLoop>& loops, const std::vector<std::string>& loopSequences, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t numThreads)
{
	auto path = getConsensusPath(loops, graph);
	std::string consensusSeq = getConsensusSequence(path, graph, pathStartClip, pathEndClip);
	consensusSeq = polishConsensus(consensusSeq, loopSequences, numThreads);
	return consensusSeq;
}

std::vector<std::vector<std::tuple<size_t, size_t, std::string>>> getEditsForPhasing(const std::string& consensusSeq, const std::vector<std::string>& loopSequences, const size_t numThreads)
{
	size_t firstMatchPos = 0;
	size_t lastMatchPos = consensusSeq.size();;
	std::vector<std::vector<std::tuple<size_t, size_t, std::string>>> editsPerRead;
	editsPerRead.resize(loopSequences.size());
	std::mutex resultMutex;
	iterateEdits(consensusSeq, loopSequences, numThreads, [&firstMatchPos, &lastMatchPos, &resultMutex, &editsPerRead](const size_t threadId, const size_t readId, const size_t firstMatchRef, const size_t lastMatchRef, const size_t firstMatchRead, const size_t lastMatchRead, const std::tuple<size_t, size_t, std::string, size_t>& edit)
	{
		auto filterEdit = std::make_tuple(std::get<0>(edit), std::get<1>(edit), std::get<2>(edit));
		editsPerRead[readId].emplace_back(filterEdit);
		std::lock_guard<std::mutex> lock { resultMutex };
		firstMatchPos = std::max(firstMatchPos, firstMatchRef);
		lastMatchPos = std::min(lastMatchPos, lastMatchRef);
	});
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (size_t j = editsPerRead[i].size()-1; j < editsPerRead[i].size(); j++)
		{
			if (std::get<0>(editsPerRead[i][j]) < firstMatchPos || std::get<1>(editsPerRead[i][j]) > lastMatchPos)
			{
				std::swap(editsPerRead[i][j], editsPerRead[i].back());
				editsPerRead[i].pop_back();
			}
		}
	}
	return editsPerRead;
}

std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>> getSNPMSA(const std::vector<std::vector<std::tuple<size_t, size_t, std::string>>>& editsPerRead)
{
	const size_t minCoverage = std::max<size_t>(5, editsPerRead.size()*0.1);
	const size_t minCoverageWithDels = std::max<size_t>(10, editsPerRead.size()*0.1);
	std::vector<std::vector<std::vector<size_t>>> readsWhichHaveAlt;
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (size_t j = 0; j < editsPerRead[i].size(); j++)
		{
			while (readsWhichHaveAlt.size() < std::get<1>(editsPerRead[i][j])) readsWhichHaveAlt.emplace_back();
			const std::tuple<size_t, size_t, std::string>& t = editsPerRead[i][j];
			if (std::get<1>(t) == std::get<0>(t)) continue;
			if (std::get<2>(t).size() >= 1 && std::get<1>(t)-std::get<0>(t) != std::get<2>(t).size()) continue;
			for (size_t k = std::get<0>(t); k < std::get<1>(t); k++)
			{
				size_t allele = 5;
				if (k-std::get<0>(t) < std::get<2>(t).size())
				{
					switch(std::get<2>(t)[k-std::get<0>(t)])
					{
						case 'A':
							allele = 0;
							break;
						case 'C':
							allele = 1;
							break;
						case 'G':
							allele = 2;
							break;
						case 'T':
							allele = 3;
							break;
						default:
							assert(false);
							break;
					}
				}
				while (readsWhichHaveAlt[k].size() <= allele) readsWhichHaveAlt[k].emplace_back(0);
				readsWhichHaveAlt[k][allele].emplace_back(i);
			}
		}
	}
	std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>> result;
	for (size_t i = 0; i < readsWhichHaveAlt.size(); i++)
	{
		size_t totalCoverage = 0;
		bool hasSmallAllele = false;
		size_t totalNonDeletionAlleles = 0;
		bool hasDeletion = readsWhichHaveAlt[i].size() > 5 && readsWhichHaveAlt[i][5].size() >= 1;
		bool hasNeighboringDeletion = (i > 0 && readsWhichHaveAlt[i-1].size() > 5 && readsWhichHaveAlt[i-1][5].size() >= 1) || (i+1 < readsWhichHaveAlt.size() && readsWhichHaveAlt[i+1].size() > 5 && readsWhichHaveAlt[i+1][5].size() >= 1);
		for (size_t allele = 0; allele < readsWhichHaveAlt[i].size(); allele++)
		{
			totalCoverage += readsWhichHaveAlt[i][allele].size();
			if (readsWhichHaveAlt[i][allele].size() == 0) continue;
			if (allele != 5) totalNonDeletionAlleles += 1;
			if (hasDeletion && hasNeighboringDeletion && readsWhichHaveAlt[i][allele].size() < minCoverageWithDels && allele != 5)
			{
				hasSmallAllele = true;
			}
			if (readsWhichHaveAlt[i][allele].size() < minCoverage && allele != 5)
			{
				hasSmallAllele = true;
			}
		}
		if (hasSmallAllele) continue;
		if (totalCoverage+minCoverage > editsPerRead.size()) continue;
		if (hasDeletion && hasNeighboringDeletion && totalCoverage+minCoverageWithDels > editsPerRead.size()) continue;
		if (totalNonDeletionAlleles < 1) continue;
		result.emplace_back();
		result.back().second = i;
		result.back().first.resize(7);
		std::vector<bool> readHasAlt;
		readHasAlt.resize(editsPerRead.size(), false);
		for (size_t allele = 0; allele < readsWhichHaveAlt[i].size(); allele++)
		{
			for (size_t read : readsWhichHaveAlt[i][allele])
			{
				result.back().first[allele].emplace_back(read);
				readHasAlt[read] = true;
			}
		}
		for (size_t j = 0; j < readHasAlt.size(); j++)
		{
			if (readHasAlt[j]) continue;
			result.back().first[6].emplace_back(j);
		}
		assert(result.back().first[6].size() >= minCoverage);
		assert(!hasDeletion || !hasNeighboringDeletion || result.back().first[6].size() >= minCoverageWithDels);
		for (size_t j = result.back().first.size()-1; j < result.back().first.size(); j--)
		{
			if (result.back().first[j].size() != 0) continue;
			std::swap(result.back().first[j], result.back().first.back());
			result.back().first.pop_back();
		}
		assert(result.back().first.size() >= 2);
		assert((!hasDeletion && result.back().first.size() == totalNonDeletionAlleles+1) || (hasDeletion && result.back().first.size() == totalNonDeletionAlleles+2));
	}
	return result;
}

std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>> getPhasableVariantInfoBiallelicAltRefSNPsBigIndels(const std::vector<std::vector<std::tuple<size_t, size_t, std::string>>>& editsPerRead, const std::string& refSequence)
{
	const size_t minCoverage = 3;
	std::map<std::tuple<size_t, size_t, std::string>, size_t> editCounts;
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (size_t j = 0; j < editsPerRead[i].size(); j++)
		{
			editCounts[editsPerRead[i][j]] += 1;
		}
	}
	std::vector<std::tuple<size_t, size_t, std::string>> distinctEdits;
	for (auto pair : editCounts)
	{
		distinctEdits.emplace_back(pair.first);
	}
	std::sort(distinctEdits.begin(), distinctEdits.end());
	std::set<std::tuple<size_t, size_t, std::string>> discardedEdits;
	size_t currentClusterStart = 0;
	size_t currentClusterLastPos = std::get<1>(distinctEdits[0]);
	std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>> result;
	std::map<std::tuple<size_t, size_t, std::string>, std::pair<size_t, size_t>> validEditIndices;
	for (size_t i = 1; i <= distinctEdits.size(); i++)
	{
		if (i == distinctEdits.size() || std::get<0>(distinctEdits[i]) > currentClusterLastPos)
		{
			if (i-currentClusterStart <= 1)
			{
				size_t countSNPs = 0;
				size_t countBigIndels = 0;
				for (size_t j = currentClusterStart; j < i; j++)
				{
					if (isSNP(distinctEdits[j])) countSNPs += 1;
					if (isBigIndel(distinctEdits[j]) && !isMicrosatelliteOrHomopolymerIndel(distinctEdits[j], refSequence)) countBigIndels += 1;
				}
				if (countSNPs >= 1 || countBigIndels >= 1)
				{
					result.emplace_back();
					result.back().first.emplace_back(); // ref allele
					size_t minPos = std::get<0>(distinctEdits[currentClusterStart]);
					size_t maxPos = std::get<1>(distinctEdits[currentClusterStart]);
					for (size_t j = currentClusterStart; j < i; j++)
					{
						result.back().first.emplace_back();
						validEditIndices[distinctEdits[j]] = std::make_pair(result.size()-1, result.back().first.size()-1);
						minPos = std::min(minPos, std::get<0>(distinctEdits[j]));
						maxPos = std::max(maxPos, std::get<1>(distinctEdits[j]));
					}
					assert(result.back().first.size() >= 2);
					result.back().second = (minPos+maxPos)/2;
				}
			}
			currentClusterStart = i;
		}
		if (i != distinctEdits.size()) currentClusterLastPos = std::max(currentClusterLastPos, std::get<1>(distinctEdits[i]));
	}
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (const auto& t : editsPerRead[i])
		{
			if (validEditIndices.count(t) == 0) continue;
			auto index = validEditIndices.at(t);
			result[index.first].first[index.second].emplace_back(i);
		}
	}
	for (size_t i = result.size()-1; i < result.size(); i--)
	{
		bool hasMultipleAllelesFromSameRead = false;
		std::vector<bool> hasAlt;
		hasAlt.resize(editsPerRead.size(), false);
		for (size_t allele = 1; allele < result[i].first.size(); allele++)
		{
			for (size_t j : result[i].first[allele])
			{
				if (hasAlt[j])
				{
					hasMultipleAllelesFromSameRead = true;
					break;
				}
				hasAlt[j] = true;
			}
			if (hasMultipleAllelesFromSameRead) break;
		}
		if (hasMultipleAllelesFromSameRead)
		{
			std::swap(result[i], result.back());
			result.pop_back();
			continue;
		}
		for (size_t j = 0; j < hasAlt.size(); j++)
		{
			if (hasAlt[j]) continue;
			result[i].first[0].emplace_back(j);
		}
	}
	for (size_t i = result.size()-1; i < result.size(); i--)
	{
		bool hasLowCoverage = false;
		for (size_t j = 0; j < result[i].first.size(); j++)
		{
			if (result[i].first[j].size() >= minCoverage) continue;
			hasLowCoverage = true;
			break;
		}
		if (hasLowCoverage)
		{
			std::swap(result[i], result.back());
			result.pop_back();
		}
	}
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) { return left.second < right.second; });
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].first.size(); j++)
		{
			std::sort(result[i].first[j].begin(), result[i].first[j].end());
		}
	}
	return result;
}

bool vectorsMatch(const std::vector<size_t>& left, const std::vector<size_t>& right)
{
	if (left.size() != right.size()) return false;
	for (size_t i = 0; i < left.size(); i++)
	{
		if (left[i] != right[i]) return false;
	}
	return true;
}

std::vector<std::vector<size_t>> trySplitTwoSites(const std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>>& phasableVariantInfo, const size_t numReads)
{
	const size_t minDistance = 100;
	const size_t minCoverage = std::max<size_t>(5, numReads*0.05);
	size_t bestSiteOne = std::numeric_limits<size_t>::max();
	size_t bestSiteTwo = std::numeric_limits<size_t>::max();
	size_t bestSiteMinorAlleleCoverage = 0;
	std::vector<size_t> readsInFirstCluster;
	for (size_t i = 0; i < phasableVariantInfo.size(); i++)
	{
		for (size_t j = i+1; j < phasableVariantInfo.size(); j++)
		{
			assert(phasableVariantInfo[i].second < phasableVariantInfo[j].second);
			if (phasableVariantInfo[j].second < phasableVariantInfo[i].second + minDistance) continue;
			std::vector<std::pair<size_t, size_t>> allelesPerRead;
			allelesPerRead.resize(numReads, std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
			for (size_t k = 0; k < phasableVariantInfo[i].first.size(); k++)
			{
				for (size_t read : phasableVariantInfo[i].first[k])
				{
					assert(allelesPerRead[read].first == std::numeric_limits<size_t>::max());
					allelesPerRead[read].first = k;
				}
			}
			for (size_t k = 0; k < phasableVariantInfo[j].first.size(); k++)
			{
				for (size_t read : phasableVariantInfo[j].first[k])
				{
					assert(allelesPerRead[read].second == std::numeric_limits<size_t>::max());
					allelesPerRead[read].second = k;
				}
			}
			phmap::flat_hash_set<std::pair<size_t, size_t>> alleles;
			phmap::flat_hash_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> parent;
			for (auto t : allelesPerRead)
			{
				alleles.insert(t);
				parent[t] = t;
			}
			for (auto t : alleles)
			{
				for (auto t2 : alleles)
				{
					if (t2 == t) continue;
					size_t edits = 0;
					if (t.first != t2.first) edits += 1;
					if (t.second != t2.second) edits += 1;
					assert(edits > 0);
					if (edits == 1) merge(parent, t, t2);
				}
			}
			phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> clusters;
			for (auto t : alleles)
			{
				auto p = find(parent, t);
				if (clusters.count(p) == 1) continue;
				size_t cluster = clusters.size();
				clusters[p] = cluster;
			}
			if (clusters.size() != 2) continue;
			phmap::flat_hash_map<size_t, size_t> clusterCoverage;
			for (auto t : allelesPerRead)
			{
				clusterCoverage[clusters.at(find(parent, t))] += 1;
			}
			bool hasSmallCluster = false;
			size_t minorAlleleCoverage = numReads;
			for (auto pair : clusterCoverage)
			{
				minorAlleleCoverage = std::min(minorAlleleCoverage, pair.second);
			}
			assert(minorAlleleCoverage*2 <= numReads);
			if (minorAlleleCoverage < minCoverage) hasSmallCluster = true;
			if (hasSmallCluster) continue;
			if (minorAlleleCoverage > bestSiteMinorAlleleCoverage)
			{
				bestSiteMinorAlleleCoverage = minorAlleleCoverage;
				bestSiteOne = i;
				bestSiteTwo = j;
				readsInFirstCluster.clear();
				size_t firstCluster = std::numeric_limits<size_t>::max();
				for (size_t read = 0; read < allelesPerRead.size(); read++)
				{
					size_t cluster = clusters.at(find(parent, allelesPerRead[read]));
					if (firstCluster == std::numeric_limits<size_t>::max()) firstCluster = cluster;
					if (cluster == firstCluster) readsInFirstCluster.emplace_back(read);
				}
			}
		}
	}
	if (bestSiteMinorAlleleCoverage == 0) return std::vector<std::vector<size_t>> {};
	assert(bestSiteOne < phasableVariantInfo.size());
	assert(bestSiteTwo < phasableVariantInfo.size());
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "split two sites " << phasableVariantInfo[bestSiteOne].second << " " << phasableVariantInfo[bestSiteTwo].second << std::endl;
	std::vector<std::vector<size_t>> result;
	result.resize(2);
	result[0] = readsInFirstCluster;
	std::vector<bool> readUsed;
	readUsed.resize(numReads, false);
	for (size_t read : readsInFirstCluster) readUsed[read] = true;
	for (size_t read = 0; read < numReads; read++)
	{
		if (readUsed[read]) continue;
		result[1].emplace_back(read);
	}
	return result;
}

std::vector<std::vector<size_t>> trySplitThreeSites(const std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>>& phasableVariantInfo, const size_t numReads)
{
	const size_t minDistance = 100;
	const size_t minCoverage = std::max<size_t>(5, numReads*0.05);
	std::vector<std::vector<size_t>> allelesPerRead;
	for (size_t i = 0; i < phasableVariantInfo.size(); i++)
	{
		for (size_t j = 0; j < phasableVariantInfo[i].first.size(); j++)
		{
			for (auto read : phasableVariantInfo[i].first[j])
			{
				while (allelesPerRead.size() <= read)
				{
					allelesPerRead.emplace_back();
					allelesPerRead.back().resize(phasableVariantInfo.size(), std::numeric_limits<size_t>::max());
				}
				assert(allelesPerRead[read][i] == std::numeric_limits<size_t>::max());
				allelesPerRead[read][i] = j;
			}
		}
	}
	for (size_t i = 0; i < allelesPerRead.size(); i++)
	{
		for (size_t j = 0; j < allelesPerRead[i].size(); j++)
		{
			assert(allelesPerRead[i][j] != std::numeric_limits<size_t>::max());
		}
	}
	for (size_t blocki = 0; blocki*10 < (phasableVariantInfo.size()+9)/10; blocki++)
	{
		for (size_t blockj = blocki; blockj*10 < (phasableVariantInfo.size()+9)/10; blockj++)
		{
			for (size_t blockk = blockj; blockk*10 < (phasableVariantInfo.size()+9)/10; blockk++)
			{
				std::set<std::vector<size_t>> distinctBlocks;
				for (size_t i = 0; i < allelesPerRead.size(); i++)
				{
					std::vector<size_t> thisBlock;
					for (size_t j = 0; j < 10; j++)
					{
						if (blocki*10+j < phasableVariantInfo.size())
						{
							thisBlock.emplace_back(allelesPerRead[i][blocki*10+j]);
						}
						else
						{
							thisBlock.emplace_back(std::numeric_limits<size_t>::max());
						}
					}
					for (size_t j = 0; j < 10; j++)
					{
						if (blockj*10+j < phasableVariantInfo.size())
						{
							thisBlock.emplace_back(allelesPerRead[i][blockj*10+j]);
						}
						else
						{
							thisBlock.emplace_back(std::numeric_limits<size_t>::max());
						}
					}
					for (size_t j = 0; j < 10; j++)
					{
						if (blockk*10+j < phasableVariantInfo.size())
						{
							thisBlock.emplace_back(allelesPerRead[i][blockk*10+j]);
						}
						else
						{
							thisBlock.emplace_back(std::numeric_limits<size_t>::max());
						}
					}
					assert(thisBlock.size() == 30);
					distinctBlocks.emplace(thisBlock);
				}
				for (size_t offi = 0; offi < 10 && blocki*10+offi < phasableVariantInfo.size(); offi++)
				{
					for (size_t offj = (blocki == blockj ? offi+1 : 0); offj < 10 && blockj*10+offj < phasableVariantInfo.size(); offj++)
					{
						for (size_t offk = (blockj == blockk ? offj+1 : 0); offk < 10 && blockk*10+offk < phasableVariantInfo.size(); offk++)
						{
							size_t i = blocki*10+offi;
							size_t j = blockj*10+offj;
							size_t k = blockk*10+offk;
							assert(i < j);
							assert(j < k);
							assert(phasableVariantInfo[i].second < phasableVariantInfo[j].second);
							assert(phasableVariantInfo[j].second < phasableVariantInfo[k].second);
							if (phasableVariantInfo[j].second < phasableVariantInfo[i].second + minDistance) continue;
							if (phasableVariantInfo[k].second < phasableVariantInfo[j].second + minDistance) continue;
							std::set<std::tuple<size_t, size_t, size_t>> alleles;
							for (const auto& block : distinctBlocks)
							{
								alleles.emplace(block[offi], block[10+offj], block[20+offk]);
							}
							std::map<std::tuple<size_t, size_t, size_t>, std::tuple<size_t, size_t, size_t>> parent;
							for (auto t : alleles)
							{
								if (parent.count(t) == 0) parent[t] = t;
								for (auto t2 : alleles)
								{
									if (t == t2) continue;
									size_t matches = 0;
									if (std::get<0>(t) == std::get<0>(t2)) matches += 1;
									if (std::get<1>(t) == std::get<1>(t2)) matches += 1;
									if (std::get<2>(t) == std::get<2>(t2)) matches += 1;
									if (matches < 2) continue;
									while (parent[parent[t]] != parent[t]) parent[t] = parent[parent[t]];
									while (parent[parent[t2]] != parent[t2]) parent[t2] = parent[parent[t2]];
									parent[parent[t2]] = parent[t];
								}
							}
							std::map<std::tuple<size_t, size_t, size_t>, size_t> newClusters;
							for (auto t : alleles)
							{
								while (parent[parent[t]] != parent[t]) parent[t] = parent[parent[t]];
								if (newClusters.count(parent[t]) == 1) continue;
								size_t cluster = newClusters.size();
								newClusters[parent[t]] = cluster;
							}
							if (newClusters.size() < 2) continue;
							assert(i < phasableVariantInfo.size());
							assert(j < phasableVariantInfo.size());
							assert(k < phasableVariantInfo.size());
							std::vector<std::vector<size_t>> result;
							result.resize(newClusters.size());
							for (size_t m = 0; m < allelesPerRead.size(); m++)
							{
								result[newClusters.at(parent.at(std::make_tuple(allelesPerRead[m][i], allelesPerRead[m][j], allelesPerRead[m][k])))].emplace_back(m);
							}
							bool allClustersCovered = true;
							for (size_t m = 0; m < result.size(); m++)
							{
								if (result[m].size() < minCoverage) allClustersCovered = false;
							}
							if (allClustersCovered)
							{
								Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "split three sites " << phasableVariantInfo[i].second << " " << phasableVariantInfo[j].second << " " << phasableVariantInfo[k].second << std::endl;
								return result;
							}
						}
					}
				}
			}
		}
	}
	return std::vector<std::vector<size_t>> {};
}

std::vector<std::vector<size_t>> splitAndAdd(const std::vector<OntLoop>& previous, const std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>>& phasableVariantInfo)
{
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "recurse phase of cluster with size " << previous.size() << std::endl;
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "count sites " << phasableVariantInfo.size() << std::endl;
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "try split two sites" << std::endl;
	auto startTime = getTime();
	std::vector<std::vector<size_t>> splitted = trySplitTwoSites(phasableVariantInfo, previous.size());
	auto endTime = getTime();
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "two sites took " << formatTime(startTime, endTime) << std::endl;
	if (splitted.size() >= 2)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "two sites splitted to clusters of sizes";
		for (size_t i = 0; i < splitted.size(); i++)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << splitted[i].size();
		}
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
		return splitted;
	}
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "try split three sites" << std::endl;
	startTime = getTime();
	splitted = trySplitThreeSites(phasableVariantInfo, previous.size());
	endTime = getTime();
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "three sites took " << formatTime(startTime, endTime) << std::endl;
	if (splitted.size() >= 2)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "three sites splitted to clusters of sizes";
		for (size_t i = 0; i < splitted.size(); i++)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << splitted[i].size();
		}
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
		return splitted;
	}
	std::vector<std::vector<size_t>> result;
	result.emplace_back();
	for (size_t i = 0; i < previous.size(); i++)
	{
		result.back().emplace_back(i);
	}
	return result;
}

void filterUniqueRegionsToSpan(std::vector<std::pair<size_t, size_t>>& uniqueRegionsInRef, const size_t minPos, const size_t maxPos)
{
	for (size_t i = uniqueRegionsInRef.size()-1; i < uniqueRegionsInRef.size(); i--)
	{
		if (uniqueRegionsInRef[i].first < minPos) uniqueRegionsInRef[i].first = minPos;
		if (uniqueRegionsInRef[i].second > maxPos) uniqueRegionsInRef[i].second = maxPos;
		if (uniqueRegionsInRef[i].first >= uniqueRegionsInRef[i].second)
		{
			std::swap(uniqueRegionsInRef[i], uniqueRegionsInRef.back());
			uniqueRegionsInRef.pop_back();
		}
	}
	std::sort(uniqueRegionsInRef.begin(), uniqueRegionsInRef.end());
}

void filterMatchRegionsToUniqueRegions(const std::vector<std::pair<size_t, size_t>>& uniqueRegionsInRef, std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& matchRegionsInReads)
{
	for (size_t i = 0; i < matchRegionsInReads.size(); i++)
	{
		phmap::flat_hash_set<std::tuple<size_t, size_t, size_t>> newResult;
		for (size_t j = matchRegionsInReads[i].size()-1; j < matchRegionsInReads[i].size(); j--)
		{
			for (auto t : uniqueRegionsInRef)
			{
				if (t.first > std::get<1>(matchRegionsInReads[i][j])) continue;
				if (t.second < std::get<0>(matchRegionsInReads[i][j])) continue;
				size_t leftClip = 0;
				size_t rightClip = 0;
				if (t.first > std::get<0>(matchRegionsInReads[i][j])) leftClip = t.first - std::get<0>(matchRegionsInReads[i][j]);
				if (t.second < std::get<1>(matchRegionsInReads[i][j])) rightClip = std::get<1>(matchRegionsInReads[i][j]) - t.second;
				if (leftClip + rightClip < std::get<1>(matchRegionsInReads[i][j]) - std::get<0>(matchRegionsInReads[i][j]))
				{
					newResult.emplace(std::get<0>(matchRegionsInReads[i][j]) + leftClip, std::get<1>(matchRegionsInReads[i][j]) - rightClip, std::get<2>(matchRegionsInReads[i][j]) + leftClip);
				}
			}
		}
		matchRegionsInReads[i].clear();
		matchRegionsInReads[i].insert(matchRegionsInReads[i].end(), newResult.begin(), newResult.end());
		std::sort(matchRegionsInReads[i].begin(), matchRegionsInReads[i].end());
	}
}

std::vector<std::pair<size_t, size_t>> getAllmatchRegions(const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& matchRegionsInReads)
{
	std::vector<std::pair<size_t, size_t>> result;
	std::vector<size_t> indexPerRead;
	indexPerRead.resize(matchRegionsInReads.size(), 0);
	while (true)
	{
		size_t maxMatchStart = 0;
		size_t minMatchEnd = std::numeric_limits<size_t>::max();
		bool shouldBreak = false;
		for (size_t i = 0; i < matchRegionsInReads.size(); i++)
		{
			if (indexPerRead[i] == matchRegionsInReads[i].size())
			{
				shouldBreak = true;
				break;
			}
			assert(indexPerRead[i] < matchRegionsInReads[i].size());
			maxMatchStart = std::max(maxMatchStart, std::get<0>(matchRegionsInReads[i][indexPerRead[i]]));
			minMatchEnd = std::min(minMatchEnd, std::get<1>(matchRegionsInReads[i][indexPerRead[i]]));
		}
		if (shouldBreak) break;
		if (minMatchEnd < maxMatchStart)
		{
			for (size_t i = 0; i < matchRegionsInReads.size(); i++)
			{
				if (indexPerRead[i] == matchRegionsInReads[i].size()) continue;
				while (indexPerRead[i] < matchRegionsInReads[i].size() && std::get<1>(matchRegionsInReads[i][indexPerRead[i]]) < maxMatchStart) indexPerRead[i] += 1;
				if (indexPerRead[i] == matchRegionsInReads[i].size())
				{
					shouldBreak = true;
					break;
				}
			}
			continue;
		}
		if (shouldBreak) break;
		result.emplace_back(maxMatchStart, minMatchEnd);
		for (size_t i = 0; i < matchRegionsInReads.size(); i++)
		{
			if (indexPerRead[i] == matchRegionsInReads[i].size()) continue;
			if (std::get<1>(matchRegionsInReads[i][indexPerRead[i]]) == minMatchEnd) indexPerRead[i] += 1;
		}
	}
	return result;
}

std::vector<std::vector<size_t>> getReadSequenceLengthsBetweenAllMatchRegions(const std::vector<std::pair<size_t, size_t>>& allmatchRegions, const std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>& matchRegionsInReads)
{
	std::vector<std::vector<size_t>> result;
	result.resize(matchRegionsInReads.size());
	for (size_t read = 0; read < matchRegionsInReads.size(); read++)
	{
		size_t readIndex = 0;
		size_t lastMatchReadPos = std::numeric_limits<size_t>::max();
		for (size_t allmatchIndex = 0; allmatchIndex < allmatchRegions.size(); allmatchIndex++)
		{
			while (std::get<1>(matchRegionsInReads[read][readIndex]) < allmatchRegions[allmatchIndex].first)
			{
				readIndex += 1;
				assert(readIndex < matchRegionsInReads[read].size());
			}
			assert(std::get<0>(matchRegionsInReads[read][readIndex]) <= allmatchRegions[allmatchIndex].first);
			assert(std::get<1>(matchRegionsInReads[read][readIndex]) >= allmatchRegions[allmatchIndex].second);
			if (allmatchIndex > 0)
			{
				assert(lastMatchReadPos != std::numeric_limits<size_t>::max());
				result[read].emplace_back(std::get<2>(matchRegionsInReads[read][readIndex]) + allmatchRegions[allmatchIndex].first - std::get<0>(matchRegionsInReads[read][readIndex]) - lastMatchReadPos);
			}
			lastMatchReadPos = std::get<2>(matchRegionsInReads[read][readIndex]) + allmatchRegions[allmatchIndex].second - std::get<0>(matchRegionsInReads[read][readIndex]);
		}
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		assert(result[i].size() == allmatchRegions.size()-1);
	}
	return result;
}

std::vector<std::vector<size_t>> getReadLengthAlleles(const std::vector<std::pair<size_t, size_t>>& allmatchRegions, const std::vector<std::vector<size_t>>& readVariantLengths)
{
	std::vector<std::vector<size_t>> result;
	result.resize(readVariantLengths.size());
	for (size_t i = 0; i+1 < allmatchRegions.size(); i++)
	{
		size_t minSeparation = std::max<size_t>(20, (std::get<0>(allmatchRegions[i+1])-std::get<1>(allmatchRegions[i])) * 0.05);
		std::vector<size_t> lengths;
		for (size_t j = 0; j < readVariantLengths.size(); j++)
		{
			lengths.emplace_back(readVariantLengths[j][i]);
		}
		std::sort(lengths.begin(), lengths.end());
		std::vector<size_t> separators;
		for (size_t i = 1; i < lengths.size(); i++)
		{
			if (lengths[i] > lengths[i-1]+minSeparation)
			{
				separators.emplace_back((lengths[i] + lengths[i-1])/2);
			}
		}
		if (separators.size() == 0) continue;
		separators.emplace_back(lengths.back()+1);
		for (size_t j = 0; j < readVariantLengths.size(); j++)
		{
			size_t index = std::upper_bound(separators.begin(), separators.end(), readVariantLengths[j][i]) - separators.begin();
			assert(index < separators.size());
			result[j].emplace_back(index);
		}
	}
	return result;
}

std::vector<std::vector<size_t>> splitByBigIndels(const std::vector<std::vector<std::tuple<size_t, size_t, std::string>>>& editsPerRead)
{
	const size_t minLengthClusterCoverage = 10;
	if (editsPerRead.size() <= 10)
	{
		std::vector<std::vector<size_t>> result;
		result.emplace_back();
		for (size_t i = 0; i < editsPerRead.size(); i++)
		{
			result.back().emplace_back(i);
		}
		return result;
	}
	std::set<std::pair<size_t, size_t>> closeToIndelSpans;
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (auto t : editsPerRead[i])
		{
			if (!isIndel(t)) continue;
			closeToIndelSpans.emplace(std::get<0>(t), std::get<1>(t));
		}
	}
	std::vector<std::pair<size_t, size_t>> closeToIndelSpansVec { closeToIndelSpans.begin(), closeToIndelSpans.end() };
	std::sort(closeToIndelSpansVec.begin(), closeToIndelSpansVec.end());
	std::vector<std::pair<size_t, size_t>> nonIndelSpans;
	size_t currentIndelSpanEnd = 0;
	for (auto pair : closeToIndelSpansVec)
	{
		if (pair.first > currentIndelSpanEnd + 10)
		{
			nonIndelSpans.emplace_back(currentIndelSpanEnd, pair.first);
		}
		currentIndelSpanEnd = std::max(currentIndelSpanEnd, pair.second);
	}
	nonIndelSpans.emplace_back(currentIndelSpanEnd, currentIndelSpanEnd+2);
	assert(nonIndelSpans.size() >= 1);
	std::vector<size_t> indelSpanSeparators;
	for (size_t i = 1; i < nonIndelSpans.size(); i++)
	{
		indelSpanSeparators.emplace_back((nonIndelSpans[i].first + nonIndelSpans[i].second)/2);
	}
	std::vector<std::vector<int>> indelLengthSums;
	indelLengthSums.resize(editsPerRead.size());
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		indelLengthSums[i].resize(indelSpanSeparators.size(), 0);
		for (auto t : editsPerRead[i])
		{
			if (!isIndel(t)) continue;
			size_t index = std::upper_bound(indelSpanSeparators.begin(), indelSpanSeparators.end(), std::get<0>(t)) - indelSpanSeparators.begin();
			assert(index < indelLengthSums[i].size());
			indelLengthSums[i][index] += std::get<2>(t).size();
			indelLengthSums[i][index] -= std::get<1>(t);
			indelLengthSums[i][index] += std::get<0>(t);
		}
	}
	std::vector<std::vector<size_t>> indelAllelesPerRead;
	indelAllelesPerRead.resize(editsPerRead.size());
	for (size_t j = 0; j < indelSpanSeparators.size(); j++)
	{
		std::vector<int> indelsPerReads;
		for (size_t i = 0; i < editsPerRead.size(); i++)
		{
			indelsPerReads.emplace_back(indelLengthSums[i][j]);
		}
		std::sort(indelsPerReads.begin(), indelsPerReads.end());
		assert(j+1 < nonIndelSpans.size());
		assert(nonIndelSpans[j+1].first >= nonIndelSpans[j].second);
		std::vector<int> separators;
		int requiredLengthDifference = std::max<int>(15, (nonIndelSpans[j+1].first - nonIndelSpans[j].second) * 0.05);
		if (indelsPerReads.back() - indelsPerReads[0] < requiredLengthDifference) continue;
		bool hasSmallSeparator = false;
		for (size_t i = 1; i < indelsPerReads.size(); i++)
		{
			if (indelsPerReads[i] < indelsPerReads[i-1] + requiredLengthDifference) continue;
			size_t lastSeparatorIndex = 0;
			if (separators.size() > 0) lastSeparatorIndex = separators.back();
			separators.emplace_back((indelsPerReads[i] + indelsPerReads[i-1]) / 2);
			if (i - lastSeparatorIndex < minLengthClusterCoverage || i+minLengthClusterCoverage > indelsPerReads.size())
			{
				hasSmallSeparator = true;
				break;
			}
			for (size_t j = lastSeparatorIndex+1; j < i; j++)
			{
				if (indelsPerReads[j] > indelsPerReads[j-1]+10)
				{
					hasSmallSeparator = true;
					break;
				}
			}
			if (hasSmallSeparator) break;
		}
		if (hasSmallSeparator) continue;
		if (separators.size() == 0) continue;
		for (size_t i = 0; i < editsPerRead.size(); i++)
		{
			size_t allele = std::upper_bound(separators.begin(), separators.end(), indelLengthSums[i][j]) - separators.begin();
			assert(allele <= separators.size());
			indelAllelesPerRead[i].emplace_back(allele);
		}
	}
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "indel sites " << indelAllelesPerRead[0].size() << std::endl;
	std::map<std::vector<size_t>, size_t> allelesToCluster;
	std::vector<std::vector<size_t>> result;
	for (size_t i = 0; i < indelAllelesPerRead.size(); i++)
	{
		if (allelesToCluster.count(indelAllelesPerRead[i]) == 0)
		{
			allelesToCluster[indelAllelesPerRead[i]] = result.size();
			result.emplace_back();
		}
		result[allelesToCluster.at(indelAllelesPerRead[i])].emplace_back(i);
	}
	return result;
}

bool minorAllelesPerfectlyLinked(const phmap::flat_hash_map<std::tuple<size_t, size_t, std::string>, std::vector<size_t>>& readsWithEdit, const std::tuple<size_t, size_t, std::string>& editLeft, const std::tuple<size_t, size_t, std::string>& editRight, const size_t countReads)
{
	const size_t minCoverage = 10;
	std::vector<std::pair<bool, bool>> allelesPerRead;
	allelesPerRead.resize(countReads);
	for (size_t read : readsWithEdit.at(editLeft))
	{
		allelesPerRead[read].first = true;
	}
	for (size_t read : readsWithEdit.at(editRight))
	{
		allelesPerRead[read].second = true;
	}
	std::map<std::pair<bool, bool>, size_t> alleleCounts;
	for (size_t i = 0; i < allelesPerRead.size(); i++)
	{
		alleleCounts[allelesPerRead[i]] += 1;
	}
	bool hasSmallAllele = false;
	for (auto pair : alleleCounts)
	{
		if (pair.second < minCoverage) hasSmallAllele = true;
	}
	if (hasSmallAllele) return false;
	assert(alleleCounts.size() >= 2);
	assert(alleleCounts.size() <= 4);
	if (alleleCounts.size() == 4) return false;
	if (alleleCounts.size() == 2)
	{
		for (auto t1 : alleleCounts)
		{
			for (auto t2 : alleleCounts)
			{
				if (t1.first == t2.first) continue;
				size_t edits = 0;
				if (t1.first.first != t2.first.first) edits += 1;
				if (t1.first.second != t2.first.second) edits += 1;
				assert(edits == 1 || edits == 2);
				if (edits == 2) Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "minor alleles linked 2 edits" << std::endl;
				return edits == 2;
			}
		}
	}
	assert(alleleCounts.size() == 3);
	bool leftIsMinor = readsWithEdit.at(editLeft).size()*2 < countReads;
	bool rightIsMinor = readsWithEdit.at(editRight).size()*2 < countReads;
	if (leftIsMinor ^ rightIsMinor)
	{
		if (alleleCounts.count(std::make_pair(false, false)) == 0) return true;
		if (alleleCounts.count(std::make_pair(true, true)) == 0) return true;
		return false;
	}
	if (alleleCounts.count(std::make_pair(true, false)) == 0) return true;
	if (alleleCounts.count(std::make_pair(false, true)) == 0) return true;
	return false;
}

bool editIntersectsWithSomething(const std::tuple<size_t, size_t, std::string>& edit, const std::vector<std::tuple<size_t, size_t, std::string>>& otherEdits)
{
	for (size_t i = 0; i < otherEdits.size(); i++)
	{
		if (std::get<0>(edit) <= std::get<1>(otherEdits[i]) && std::get<1>(edit) >= std::get<0>(otherEdits[i])) return true;
	}
	return false;
}

size_t intersectSize(const std::vector<size_t>& left, const std::vector<size_t>& right)
{
	size_t result = 0;
	size_t leftIndex = 0;
	size_t rightIndex = 0;
	while (leftIndex < left.size() && rightIndex < right.size())
	{
		if (left[leftIndex] == right[rightIndex])
		{
			result += 1;
			leftIndex += 1;
			rightIndex += 1;
			continue;
		}
		if (left[leftIndex] > right[rightIndex])
		{
			rightIndex += 1;
			continue;
		}
		assert(left[leftIndex] < right[rightIndex]);
		leftIndex += 1;
	}
	return result;
}

std::vector<double> pascalTriangleFractionLog(const size_t row)
{
	// https://stackoverflow.com/questions/15580291/how-to-efficiently-calculate-a-row-in-pascals-triangle
	std::vector<double> result;
	result.emplace_back(0);
	for (size_t i = 0; i < row; i++)
	{
		result.emplace_back(result.back() + log(row-i) - log(i+1));
	}
	return result;
}

double binomialPValue(const double p, const size_t success, const size_t trials)
{
	std::vector<double> pascalLogRow = pascalTriangleFractionLog(trials);
	double normalizer = 0;
	for (size_t i = 0; i <= trials; i++)
	{
		normalizer += exp(pascalLogRow[i] + log(p) * success + log(1-p) * (trials-success));
	}
	if (success > trials * p)
	{
		double result = 0;
		for (size_t i = success; i <= trials; i++)
		{
			result += exp(pascalLogRow[i] + log(p) * success + log(1-p) * (trials-success));
		}
		return result / normalizer;
	}
	double result = 0;
	for (size_t i = 0; i <= success; i++)
	{
		result += exp(pascalLogRow[i] + log(p) * success + log(1-p) * (trials-success));
	}
	return result / normalizer;
}

bool allelesMatchWellEnough(const std::vector<std::vector<size_t>>& readsWithEdit, const std::vector<std::vector<size_t>>& readsWithRef, const size_t left, const size_t right)
{
	size_t countRefRef = intersectSize(readsWithRef[left], readsWithRef[right]);
	size_t countRefAlt = intersectSize(readsWithRef[left], readsWithEdit[right]);
	size_t countAltRef = intersectSize(readsWithEdit[left], readsWithRef[right]);
	size_t countAltAlt = intersectSize(readsWithEdit[left], readsWithEdit[right]);
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "check sites " << countRefRef << " " << countRefAlt << " " << countAltRef << " " << countAltAlt << std::endl;
	double refFractionLeft = (double)(countRefRef + countRefAlt) / (double)(countRefRef + countRefAlt + countAltRef + countAltAlt);
	double refFractionRight = (double)(countRefRef + countAltRef) / (double)(countRefRef + countRefAlt + countAltRef + countAltAlt);
	double requiredPValue = 0.000001;
	double leftPpart1 = binomialPValue(refFractionLeft, countRefRef, countRefRef + countAltRef);
	double leftPpart2 = binomialPValue(refFractionLeft, countRefAlt, countRefAlt + countAltAlt);
	double rightPpart1 = binomialPValue(refFractionRight, countRefRef, countRefRef + countRefAlt);
	double rightPpart2 = binomialPValue(refFractionRight, countAltRef, countAltRef + countAltAlt);
	bool result = false;
	if (leftPpart1 < requiredPValue) result = true;
	if (leftPpart2 < requiredPValue) result = true;
	if (rightPpart1 < requiredPValue) result = true;
	if (rightPpart2 < requiredPValue) result = true;
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "P-values " << leftPpart1 << " " << leftPpart2 << " " << rightPpart1 << " " << rightPpart2 << " result " << (result ? "yes" : "no") << std::endl;
	return result;
}

std::vector<std::vector<size_t>> splitBySNPCorrelation(const std::vector<std::vector<std::tuple<size_t, size_t, std::string>>>& editsPerRead, const std::string& refSequence)
{
	const size_t minCoverage = std::max<size_t>(10, editsPerRead.size()*0.2);
	std::map<std::tuple<size_t, size_t, std::string>, size_t> editCoverage;
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (const auto& edit : editsPerRead[i])
		{
			editCoverage[edit] += 1;
		}
	}
	std::set<std::tuple<size_t, size_t, std::string>> coveredEdits;
	for (auto pair : editCoverage)
	{
		if (pair.second < minCoverage) continue;
		if (pair.second+minCoverage > editsPerRead.size()) continue;
		coveredEdits.insert(pair.first);
	}
	std::map<std::tuple<size_t, size_t, std::string>, size_t> refCoverage;
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (const auto& edit : coveredEdits)
		{
			if (editIntersectsWithSomething(edit, editsPerRead[i])) continue;
			refCoverage[edit] += 1;
		}
	}
	std::vector<std::tuple<size_t, size_t, std::string>> goodishEdits;
	std::vector<size_t> veryGoodEditIndices;
	std::vector<size_t> perfectEditIndices;
	for (const auto& t : coveredEdits)
	{
		if (refCoverage.count(t) == 0) continue;
		if (refCoverage.at(t) < minCoverage) continue;
		if (editCoverage.at(t) + refCoverage.at(t) < editsPerRead.size() * 0.95) continue;
		assert(editCoverage.at(t) >= minCoverage);
		assert(editCoverage.at(t) + refCoverage.at(t) <= editsPerRead.size());
		if (isSNP(t) || isBigIndel(t))
		{
			if (!isMicrosatelliteOrHomopolymerIndel(t, refSequence))
			{
				if (editCoverage.at(t) + refCoverage.at(t) == editsPerRead.size())
				{
					Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "add perfect index " << goodishEdits.size() << " " << editCoverage.at(t) << " " << refCoverage.at(t) << std::endl;
					perfectEditIndices.emplace_back(goodishEdits.size());
				}
				else if (editCoverage.at(t) + refCoverage.at(t) + 3 > editsPerRead.size())
				{
					veryGoodEditIndices.emplace_back(goodishEdits.size());
				}
			}
		}
		goodishEdits.emplace_back(t);
	}
	std::sort(veryGoodEditIndices.begin(), veryGoodEditIndices.end(), [&editCoverage, &refCoverage, &goodishEdits](size_t left, size_t right)
	{
		return editCoverage.at(goodishEdits[left]) + refCoverage.at(goodishEdits[left]) < editCoverage.at(goodishEdits[right]) + refCoverage.at(goodishEdits[right]);
	});
	std::map<std::tuple<size_t, size_t, std::string>, size_t> goodishEditIndex;
	for (size_t i = 0; i < goodishEdits.size(); i++)
	{
		goodishEditIndex[goodishEdits[i]] = i;
	}
	std::vector<std::vector<size_t>> readsWithEdit;
	std::vector<std::vector<size_t>> readsWithRef;
	readsWithEdit.resize(goodishEdits.size());
	readsWithRef.resize(goodishEdits.size());
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		std::vector<bool> editFound;
		editFound.resize(goodishEdits.size(), false);
		for (const auto& edit : editsPerRead[i])
		{
			if (goodishEditIndex.count(edit) == 0) continue;
			readsWithEdit[goodishEditIndex.at(edit)].emplace_back(i);
			editFound[goodishEditIndex.at(edit)] = true;
		}
		for (size_t j = 0; j < editFound.size(); j++)
		{
			if (editFound[j]) continue;
			if (editIntersectsWithSomething(goodishEdits[j], editsPerRead[i])) continue;
			readsWithRef[j].emplace_back(i);
		}
	}
	for (size_t i = 0; i < goodishEdits.size(); i++)
	{
		assert(readsWithEdit[i].size() == editCoverage.at(goodishEdits[i]));
		assert(readsWithRef[i].size() == refCoverage.at(goodishEdits[i]));
		assert(intersectSize(readsWithEdit[i], readsWithRef[i]) == 0);
	}
	for (size_t i : perfectEditIndices)
	{
		assert(i < readsWithEdit.size());
		assert(i < readsWithRef.size());
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::get<0>(goodishEdits[i]) << "-" << std::get<1>(goodishEdits[i]) << "\"" << std::get<2>(goodishEdits[i]) << "\" " << i << " " << readsWithEdit[i].size() << " " << readsWithRef[i].size() << " " << editsPerRead.size() << std::endl;
		assert(readsWithEdit[i].size() + readsWithRef[i].size() == editsPerRead.size());
		for (size_t j = 0; j < goodishEdits.size(); j++)
		{
			if (i == j) continue;
			if (std::get<0>(goodishEdits[j]) < std::get<1>(goodishEdits[i])+100 && std::get<0>(goodishEdits[i]) < std::get<1>(goodishEdits[j])+100) continue;
			if (!allelesMatchWellEnough(readsWithEdit, readsWithRef, i, j)) continue;
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "split by perfect edits: " << std::get<0>(goodishEdits[i]) << "-" << std::get<1>(goodishEdits[i]) << "\"" << std::get<2>(goodishEdits[i]) << "\" " << std::get<0>(goodishEdits[j]) << "-" << std::get<1>(goodishEdits[j]) << "\"" << std::get<2>(goodishEdits[j]) << std::endl;
			std::vector<std::vector<size_t>> result;
			result.resize(2);
			for (size_t read : readsWithRef[i])
			{
				result[0].emplace_back(read);
			}
			for (size_t read : readsWithEdit[i])
			{
				result[1].emplace_back(read);
			}
			return result;
		}
	}
	for (size_t i : veryGoodEditIndices)
	{
		assert(readsWithEdit[i].size() + readsWithRef[i].size() < editsPerRead.size());
		assert(readsWithEdit[i].size() + readsWithRef[i].size() + 3 >= editsPerRead.size());
		for (size_t j = 0; j < goodishEdits.size(); j++)
		{
			if (i == j) continue;
			if (std::get<0>(goodishEdits[j]) < std::get<1>(goodishEdits[i])+100 && std::get<0>(goodishEdits[i]) < std::get<1>(goodishEdits[j])+100) continue;
			if (!allelesMatchWellEnough(readsWithEdit, readsWithRef, i, j)) continue;
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "split by edits: " << std::get<0>(goodishEdits[i]) << "-" << std::get<1>(goodishEdits[i]) << "\"" << std::get<2>(goodishEdits[i]) << "\" " << std::get<0>(goodishEdits[j]) << "-" << std::get<1>(goodishEdits[j]) << "\"" << std::get<2>(goodishEdits[j]) << std::endl;
			std::vector<std::vector<size_t>> result;
			result.resize(2);
			std::vector<bool> readFound;
			readFound.resize(editsPerRead.size(), false);
			for (size_t read : readsWithRef[i])
			{
				result[0].emplace_back(read);
				assert(!readFound[read]);
				readFound[read] = true;
			}
			for (size_t read : readsWithEdit[i])
			{
				result[1].emplace_back(read);
				assert(!readFound[read]);
				readFound[read] = true;
			}
			for (size_t read = 0; read < editsPerRead.size(); read++)
			{
				if (readFound[read]) continue;
				result.emplace_back();
				result.back().emplace_back(read);
			}
			return result;
		}
	}
	std::vector<std::vector<size_t>> result;
	result.resize(1);
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		result[0].emplace_back(i);
	}
	return result;
}

std::vector<std::vector<size_t>> splitByLinkedMinorAlleles(const std::vector<std::vector<std::tuple<size_t, size_t, std::string>>>& editsPerRead, const std::string& refSequence)
{
	const size_t minCoverage = std::max<size_t>(10, editsPerRead.size()*0.1);
	const size_t minDistance = 100;
	if (editsPerRead.size() < 2*minCoverage)
	{
		std::vector<std::vector<size_t>> result;
		result.emplace_back();
		for (size_t i = 0; i < editsPerRead.size(); i++)
		{
			result.back().emplace_back(i);
		}
		return result;
	}
	phmap::flat_hash_map<std::tuple<size_t, size_t, std::string>, std::vector<size_t>> readsWithEdit;
	for (size_t i = 0; i < editsPerRead.size(); i++)
	{
		for (const auto& t : editsPerRead[i])
		{
			if (!isSNP(t) && !isBigIndel(t)) continue;
			if (isMicrosatelliteOrHomopolymerIndel(t, refSequence)) continue;
			readsWithEdit[t].emplace_back(i);
		}
	}
	std::vector<std::tuple<size_t, size_t, std::string>> coveredEdits;
	for (const auto& pair : readsWithEdit)
	{
		if (pair.second.size() < minCoverage) continue;
		if (pair.second.size()+minCoverage > editsPerRead.size()) continue;
		coveredEdits.emplace_back(pair.first);
	}
	std::sort(coveredEdits.begin(), coveredEdits.end());
	std::vector<std::vector<size_t>> allelesPerRead;
	allelesPerRead.resize(editsPerRead.size());
	for (size_t i = 0; i < coveredEdits.size(); i++)
	{
		for (size_t j = i+1; j < coveredEdits.size(); j++)
		{
			if (std::get<0>(coveredEdits[j]) < std::get<1>(coveredEdits[i])+minDistance) continue;
			if (!minorAllelesPerfectlyLinked(readsWithEdit, coveredEdits[i], coveredEdits[j], editsPerRead.size())) continue;
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "minor alleles perfectly linked poses " << std::get<0>(coveredEdits[i]) << "-" << std::get<1>(coveredEdits[i]) << "\"" << std::get<2>(coveredEdits[i]) << "\"" << " " << std::get<0>(coveredEdits[j]) << "-" << std::get<1>(coveredEdits[j]) << "\"" << std::get<2>(coveredEdits[j]) << "\"" << std::endl;
			std::vector<bool> hasFirstEdit;
			hasFirstEdit.resize(editsPerRead.size(), false);
			for (size_t read : readsWithEdit.at(coveredEdits[i]))
			{
				hasFirstEdit[read] = true;
				allelesPerRead[read].emplace_back(1);
			}
			for (size_t read = 0; read < hasFirstEdit.size(); read++)
			{
				if (hasFirstEdit[read]) continue;
				allelesPerRead[read].emplace_back(0);
			}
			std::vector<bool> hasSecondEdit;
			hasSecondEdit.resize(editsPerRead.size(), false);
			for (size_t read : readsWithEdit.at(coveredEdits[j]))
			{
				hasSecondEdit[read] = true;
				allelesPerRead[read].emplace_back(1);
			}
			for (size_t read = 0; read < hasSecondEdit.size(); read++)
			{
				if (hasSecondEdit[read]) continue;
				allelesPerRead[read].emplace_back(0);
			}
			break;
		}
		if (allelesPerRead[0].size() > 0) break;
	}
	std::map<std::vector<size_t>, size_t> allelesToCluster;
	std::vector<std::vector<size_t>> result;
	for (size_t i = 0; i < allelesPerRead.size(); i++)
	{
		if (allelesToCluster.count(allelesPerRead[i]) == 0)
		{
			size_t cluster = result.size();
			allelesToCluster[allelesPerRead[i]] = cluster;
			result.emplace_back();
		}
		result[allelesToCluster.at(allelesPerRead[i])].emplace_back(i);
	}
	return result;
}

std::vector<std::vector<size_t>> splitByEditLinkage(const std::vector<std::vector<std::tuple<size_t, size_t, std::string>>>& editsPerRead, const std::string& refSequence)
{
	const size_t minEditsWithinCluster = (editsPerRead.size() < 20) ? 4 : 5;
	const size_t minReadsWithinCluster = (editsPerRead.size() < 20) ? 4 : 5;
	const size_t minDistanceBetweenEdits = 100;
	phmap::flat_hash_map<std::tuple<size_t, size_t, std::string>, size_t> editToIndex;
	std::vector<std::tuple<size_t, size_t, std::string>> edits;
	std::vector<std::vector<size_t>> readsWithEdit;
	std::vector<size_t> editPosition;
	for (size_t read = 0; read < editsPerRead.size(); read++)
	{
		for (const auto& edit : editsPerRead[read])
		{
			if (editToIndex.count(edit) == 0)
			{
				edits.emplace_back(edit);
				edits.emplace_back(edit);
				editToIndex[edit] = readsWithEdit.size();
				readsWithEdit.emplace_back();
				readsWithEdit.emplace_back();
				editPosition.emplace_back((std::get<0>(edit) + std::get<1>(edit)) / 2);
			}
			readsWithEdit[editToIndex.at(edit)].emplace_back(read);
		}
	}
	assert(readsWithEdit.size() % 2 == 0);
	for (size_t i = 0; i < readsWithEdit.size(); i += 2)
	{
		std::vector<bool> hasEdit;
		hasEdit.resize(editsPerRead.size(), false);
		for (size_t read : readsWithEdit[i])
		{
			hasEdit[read] = true;
		}
		for (size_t read = 0; read < hasEdit.size(); read++)
		{
			if (hasEdit[read]) continue;
			readsWithEdit[i+1].emplace_back(read);
		}
	}
	for (size_t i = 0; i < readsWithEdit.size(); i++)
	{
		std::sort(readsWithEdit[i].begin(), readsWithEdit[i].end());
	}
	std::vector<size_t> order;
	for (size_t i = 0; i < readsWithEdit.size(); i++)
	{
		order.emplace_back(i);
	}
	std::sort(order.begin(), order.end(), [&readsWithEdit](const size_t left, const size_t right)
	{
		return readsWithEdit[left] < readsWithEdit[right];
	});
	size_t clusterStart = 0;
	std::vector<std::vector<size_t>> readAlleles;
	readAlleles.resize(editsPerRead.size());
	for (size_t i = 1; i <= order.size(); i++)
	{
		if (i < order.size() && readsWithEdit[order[i-1]] == readsWithEdit[order[i]]) continue;
		if (i-clusterStart >= minEditsWithinCluster)
		{
			if (readsWithEdit[order[clusterStart]].size() >= minReadsWithinCluster && readsWithEdit[order[clusterStart]].size()+minReadsWithinCluster <= editsPerRead.size())
			{
				std::vector<size_t> positions;
				for (size_t j = clusterStart; j < i; j++)
				{
					positions.emplace_back(editPosition[order[j]]);
				}
				std::sort(positions.begin(), positions.end());
				size_t reallyDifferentLocations = 1;
				for (size_t j = 1; j < positions.size(); j++)
				{
					if (positions[j] < positions[j-1] + minDistanceBetweenEdits) continue;
					reallyDifferentLocations += 1;
				}
				bool hasGoodEdit = false;
				for (size_t j = clusterStart; j < i; j++)
				{
					if (isSNP(edits[order[j]]))
					{
						hasGoodEdit = true;
						break;
					}
					if (isBigIndel(edits[order[j]]) && !isMicrosatelliteOrHomopolymerIndel(edits[order[j]], refSequence))
					{
						hasGoodEdit = true;
						break;
					}
				}
				if (hasGoodEdit && reallyDifferentLocations >= minEditsWithinCluster)
				{
					std::vector<bool> hasAlt;
					hasAlt.resize(editsPerRead.size(), false);
					for (size_t read : readsWithEdit[order[clusterStart]])
					{
						hasAlt[read] = true;
						readAlleles[read].emplace_back(0);
					}
					for (size_t j = 0; j < editsPerRead.size(); j++)
					{
						if (hasAlt[j]) continue;
						readAlleles[j].emplace_back(1);
					}
				}
			}
		}
		clusterStart = i;
	}
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "edit linkage clusters " << readAlleles[0].size() << std::endl;
	std::map<std::vector<size_t>, size_t> allelesToCluster;
	std::vector<std::vector<size_t>> result;
	for (size_t i = 0; i < readAlleles.size(); i++)
	{
		if (allelesToCluster.count(readAlleles[i]) == 0)
		{
			allelesToCluster[readAlleles[i]] = result.size();
			result.emplace_back();
		}
		result[allelesToCluster.at(readAlleles[i])].emplace_back(i);
	}
	return result;
}

std::pair<size_t, size_t> getBubble(const GfaGraph& graph, const Node startNode)
{
	const std::pair<size_t, size_t> noBubble { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
	// std::numeric_limits<size_t>::max() cannot be used for undefined because Node can't contain it
	// break convention and use graph.numNodes() as undefined
	const Node undefinedNode { graph.numNodes(), true };
	if (graph.edges.count(startNode) == 0) return noBubble;
	if (graph.edges.at(startNode).size() != 2) return noBubble;
	Node otherSide = undefinedNode;
	size_t alleleOne = std::numeric_limits<size_t>::max();
	size_t alleleTwo = std::numeric_limits<size_t>::max();
	for (auto edge : graph.edges.at(startNode))
	{
		Node target = std::get<0>(edge);
		if (graph.edges.count(target) == 0) return noBubble;
		if (graph.edges.at(target).size() != 1) return noBubble;
		assert(graph.edges.count(reverse(target)) == 1);
		if (graph.edges.at(reverse(target)).size() != 1) return noBubble;
		if (otherSide == undefinedNode)
		{
			assert(alleleOne == std::numeric_limits<size_t>::max());
			alleleOne = target.id();
			otherSide = std::get<0>(*graph.edges.at(target).begin());
		}
		else
		{
			assert(alleleOne != std::numeric_limits<size_t>::max());
			assert(alleleTwo == std::numeric_limits<size_t>::max());
			alleleTwo = target.id();
			if (std::get<0>(*graph.edges.at(target).begin()) != otherSide) return noBubble;
		}
	}
	assert(otherSide != undefinedNode);
	assert(graph.edges.count(reverse(otherSide)) == 1);
	if (graph.edges.at(reverse(otherSide)).size() != 2) return noBubble;
	if (alleleOne == startNode.id()) return noBubble;
	if (alleleOne == otherSide.id()) return noBubble;
	if (alleleTwo == startNode.id()) return noBubble;
	if (alleleTwo == otherSide.id()) return noBubble;
	if (alleleOne == alleleTwo) return noBubble;
	return std::make_pair(alleleOne, alleleTwo);
}

std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>> getDBGvariants(const std::vector<OntLoop>& cluster, const std::string MBGPath, const std::string tmpPath)
{
	size_t minCoverage = std::max<size_t>(5, cluster.size()*0.05);
	std::string graphFile = tmpPath + "/tmpgraph.gfa";
	std::string pathsFile = tmpPath + "/tmppaths.gaf";
	std::string readsFile = tmpPath + "/tmpreads.fa";
	{
		std::ofstream file { readsFile };
		for (size_t i = 0; i < cluster.size(); i++)
		{
			file << ">" << i << std::endl;
			file << cluster[i].rawSequence << std::endl;
		}
	}
	std::string mbgCommand = MBGPath + " -i " + readsFile + " -o " + graphFile + " --output-sequence-paths " + pathsFile + " --error-masking=no -a 5 -u 5 --no-multiplex-cleaning -r 18446744073709551615 -R 18446744073709551615 -k 31 --only-local-resolve --do-unsafe-guesswork-resolutions 1> " + tmpPath + "/tmpstdout.txt 2> " + tmpPath + "/tmpstderr.txt";
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "MBG command:" << std::endl;
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << mbgCommand << std::endl;
	int runresult = system(mbgCommand.c_str());
	if (runresult != 0)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "MBG did not run successfully" << std::endl;
		std::abort();
	}
	GfaGraph graph;
	graph.loadFromFile(graphFile);
	std::vector<ReadPath> readPaths = loadReadPaths(pathsFile, graph);
	std::vector<phmap::flat_hash_set<size_t>> readsTouchingNode;
	readsTouchingNode.resize(graph.numNodes());
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		size_t read = std::stoull(readPaths[i].readName);
		for (Node node : readPaths[i].path)
		{
			readsTouchingNode[node.id()].insert(read);
		}
	}
	phmap::flat_hash_map<size_t, std::pair<size_t, size_t>> nodeBelongsToBubble;
	size_t countBubbles = 0;
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (readsTouchingNode[i].size() != cluster.size()) continue;
		std::pair<size_t, size_t> bubble = getBubble(graph, Node { i, true });
		if (bubble.first != std::numeric_limits<size_t>::max())
		{
			assert(bubble.second != std::numeric_limits<size_t>::max());
			assert(bubble.first != bubble.second);
			if (nodeBelongsToBubble.count(std::min(bubble.first, bubble.second)) == 1)
			{
				assert(nodeBelongsToBubble.count(std::max(bubble.first, bubble.second)) == 1);
			}
			else
			{
				nodeBelongsToBubble[std::min(bubble.first, bubble.second)] = std::make_pair(countBubbles, 0);
				nodeBelongsToBubble[std::max(bubble.first, bubble.second)] = std::make_pair(countBubbles, 1);
				countBubbles += 1;
			}
		}
		bubble = getBubble(graph, Node { i, false });
		if (bubble.first != std::numeric_limits<size_t>::max())
		{
			assert(bubble.second != std::numeric_limits<size_t>::max());
			assert(bubble.first != bubble.second);
			if (nodeBelongsToBubble.count(std::min(bubble.first, bubble.second)) == 1)
			{
				assert(nodeBelongsToBubble.count(std::max(bubble.first, bubble.second)) == 1);
			}
			else
			{
				nodeBelongsToBubble[std::min(bubble.first, bubble.second)] = std::make_pair(countBubbles, 0);
				nodeBelongsToBubble[std::max(bubble.first, bubble.second)] = std::make_pair(countBubbles, 1);
				countBubbles += 1;
			}
		}
	}
	std::vector<std::pair<std::vector<std::vector<size_t>>, size_t>> result;
	result.resize(countBubbles);
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i].second = i*1000; // fake positions, graph always has good good sites and doesn't falsely split a single variant into multiple edits
		result[i].first.resize(3); // always simple biallelic bubbles plus one for reads with no bubble allele because of sequencing errors
	}
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (nodeBelongsToBubble.count(i) == 0) continue;
		size_t bubbleIndex = nodeBelongsToBubble.at(i).first;
		size_t alleleIndex = nodeBelongsToBubble.at(i).second;
		for (size_t read : readsTouchingNode[i])
		{
			result[bubbleIndex].first[alleleIndex].emplace_back(read);
		}
	}
	for (size_t i = result.size()-1; i < result.size(); i--)
	{
		std::vector<bool> hasAllele;
		hasAllele.resize(cluster.size());
		assert(result[i].first.size() == 3);
		assert(result[i].first[2].size() == 0);
		bool good = true;
		if (result[i].first[0].size() < minCoverage) good = false;
		if (result[i].first[1].size() < minCoverage) good = false;
		for (size_t j = 0; j < 2; j++)
		{
			for (size_t read : result[i].first[j])
			{
				if (hasAllele[read]) good = false;
				hasAllele[read] = true;
			}
		}
		for (size_t j = 0; j < hasAllele.size(); j++)
		{
			if (hasAllele[j]) continue;
			result[i].first[2].emplace_back(j);
		}
		if (result[i].first[2].size() >= 3) good = false;
		// equality valid for breaking because truncation rounds down, only allow missing values once coverage is 20
		if (result[i].first[0].size() + result[i].first[1].size() <= cluster.size()*0.95) good = false;
		if (!good)
		{
			std::swap(result[i], result.back());
			result.pop_back();
		}
	}
	std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) { return left.second < right.second; });
	return result;
}

void callVariantsAndSplitRecursively(std::vector<std::vector<OntLoop>>& result, const std::vector<OntLoop>& cluster, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t numThreads, const std::string MBGPath, const std::string tmpPath, const bool correctBeforePhasing)
{
	if (cluster.size() <= 5)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "skip small cluster with size " << cluster.size() << std::endl;
		result.emplace_back(cluster);
		return;
	}
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "begin phasing cluster with size " << cluster.size() << std::endl;
	std::vector<std::string> sequences;
	if (correctBeforePhasing)
	{
		sequences = getSelfCorrectedLoops(cluster);
	}
	else
	{
		for (size_t i = 0; i < cluster.size(); i++)
		{
			sequences.emplace_back(cluster[i].rawSequence);
		}
	}
	auto startTime = getTime();
	auto refSequence = getConsensus(cluster, sequences, graph, pathStartClip, pathEndClip, numThreads);
	auto edits = getEditsForPhasing(refSequence, sequences, numThreads);
	auto phasableVariantInfo = getPhasableVariantInfoBiallelicAltRefSNPsBigIndels(edits, refSequence);
	auto endTime = getTime();
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "getting variants took " << formatTime(startTime, endTime) << std::endl;
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "got variant info of cluster with size " << cluster.size() << std::endl;
	std::vector<std::vector<size_t>> resultHere;
	std::vector<size_t> allIndices;
	resultHere = splitAndAdd(cluster, phasableVariantInfo);
	assert(resultHere.size() >= 1);
	if (resultHere.size() == 1)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "split by bubble alleles size " << cluster.size() << std::endl;
		resultHere = splitByBubbleAlleles(cluster, sequences, graph, pathStartClip, pathEndClip, numThreads);
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "bubble alleles splitted to " << resultHere.size() << " clusters:";
		for (size_t i = 0; i < resultHere.size(); i++)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << resultHere[i].size();
		}
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
	}
	if (resultHere.size() == 1)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "split by SNP MSA size " << cluster.size() << std::endl;
		auto SNPMSA = getSNPMSA(edits);
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << SNPMSA.size() << " sites in SNP MSA" << std::endl;
		resultHere.clear();
		resultHere = splitAndAdd(cluster, SNPMSA);
		if (resultHere.size() > 1)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "SNP MSA splitted to " << resultHere.size() << " clusters:";
			for (size_t i = 0; i < resultHere.size(); i++)
			{
				Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << resultHere[i].size();
			}
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
		}
	}
	if (resultHere.size() == 1)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "split by DBG variants size " << cluster.size() << std::endl;
		auto DBGvariants = getDBGvariants(cluster, MBGPath, tmpPath);
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << DBGvariants.size() << " sites in DBG variants" << std::endl;
		resultHere.clear();
		resultHere = splitAndAdd(cluster, DBGvariants);
		if (resultHere.size() > 1)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "DBG variants splitted to " << resultHere.size() << " clusters:";
			for (size_t i = 0; i < resultHere.size(); i++)
			{
				Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << resultHere[i].size();
			}
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
		}
	}
/*	if (resultHere.size() == 1)
	{
		std::cerr << "split by big indels size " << cluster.size() << std::endl;
		resultHere = splitByBigIndels(edits);
		std::cerr << "big indels splitted to " << resultHere.size() << " clusters:";
		for (size_t i = 0; i < resultHere.size(); i++)
		{
			std::cerr << " " << resultHere[i].size();
		}
		std::cerr << std::endl;
	}*/
/*	if (resultHere.size() == 1)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "split by edit linkage size " << cluster.size() << std::endl;
		resultHere = splitByEditLinkage(edits, refSequence);
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "edit linkage splitted to " << resultHere.size() << " clusters:";
		for (size_t i = 0; i < resultHere.size(); i++)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << resultHere[i].size();
		}
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
	}*/
/*	if (resultHere.size() == 1)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "split by linked minor alleles size " << cluster.size() << std::endl;
		resultHere = splitByLinkedMinorAlleles(edits, refSequence);
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "linked minor alleles splitted to " << resultHere.size() << " clusters:";
		for (size_t i = 0; i < resultHere.size(); i++)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << resultHere[i].size();
		}
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
	}*/
/*	if (resultHere.size() == 1 && correctBeforePhasing)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "split by SNP correlation size " << cluster.size() << std::endl;
		resultHere = splitBySNPCorrelation(edits, refSequence);
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "SNP correlation splitted to " << resultHere.size() << " clusters:";
		for (size_t i = 0; i < resultHere.size(); i++)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << resultHere[i].size();
		}
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
	}*/
	if (resultHere.size() == 1)
	{
		if (correctBeforePhasing)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "final cluster with size " << cluster.size() << " with reads:";
			for (const auto& read : cluster)
			{
				Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << read.readName << "_" << read.approxStart << "_" << read.approxEnd;
			}
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
			result.emplace_back(cluster);
			return;
		}
		else
		{
			callVariantsAndSplitRecursively(result, cluster, graph, pathStartClip, pathEndClip, numThreads, MBGPath, tmpPath, true);
			return;
		}
	}
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "split cluster with size " << cluster.size() << " corrected " << (correctBeforePhasing ? "yes" : "no")  << " into " << resultHere.size() << " clusters with reads:";
	for (const auto& read : cluster)
	{
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << read.readName << "_" << read.approxStart << "_" << read.approxEnd;
	}
	Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
	for (size_t i = 0; i < resultHere.size(); i++)
	{
		std::vector<OntLoop> filteredClusters;
		for (size_t j : resultHere[i])
		{
			filteredClusters.emplace_back(cluster[j]);
		}
		callVariantsAndSplitRecursively(result, filteredClusters, graph, pathStartClip, pathEndClip, numThreads, MBGPath, tmpPath, false);
	}
}

std::vector<std::vector<OntLoop>> editSplitClusters(const std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t numThreads, const std::string MBGPath, const std::string tmpPath)
{
	std::vector<std::vector<OntLoop>> result;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		if (clusters[i].size() <= 5)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "skip low coverage cluster with size " << clusters[i].size() << std::endl;
			result.emplace_back(clusters[i]);
			continue;
		}
		callVariantsAndSplitRecursively(result, clusters[i], graph, pathStartClip, pathEndClip, numThreads, MBGPath, tmpPath, false);
	}
	return result;
}
