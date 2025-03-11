#ifndef ConsensusHelper_h
#define ConsensusHelper_h

#include <thread>
#include <vector>
#include <string>
#include "ClusterMisc.h"
#include "phmap.h"

phmap::flat_hash_map<uint64_t, size_t> getRefKmers(const std::string_view& ref, const size_t k);
std::vector<std::tuple<size_t, size_t, std::string, size_t>> getEditsRec(const std::string_view& refSeq, const std::string_view& querySeq, const size_t k);
std::vector<std::pair<size_t, size_t>> getKmerAnchors(const std::string_view& ref, const phmap::flat_hash_map<uint64_t, size_t>& refKmers, const std::string_view& query, const size_t k);
std::string getConsensusSequence(const std::vector<Node>& consensusPath, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip);
std::vector<std::string> getSelfCorrectedLoops(const std::vector<OntLoop>& cluster);
std::string polishConsensus(const std::string rawConsensus, const std::vector<std::string>& seqs, const size_t numThreads);
std::vector<Node> getConsensusPath(const std::vector<OntLoop>& rawPaths, const GfaGraph& graph);
std::string getSequence(const std::vector<Node>& nodes, const std::vector<std::string>& nodeSeqs, const std::vector<std::string>& revCompNodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges);
std::vector<std::vector<size_t>> getAllPresentKmerChain(const std::string consensusSeq, const std::vector<std::string>& loopSequences, const size_t k);
void removeRepeatKmers(std::vector<std::pair<uint64_t, uint64_t>>& matches);
void removeRepeatKmersEitherRepeat(std::vector<std::pair<uint64_t, uint64_t>>& matches);
void removeNonDiagonalKmers(std::vector<std::pair<size_t, size_t>>& kmerMatches);
std::vector<std::pair<size_t, size_t>> getIncreasingChain(const std::vector<std::pair<size_t, size_t>>& kmerMatches);

template <typename F>
void iterateEdits(const std::string& rawConsensus, const phmap::flat_hash_map<uint64_t, size_t>& refKmers, const std::string& sequence, const size_t k, const size_t threadIndex, const size_t readIndex, F callback)
{
	auto anchors = getKmerAnchors(std::string_view { rawConsensus }, refKmers, std::string_view { sequence }, k);
	if (anchors.size() == 0) return;
	for (size_t i = 1; i < anchors.size(); i++)
	{
		assert(anchors[i-1].first < anchors[i].first);
		assert(anchors[i-1].second < anchors[i].second);
		assert(anchors[i].first+k <= rawConsensus.size());
		assert(anchors[i].second+k <= sequence.size());
		assert(rawConsensus.substr(anchors[i].first, k) == sequence.substr(anchors[i].second, k));
		assert(rawConsensus.substr(anchors[i-1].first, k) == sequence.substr(anchors[i-1].second, k));
		if (anchors[i].first - anchors[i-1].first == anchors[i].second - anchors[i-1].second && anchors[i].first - anchors[i-1].first < k) continue;
		std::string_view refSubstr { rawConsensus.data()+anchors[i-1].first, anchors[i].first - anchors[i-1].first + k };
		std::string_view querySubstr { sequence.data()+anchors[i-1].second, anchors[i].second - anchors[i-1].second + k };
		auto edits = getEditsRec(refSubstr, querySubstr, 25);
		for (auto edit : edits)
		{
			assert(std::get<1>(edit) >= std::get<0>(edit));
			assert(std::get<1>(edit) < refSubstr.size());
			std::get<0>(edit) += anchors[i-1].first;
			std::get<1>(edit) += anchors[i-1].first;
			assert(std::get<1>(edit) < rawConsensus.size());
			std::get<3>(edit) += anchors[i-1].second;
			callback(threadIndex, readIndex, anchors[0].first, anchors.back().first+k, anchors[0].second, anchors.back().second+k, edit);
		}
	}
}

template <typename F>
void iterateEdits(const std::string& rawConsensus, const std::vector<OntLoop>& rawLoops, const size_t numThreads, F callback)
{
	const size_t k = 31;
	auto refKmers = getRefKmers(std::string_view { rawConsensus }, k);
	std::atomic<size_t> nextIndex;
	nextIndex = 0;
	std::vector<std::thread> threads;
	for (size_t threadIndex = 0; threadIndex < numThreads; threadIndex++)
	{
		threads.emplace_back([&rawConsensus, &rawLoops, &nextIndex, &refKmers, threadIndex, callback]()
		{
			while (true)
			{
				size_t readIndex = nextIndex++;
				if (readIndex >= rawLoops.size()) break;
				iterateEdits(rawConsensus, refKmers, rawLoops[readIndex].rawSequence, k, threadIndex, readIndex, callback);
			}
		});
	}
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
}

template <typename F>
void iterateEdits(const std::string& rawConsensus, const std::vector<std::string>& rawLoops, const size_t numThreads, F callback)
{
	const size_t k = 31;
	auto refKmers = getRefKmers(std::string_view { rawConsensus }, k);
	std::atomic<size_t> nextIndex;
	nextIndex = 0;
	std::vector<std::thread> threads;
	for (size_t threadIndex = 0; threadIndex < numThreads; threadIndex++)
	{
		threads.emplace_back([&rawConsensus, &rawLoops, &nextIndex, &refKmers, threadIndex, callback]()
		{
			while (true)
			{
				size_t readIndex = nextIndex++;
				if (readIndex >= rawLoops.size()) break;
				iterateEdits(rawConsensus, refKmers, rawLoops[readIndex], k, threadIndex, readIndex, callback);
			}
		});
	}
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
}

#endif
