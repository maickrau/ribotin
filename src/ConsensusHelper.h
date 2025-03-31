#ifndef ConsensusHelper_h
#define ConsensusHelper_h

#include <thread>
#include <vector>
#include <string>
#include "ClusterMisc.h"
#include "phmap.h"

phmap::flat_hash_map<uint64_t, size_t> getRefKmers(const std::string_view& ref, const size_t k);
std::vector<std::tuple<size_t, size_t, std::string, size_t>> getEdits(const std::string_view& refSequence, const std::string_view& querySequence, size_t maxEdits);
std::vector<std::pair<size_t, size_t>> getKmerAnchors(const std::string_view& ref, const phmap::flat_hash_map<uint64_t, size_t>& refKmers, const std::string_view& query, const size_t k);
std::string getConsensusSequence(const std::vector<Node>& consensusPath, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip);
std::vector<std::string> getSelfCorrectedLoops(const std::vector<OntLoop>& cluster);
std::string polishConsensus(const std::string rawConsensus, const std::vector<std::string>& seqs, const size_t numThreads);
std::vector<Node> getConsensusPath(const std::vector<OntLoop>& rawPaths, const GfaGraph& graph);
std::string getSequence(const std::vector<Node>& nodes, const std::vector<std::string>& nodeSeqs, const std::vector<std::string>& revCompNodeSeqs, const std::unordered_map<Node, phmap::flat_hash_set<std::tuple<Node, size_t, size_t>>>& edges);
void removeRepeatKmers(std::vector<std::pair<uint64_t, uint64_t>>& matches);
void removeRepeatKmersEitherRepeat(std::vector<std::pair<uint64_t, uint64_t>>& matches);
void removeNonDiagonalKmers(std::vector<std::pair<size_t, size_t>>& kmerMatches);
std::vector<std::pair<size_t, size_t>> getIncreasingChain(const std::vector<std::pair<size_t, size_t>>& kmerMatches);
std::vector<std::tuple<size_t, size_t, size_t>> getMatchRegions(const std::string& refSequence, const std::string& querySequence);
std::vector<std::vector<size_t>> getMatchBases(const std::string& consensusSeq, const std::vector<std::string>& loopSequences);
std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<std::vector<std::tuple<size_t, size_t, size_t>>>> getRegionsConservedInAllSequences(const std::string& consensusSeq, const std::vector<std::string>& loopSequences);

template <typename F>
void iterateEdits(const std::string& rawConsensus, const std::string& sequence, const size_t threadIndex, const size_t readIndex, F callback)
{
	auto anchors = getMatchRegions(rawConsensus, sequence);
	if (anchors.size() == 0) return;
	for (size_t i = 1; i < anchors.size(); i++)
	{
		assert(std::get<0>(anchors[i-1]) < std::get<0>(anchors[i]));
		assert(std::get<1>(anchors[i-1]) < std::get<1>(anchors[i]));
		assert(std::get<0>(anchors[i])+std::get<2>(anchors[i]) <= rawConsensus.size());
		assert(std::get<1>(anchors[i])+std::get<2>(anchors[i]) <= sequence.size());
		assert(rawConsensus.substr(std::get<0>(anchors[i]), std::get<2>(anchors[i])) == sequence.substr(std::get<1>(anchors[i]), std::get<2>(anchors[i])));
		assert(rawConsensus.substr(std::get<0>(anchors[i-1]), std::get<2>(anchors[i-1])) == sequence.substr(std::get<1>(anchors[i-1]), std::get<2>(anchors[i-1])));
		if (std::get<0>(anchors[i]) - std::get<0>(anchors[i-1]) == std::get<1>(anchors[i]) - std::get<1>(anchors[i-1]) && std::get<0>(anchors[i]) - std::get<0>(anchors[i-1]) < std::get<2>(anchors[i-1])) continue;
		std::string_view refSubstr { rawConsensus.data()+std::get<0>(anchors[i-1]), std::get<0>(anchors[i]) - std::get<0>(anchors[i-1]) + std::get<2>(anchors[i]) };
		std::string_view querySubstr { sequence.data()+std::get<1>(anchors[i-1]), std::get<1>(anchors[i]) - std::get<1>(anchors[i-1]) + std::get<2>(anchors[i]) };
		auto edits = getEdits(refSubstr, querySubstr, (refSubstr.size()+querySubstr.size())*2);
		for (auto edit : edits)
		{
			assert(std::get<1>(edit) >= std::get<0>(edit));
			assert(std::get<1>(edit) < refSubstr.size());
			std::get<0>(edit) += std::get<0>(anchors[i-1]);
			std::get<1>(edit) += std::get<0>(anchors[i-1]);
			assert(std::get<1>(edit) < rawConsensus.size());
			std::get<3>(edit) += std::get<1>(anchors[i-1]);
			callback(threadIndex, readIndex, std::get<0>(anchors[0]), std::get<0>(anchors.back())+std::get<2>(anchors.back()), std::get<1>(anchors[0]), std::get<1>(anchors.back())+std::get<2>(anchors.back()), edit);
		}
	}
}

template <typename F>
void iterateEdits(const std::string& rawConsensus, const std::vector<OntLoop>& rawLoops, const size_t numThreads, F callback)
{
	std::atomic<size_t> nextIndex;
	nextIndex = 0;
	std::vector<std::thread> threads;
	for (size_t threadIndex = 0; threadIndex < numThreads; threadIndex++)
	{
		threads.emplace_back([&rawConsensus, &rawLoops, &nextIndex, threadIndex, callback]()
		{
			while (true)
			{
				size_t readIndex = nextIndex++;
				if (readIndex >= rawLoops.size()) break;
				iterateEdits(rawConsensus, rawLoops[readIndex].rawSequence, threadIndex, readIndex, callback);
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
	std::atomic<size_t> nextIndex;
	nextIndex = 0;
	std::vector<std::thread> threads;
	for (size_t threadIndex = 0; threadIndex < numThreads; threadIndex++)
	{
		threads.emplace_back([&rawConsensus, &rawLoops, &nextIndex, threadIndex, callback]()
		{
			while (true)
			{
				size_t readIndex = nextIndex++;
				if (readIndex >= rawLoops.size()) break;
				iterateEdits(rawConsensus, rawLoops[readIndex], threadIndex, readIndex, callback);
			}
		});
	}
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
}

#endif
