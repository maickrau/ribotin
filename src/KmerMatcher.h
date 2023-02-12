#ifndef KmerMatcher_h
#define KmerMatcher_h

#include <string>
#include <phmap.h>
#include "fastqloader.h"
#include "FastHasher.h"

size_t charToInt(char c);

template <typename F>
void findSyncmerPositions(const std::string& sequence, size_t kmerSize, size_t smerSize, std::vector<std::tuple<size_t, uint64_t>>& smerOrder, F callback)
{
	if (sequence.size() < kmerSize) return;
	assert(smerSize <= kmerSize);
	size_t windowSize = kmerSize - smerSize + 1;
	assert(windowSize >= 1);
	FastHasher fwkmerHasher { smerSize };
	for (size_t i = 0; i < smerSize; i++)
	{
		fwkmerHasher.addChar(sequence[i]);
	}
	auto thisHash = fwkmerHasher.hash();
	smerOrder.emplace_back(0, thisHash);
	for (size_t i = 1; i < windowSize; i++)
	{
		size_t seqPos = smerSize+i-1;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-smerSize]);
		uint64_t hash = fwkmerHasher.hash();
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > hash) smerOrder.pop_back();
		smerOrder.emplace_back(i, hash);
	}
	if ((std::get<0>(smerOrder.front()) == 0) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == windowSize-1))
	{
		callback(0);
	}
	for (size_t i = windowSize; smerSize+i-1 < sequence.size(); i++)
	{
		size_t seqPos = smerSize+i-1;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-smerSize]);
		uint64_t hash = fwkmerHasher.hash();
		// even though pop_front is used it turns out std::vector is faster than std::deque ?!
		// because pop_front is O(w), but it is only called in O(1/w) fraction of loops
		// so the performace penalty of pop_front does not scale with w!
		// and std::vector's speed in normal, non-popfront operation outweighs the slow pop_front
		while (smerOrder.size() > 0 && std::get<0>(smerOrder.front()) <= i - windowSize) smerOrder.erase(smerOrder.begin());
		while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > hash) smerOrder.pop_back();
		smerOrder.emplace_back(i, hash);
		if ((std::get<0>(smerOrder.front()) == i-windowSize+1) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == i))
		{
			callback(i-windowSize+1);
		}
	}
}

template <typename F>
void iterateSyncmersFast(const std::string& seq, const size_t k, const size_t s, F callback)
{
	static std::vector<std::tuple<size_t, uint64_t>> smerOrder;
	std::vector<size_t> startPoses;
	findSyncmerPositions(seq, k, s, smerOrder, [&startPoses](size_t pos) { startPoses.push_back(pos); });
	FastHasher fwkmerHasher { k };
	for (size_t i = 0; i < k; i++)
	{
		fwkmerHasher.addChar(seq[i]);
	}
	size_t pospos = 0;
	if (startPoses.size() == 0) return;
	if (startPoses[0] == 0)
	{
		callback(0, fwkmerHasher.hash());
		pospos += 1;
	}
	for (size_t i = k; i < seq.size(); i++)
	{
		assert(pospos < startPoses.size());
		assert(startPoses[pospos] >= i-k+1);
		fwkmerHasher.addChar(seq[i]);
		fwkmerHasher.removeChar(seq[i-k]);
		if (startPoses[pospos] > i-k+1) continue;
		auto thisHash = fwkmerHasher.hash();
		callback(i-k+1, thisHash);
		pospos += 1;
		if (pospos == startPoses.size()) break;
	}
	assert(pospos == startPoses.size());
}

template <typename F>
void iterateHashesFast(const std::string& seq, const size_t k, const phmap::flat_hash_set<uint64_t>& hashes, F callback)
{
	FastHasher fwkmerHasher { k };
	for (size_t i = 0; i < k; i++)
	{
		fwkmerHasher.addChar(seq[i]);
	}
	for (size_t i = k; i < seq.size(); i++)
	{
		auto thisHash = fwkmerHasher.hash();
		if (hashes.count(thisHash) == 1)
		{
			callback(i-k, thisHash);
		}
		fwkmerHasher.addChar(seq[i]);
		fwkmerHasher.removeChar(seq[i-k]);
	}
}

template <typename F>
void iterateMatchingReads(std::string refFile, const std::vector<std::string>& readFiles, size_t k, size_t minMatch, F callback)
{
	size_t s = k-10;
	phmap::flat_hash_set<size_t> kmers;
	size_t count = 0;
	FastQ::streamFastqFromFile(refFile, false, [k, s, &kmers, &count](FastQ& fastq)
	{
		auto revcomp = fastq.reverseComplement();
		for (size_t i = 0; i < fastq.sequence.size(); i++)
		{
			fastq.sequence[i] = charToInt(fastq.sequence[i]);
		}
		for (size_t i = 0; i < fastq.sequence.size(); i++)
		{
			revcomp.sequence[i] = charToInt(revcomp.sequence[i]);
		}
		iterateSyncmersFast(fastq.sequence, k, s, [&kmers, &count](size_t pos, size_t hash)
		{
			kmers.insert(hash);
			count += 1;
		});
		iterateSyncmersFast(revcomp.sequence, k, s, [&kmers, &count](size_t pos, size_t hash)
		{
			kmers.insert(hash);
			count += 1;
		});
	});
	for (const auto& file : readFiles)
	{
		FastQ::streamFastqFromFile(file, false, [k, s, minMatch, &kmers, callback](FastQ& fastq)
		{
			for (size_t i = 0; i < fastq.sequence.size(); i++)
			{
				fastq.sequence[i] = charToInt(fastq.sequence[i]);
			}
			bool include = false;
			size_t lastMatch = 0;
			size_t matchSum = 0;
			iterateHashesFast(fastq.sequence, k, kmers, [&matchSum, &lastMatch, k](size_t pos, size_t hash)
			{
				if (pos-lastMatch < k)
				{
					matchSum += pos-lastMatch;
					lastMatch = pos;
				}
				else
				{
					matchSum += k;
					lastMatch = pos;
				}
			});
			if (matchSum >= minMatch) include = true;
			if (!include) return;
			for (size_t i = 0; i < fastq.sequence.size(); i++)
			{
				fastq.sequence[i] = "ACGT"[fastq.sequence[i]];
			}
			callback(fastq);
		});
	}
}

#endif
