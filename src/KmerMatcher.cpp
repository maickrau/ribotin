#include "KmerMatcher.h"

size_t charToInt(char c)
{
	switch(c)
	{
		case 'a':
		case 'A':
			return 0;
		case 'c':
		case 'C':
			return 1;
		case 'g':
		case 'G':
			return 2;
		case 't':
		case 'T':
			return 3;
	}
	return 0;
}

std::string reverseComplementInt(std::string fwseq)
{
	std::string result { fwseq.rbegin(), fwseq.rend() };
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = 3-result[i];
	}
	return result;
}

KmerMatcher::KmerMatcher(size_t k) :
k(k),
s(k-10),
kmers()
{
}

void KmerMatcher::addReferenceKmers(std::string seq)
{
	for (size_t i = 0; i < seq.size(); i++)
	{
		seq[i] = charToInt(seq[i]);
	}
	auto revcomp = reverseComplementInt(seq);
	iterateSyncmersFast(seq, k, s, [this](size_t pos, size_t hash)
	{
		kmers.insert(hash);
	});
	iterateSyncmersFast(revcomp, k, s, [this](size_t pos, size_t hash)
	{
		kmers.insert(hash);
	});
}

size_t KmerMatcher::getMatchLength(std::string seq) const
{
	for (size_t i = 0; i < seq.size(); i++)
	{
		seq[i] = charToInt(seq[i]);
	}
	size_t lastMatch = 0;
	size_t matchSum = 0;
	iterateHashesFast(seq, k, kmers, [this, &matchSum, &lastMatch](size_t pos, size_t hash)
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
	return matchSum;
}
