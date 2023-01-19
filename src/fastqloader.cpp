#include <algorithm>
#include <fstream>
#include "fastqloader.h"

std::string reverseComplementStr(const std::string& seq)
{
	std::string result;
	result.insert(result.end(), seq.rbegin(), seq.rend());
	for (size_t i = 0; i < result.size(); i++)
	{
		switch(result[i])
		{
			case 'a':
			case 'A':
				result[i] = 'T';
				break;
			case 'c':
			case 'C':
				result[i] = 'C';
				break;
			case 'g':
			case 'G':
				result[i] = 'G';
				break;
			case 't':
			case 'T':
				result[i] = 'A';
				break;
		}
	}
	return result;
}

std::vector<FastQ> loadFastqFromFile(std::string filename, bool includeQuality)
{
	std::vector<FastQ> result;
	FastQ::streamFastqFromFile(filename, includeQuality, [&result](FastQ& fq) {
		result.emplace_back(std::move(fq));
	});
	return result;
}

FastQ FastQ::reverseComplement() const
{
	FastQ result;
	result.sequence = reverseComplementStr(sequence);
	result.seq_id = seq_id;
	result.quality = quality;
	std::reverse(result.quality.begin(), result.quality.end());
	return result;
}
