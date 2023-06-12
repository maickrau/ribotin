#include "RibotinUtils.h"
#include "fastqloader.h"

std::vector<std::string> split(const std::string& str, char separator)
{
	std::vector<std::string> result;
	size_t lastBreak = 0;
	for (size_t i = 0; i < str.size(); i++)
	{
		if (str[i] != separator) continue;
		if (i > lastBreak) result.emplace_back(str.begin()+lastBreak, str.begin()+i);
		lastBreak = i+1;
	}
	if (lastBreak < str.size()) result.emplace_back(str.begin()+lastBreak, str.end());
	return result;
}

size_t getSequenceLength(const std::string& filename)
{
	size_t result = 0;
	FastQ::streamFastqFromFile(filename, false, [&result](const FastQ& seq){ result = seq.sequence.size(); });
	return result;
}

