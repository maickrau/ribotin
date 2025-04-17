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

bool fileExists(std::string filename)
{
	std::ifstream file { filename };
	return file.good();
}

std::chrono::time_point<std::chrono::steady_clock> getTime()
{
	return std::chrono::steady_clock::now();
}

std::string formatTime(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end)
{
	size_t milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
	return std::to_string(milliseconds / 1000) + "," + std::to_string(milliseconds % 1000) + " s";
}

std::string nameWithoutTags(const std::string& rawname)
{
	std::string result = rawname;
	auto spacepos = result.find(' ');
	if (spacepos != std::string::npos) result = result.substr(0, spacepos);
	auto tabpos = result.find('\t')
	if (tabpos != std::string::npos) result = result.substr(0, tabpos);
	return result;
}
