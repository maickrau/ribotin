#include <string>
#include <iostream>
#include <phmap.h>
#include "KmerMatcher.h"

int main(int argc, char** argv)
{
	size_t k = std::stoi(argv[1]);
	size_t minMatch = std::stoi(argv[2]);
	std::string refFile { argv[3] };
	std::vector<std::string> readFiles;
	for (size_t i = 4; i < argc; i++)
	{
		readFiles.emplace_back(argv[i]);
	}
	iterateMatchingReads(refFile, readFiles, k, minMatch, [](const FastQ& seq)
	{
		std::cout << ">" << seq.seq_id << std::endl;
		std::cout << seq.sequence << std::endl;
	});
}
