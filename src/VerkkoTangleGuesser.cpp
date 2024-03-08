#include "fastqloader.h"
#include "KmerMatcher.h"
#include "VerkkoTangleGuesser.h"
#include "TangleGuesser.h"
#include "RibotinUtils.h"

std::string homopolymerCompress(std::string original)
{
	std::string result;
	result.push_back(original[0]);
	for (size_t i = 1; i < original.size(); i++)
	{
		if (original[i] != original[i-1]) result.push_back(original[i]);
	}
	return result;
}

std::vector<std::vector<std::string>> guessVerkkoRDNATangles(std::string verkkoBasePath, const std::vector<std::string>& referencePath)
{
	std::string graphFilePath = verkkoBasePath + "/assembly.homopolymer-compressed.gfa";
	if (!fileExists(graphFilePath))
	{
		std::cerr << "ERROR: could not find assembly graph file in the verkko assembly folder!" << std::endl;
		std::abort();
	}
	KmerMatcher matcher { 101 };
	for (auto file : referencePath)
	{
		FastQ::streamFastqFromFile(file, false, [&matcher](FastQ& fastq)
		{
			matcher.addReferenceKmers(homopolymerCompress(fastq.sequence));
		});
	}
	auto result = guessTangles(matcher, graphFilePath);
	return result;
}
