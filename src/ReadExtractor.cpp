#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "ReadExtractor.h"
#include "fastqloader.h"

std::string nameWithoutTags(const std::string& rawname)
{
	auto spacepos = rawname.find(' ');
	if (spacepos == std::string::npos) return rawname;
	return rawname.substr(0, spacepos);
}

void splitReads(std::vector<std::string> readFiles, const std::vector<std::vector<std::string>>& readsPerCluster, const std::vector<std::string>& outputFileNames)
{
	assert(readsPerCluster.size() == outputFileNames.size());
	std::unordered_map<std::string, std::unordered_set<size_t>> readBelongsToCluster;
	for (size_t i = 0; i < readsPerCluster.size(); i++)
	{
		for (const auto& read : readsPerCluster[i])
		{
			readBelongsToCluster[nameWithoutTags(read)].insert(i);
		}
	}
	std::vector<std::ofstream> outfiles;
	outfiles.reserve(outputFileNames.size());
	for (size_t i = 0; i < outputFileNames.size(); i++)
	{
		outfiles.emplace_back(outputFileNames[i]);
	}
	for (const auto& file : readFiles)
	{
		std::cerr << "extracting reads from " << file << std::endl;
		FastQ::streamFastqFromFile(std::string{ file }, false, [&readBelongsToCluster, &outfiles](FastQ& fastq)
		{
			if (readBelongsToCluster.count(nameWithoutTags(fastq.seq_id)) == 0) return;
			for (auto cluster : readBelongsToCluster.at(nameWithoutTags(fastq.seq_id)))
			{
				assert(cluster < outfiles.size());
				outfiles[cluster] << ">" << fastq.seq_id << std::endl;
				outfiles[cluster] << fastq.sequence << std::endl;
			}
		});
	}
}
