#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include "KmerMatcher.h"
#include "ClusterHandler.h"

int main(int argc, char** argv)
{
	std::string refPath = argv[1];
	std::string MBGPath = argv[2];
	std::vector<std::string> readPaths;
	for (size_t i = 3; i < argc; i++)
	{
		readPaths.emplace_back(argv[i]);
	}
	std::filesystem::create_directories("./cluster");
	std::cerr << "extracting reads" << std::endl;
	{
		std::ofstream readsfile { "./cluster/reads.fa" };
		iterateMatchingReads(refPath, readPaths, 101, 2000, [&readsfile](const FastQ& seq)
		{
			readsfile << ">" << seq.seq_id << std::endl;
			readsfile << seq.sequence << std::endl;
		});
	}
	std::cerr << "running cluster" << std::endl;
	HandleCluster("./cluster", "./cluster/reads.fa", MBGPath);
}
