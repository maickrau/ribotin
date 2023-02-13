#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include <cxxopts.hpp>
#include "KmerMatcher.h"
#include "ClusterHandler.h"

int main(int argc, char** argv)
{
	std::cerr << "rdnaConsensus-ref version " << VERSION << std::endl;
	cxxopts::Options options { "rdnaConsensus-ref" };
	options.add_options()
		("h,help", "Print help")
		("v,version", "Print version")
		("i,in", "Input reads. Multiple files can be input with -i file1.fa -i file2.fa etc (required)", cxxopts::value<std::vector<std::string>>())
		("o,out", "Output folder", cxxopts::value<std::string>()->default_value("./result"))
		("mbg", "MBG path (required)", cxxopts::value<std::string>())
		("r,ref", "Reference used for recruiting reads (required)", cxxopts::value<std::string>())
	;
	auto params = options.parse(argc, argv);
	if (params.count("v") == 1)
	{
		std::cerr << "Version: " << VERSION << std::endl;
		std::exit(0);
	}
	if (params.count("h") == 1)
	{
		std::cerr << options.help() << std::endl;
		std::exit(0);
	}
	bool paramError = false;
	if (params.count("i") == 0)
	{
		std::cerr << "Input reads (-i) are required" << std::endl;
		paramError = true;
	}
	if (params.count("r") == 0)
	{
		std::cerr << "Input reference (-r) is required" << std::endl;
		paramError = true;
	}
	if (params.count("mbg") == 0)
	{
		std::cerr << "MBG path (--mbg) is required" << std::endl;
		paramError = true;
	}
	if (paramError)
	{
		std::abort();
	}
	std::string refPath = params["r"].as<std::string>();
	std::string MBGPath = params["mbg"].as<std::string>();
	std::string outputPath = params["o"].as<std::string>();
	std::vector<std::string> readPaths = params["i"].as<std::vector<std::string>>();
	std::cerr << "output folder: " << outputPath << std::endl;
	std::filesystem::create_directories(outputPath);
	std::cerr << "extracting reads" << std::endl;
	{
		std::ofstream readsfile { outputPath + "/reads.fa" };
		iterateMatchingReads(refPath, readPaths, 101, 2000, [&readsfile](const FastQ& seq)
		{
			readsfile << ">" << seq.seq_id << std::endl;
			readsfile << seq.sequence << std::endl;
		});
	}
	std::cerr << "running" << std::endl;
	HandleCluster(outputPath, outputPath + "/reads.fa", MBGPath);
}
