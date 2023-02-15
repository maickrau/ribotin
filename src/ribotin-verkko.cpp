#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include <cxxopts.hpp>
#include "VerkkoReadAssignment.h"
#include "ReadExtractor.h"
#include "ClusterHandler.h"
#include "VerkkoClusterGuesser.h"

std::vector<std::string> getNodesFromFile(std::string filename)
{
	std::vector<std::string> result;
	std::ifstream file { filename };
	while (file.good())
	{
		std::string node;
		file >> node;
		if (node.size() > 1 && node.back() == ',') node.pop_back();
		result.push_back(node);
	}
	return result;
}

std::vector<std::string> getRawReadFilenames(std::string configPath)
{
	// terrible way, but works for now
	std::vector<std::string> result;
	std::ifstream file { configPath };
	bool nowHifiFiles = false;
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (line.find(":") != std::string::npos)
		{
			nowHifiFiles = false;
			if (line.find("HIFI_READS") != std::string::npos)
			{
				nowHifiFiles = true;
			}
		}
		if (nowHifiFiles)
		{
			if (line.size() >= 5 && line[1] == '-')
			{
				result.push_back(line.substr(4, line.size()-5));
			}
		}
	}
	return result;
}

void writeNodes(std::string filename, const std::vector<std::string>& nodes)
{
	std::ofstream file { filename };
	for (auto node : nodes)
	{
		file << node << std::endl;
	}
}

void writeName(std::string filename, const std::string& name)
{
	std::ofstream file { filename };
	file << name;
}

int main(int argc, char** argv)
{
	std::cerr << "ribotin-verkko version " << VERSION << std::endl;
	cxxopts::Options options { "ribotin-verkko" };
	options.add_options()
		("h,help", "Print help")
		("v,version", "Print version")
		("i,in", "Input verkko folder (required)", cxxopts::value<std::string>())
		("o,out", "Output folder prefix", cxxopts::value<std::string>()->default_value("./result"))
		("c,cluster", "Input files for node clusters. Multiple files may be inputed with -c file1.txt -c file2.txt ... (required)", cxxopts::value<std::vector<std::string>>())
		("guess-clusters-using-reference", "Guess the rDNA clusters using k-mer matches to given reference sequence (required)", cxxopts::value<std::vector<std::string>>())
		("mbg", "MBG path (required)", cxxopts::value<std::string>())
		("k", "k-mer size", cxxopts::value<size_t>()->default_value("101"))
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
		std::cerr << "Input verkko folder (-i) is required" << std::endl;
		paramError = true;
	}
	if (params.count("c") == 0 && params.count("guess-clusters-using-reference") == 0)
	{
		std::cerr << "Either node clusters (-c) or reference used for guessing (--guess-clusters-using-reference) are required" << std::endl;
		paramError = true;
	}
	if (params.count("c") == 1 && params.count("guess-clusters-using-reference") == 1)
	{
		std::cerr << "Only one of node clusters (-c) or reference used for guessing (--guess-clusters-using-reference) can be used" << std::endl;
		paramError = true;
	}
	if (params.count("mbg") == 0)
	{
		std::cerr << "MBG path (--mbg) is required" << std::endl;
		paramError = true;
	}
	if (params.count("k") == 1 && params["k"].as<size_t>() < 31)
	{
		std::cerr << "k must be at least 31" << std::endl;
		paramError = true;
	}
	if (paramError)
	{
		std::abort();
	}
	std::string verkkoBasePath = params["i"].as<std::string>();
	std::string MBGPath = params["mbg"].as<std::string>();
	std::string outputPrefix = params["o"].as<std::string>();
	size_t k = params["k"].as<size_t>();
	std::cerr << "output prefix: " << outputPrefix << std::endl;
	std::vector<std::vector<std::string>> clusterNodes;
	if (params.count("c") == 1)
	{
		std::vector<std::string> clusterNodeFiles = params["c"].as<std::vector<std::string>>();
		std::cerr << "reading nodes per cluster" << std::endl;
		for (size_t i = 0; i < clusterNodeFiles.size(); i++)
		{
			clusterNodes.push_back(getNodesFromFile(clusterNodeFiles[i]));
		}
	}
	else
	{
		std::cerr << "guessing clusters" << std::endl;
		clusterNodes = guessVerkkoRDNAClusters(verkkoBasePath, params["guess-clusters-using-reference"].as<std::vector<std::string>>());
		std::cerr << "resulted in " << clusterNodes.size() << " clusters" << std::endl;
	}
	size_t numClusters = clusterNodes.size();
	std::cerr << "assigning reads per cluster" << std::endl;
	auto reads = getReadNamesPerCluster(verkkoBasePath, clusterNodes);
	for (size_t i = 0; i < reads.size(); i++)
	{
		std::cerr << "cluster " << i << " has " << reads[i].size() << " reads" << std::endl;
	}
	std::vector<std::string> readFileNames;
	for (size_t i = 0; i < numClusters; i++)
	{
		std::filesystem::create_directories(outputPrefix + std::to_string(i));
		readFileNames.push_back(outputPrefix + std::to_string(i) + "/reads.fa");
		writeNodes(outputPrefix + std::to_string(i) + "/nodes.txt", clusterNodes[i]);
	}
	std::cerr << "extracting reads per cluster" << std::endl;
	splitReads(getRawReadFilenames(verkkoBasePath + "/verkko.yml"), reads, readFileNames);
	std::vector<size_t> clustersWithoutReads;
	for (size_t i = 0; i < numClusters; i++)
	{
		if (reads[i].size() == 0)
		{
			std::cerr << "WARNING: cluster " << i << " has no reads, skipping" << std::endl;
			clustersWithoutReads.push_back(i);
			continue;
		}
		std::cerr << "running cluster " << i << " in folder " << outputPrefix + std::to_string(i) << std::endl;
		HandleCluster(outputPrefix + std::to_string(i), outputPrefix + std::to_string(i) + "/reads.fa", MBGPath, k);
	}
	if (clustersWithoutReads.size() > 0)
	{
		std::cerr << "WARNING: some clusters did not have any reads assigned:";
		for (auto cluster : clustersWithoutReads) std::cerr << " " << cluster;
		std::cerr << ", something likely went wrong." << std::endl;
	}
}
