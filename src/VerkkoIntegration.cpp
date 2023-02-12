#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include "VerkkoReadAssignment.h"
#include "ReadExtractor.h"
#include "ClusterHandler.h"

std::vector<std::string> getNodesFromFile(std::string filename)
{
	std::vector<std::string> result;
	std::ifstream file { filename };
	while (file.good())
	{
		std::string node;
		file >> node;
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

int main(int argc, char** argv)
{
	std::string verkkoBasePath = argv[1];
	std::string MBGPath = argv[2];
	std::vector<std::vector<std::string>> clusterNodes;
	size_t numClusters = argc - 3;
	std::cerr << "reading nodes per cluster" << std::endl;
	for (size_t i = 3; i < argc; i++)
	{
		clusterNodes.push_back(getNodesFromFile(std::string { argv[i] }));
	}
	std::cerr << "assigning reads per cluster" << std::endl;
	auto reads = getReadNamesPerCluster(verkkoBasePath, clusterNodes);
	for (size_t i = 0; i < reads.size(); i++)
	{
		std::cerr << "cluster " << i << " has " << reads[i].size() << " reads" << std::endl;
	}
	std::vector<std::string> readFileNames;
	for (size_t i = 0; i < numClusters; i++)
	{
		std::filesystem::create_directories("./cluster" + std::to_string(i));
		readFileNames.push_back("./cluster" + std::to_string(i) + "/reads.fa");
	}
	std::cerr << "extracting reads per cluster" << std::endl;
	splitReads(getRawReadFilenames(verkkoBasePath + "/verkko.yml"), reads, readFileNames);
	for (size_t i = 0; i < numClusters; i++)
	{
		std::cerr << "running cluster " << i << std::endl;
		HandleCluster("./cluster" + std::to_string(i), "./cluster" + std::to_string(i) + "/reads.fa", MBGPath);
	}
}
