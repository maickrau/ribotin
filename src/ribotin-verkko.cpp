#include <cassert>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include <regex>
#include <cxxopts.hpp>
#include "fastqloader.h"
#include "VerkkoReadAssignment.h"
#include "ReadExtractor.h"
#include "ClusterHandler.h"
#include "VerkkoTangleGuesser.h"
#include "KmerMatcher.h"
#include "RibotinUtils.h"
#include "CommonParams.h"

std::vector<std::string> getNodesFromFile(std::string filename)
{
	std::vector<std::string> result;
	std::ifstream file { filename };
	while (file.good())
	{
		std::string node;
		file >> node;
		if (!file.good()) break;
		if (node.size() > 1 && node.back() == ',') node.pop_back();
		result.push_back(node);
	}
	return result;
}

std::vector<std::string> getRawReadFilenames(std::string configPath, std::string readTypeLine)
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
			if (line.find(readTypeLine) != std::string::npos)
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

bool checkAssemblyHasNanopore(std::string configPath)
{
	// terrible way, but works for now
	std::ifstream file { configPath };
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (line.find("withONT:") == std::string::npos) continue;
		if (line.find("True") != std::string::npos) return true;
		return false;
	}
	return false;
}

void writeNodes(std::string filename, const std::vector<std::string>& nodes)
{
	std::ofstream file { filename };
	for (auto node : nodes)
	{
		file << node << std::endl;
	}
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
		("c,tangles", "Input files for node tangles. Multiple files may be inputed with -c file1.txt -c file2.txt ... (required)", cxxopts::value<std::vector<std::string>>())
		("guess-tangles-using-reference", "Guess the rDNA tangles using k-mer matches to given reference sequence (required)", cxxopts::value<std::vector<std::string>>())
		("do-ul", "Do ONT read analysis (yes/no/as_verkko)", cxxopts::value<std::string>()->default_value("as_verkko"))
		("ul-tmp-folder", "Temporary folder for ultralong ONT read analysis", cxxopts::value<std::string>()->default_value("./tmp"))
	;
	CommonParams commonParams;
	commonParams.addParamsToOptions(options);
	auto params = options.parse(argc, argv);
	if (params.count("v") == 1)
	{
		std::cerr << "Version: " << VERSION << std::endl;
		std::exit(0);
	}
	if (params.count("h") == 1)
	{
		std::cerr << options.help() << std::endl;
		std::cerr << "Valid presets for -x: human" << std::endl;
		std::exit(0);
	}
	bool paramError = false;
	bool commonParamsParseSuccessful = commonParams.parseParamsAndPrintErrors(params);
	if (!commonParamsParseSuccessful) paramError = true;
	if (params.count("do-ul") == 1 && (params["do-ul"].as<std::string>() != "no" && params["do-ul"].as<std::string>() != "yes" && params["do-ul"].as<std::string>() != "as_verkko"))
	{
		std::cerr << "Unknown option for --do-ul: \"" << params["do-ul"].as<std::string>() << "\"" << std::endl;
		paramError = true;
	}
	if ((params["do-ul"].as<std::string>() == "yes" || params["do-ul"].as<std::string>() == "as_verkko") && !commonParams.hasGraphAligner())
	{
		std::cerr << "--graphaligner is required when using --do-ul=yes or --do-ul=as_verkko" << std::endl;
		paramError = true;
	}
	if (params.count("i") == 0)
	{
		std::cerr << "Input verkko folder (-i) is required" << std::endl;
		paramError = true;
	}
	if (params.count("c") == 0 && params.count("guess-tangles-using-reference") == 0 && params.count("x") == 0)
	{
		std::cerr << "Either node tangles (-c) or reference used for guessing (--guess-tangles-using-reference) are required" << std::endl;
		paramError = true;
	}
	if (params.count("c") == 1 && params.count("guess-tangles-using-reference") >= 1)
	{
		std::cerr << "Only one of node tangles (-c) or reference used for guessing (--guess-tangles-using-reference) can be used" << std::endl;
		paramError = true;
	}
	if (paramError)
	{
		std::abort();
	}
	std::string verkkoBasePath = params["i"].as<std::string>();
	if (!fileExists(verkkoBasePath + "/verkko.yml"))
	{
		std::cerr << "ERROR: could not find configuration file in the verkko assembly folder!" << std::endl;
		std::abort();
	}
	bool doUL = false;
	if (params["do-ul"].as<std::string>() == "yes")
	{
		doUL = true;
	}
	if (params["do-ul"].as<std::string>() == "as_verkko" && checkAssemblyHasNanopore(verkkoBasePath + "/verkko.yml"))
	{
		doUL = true;
	}
	if (doUL && !checkAssemblyHasNanopore(verkkoBasePath + "/verkko.yml"))
	{
		std::cerr << "Assembly did not use ultralong ONT reads, cannot do ultralong ONT analysis." << std::endl;
		std::cerr << "Try running with --do-ul=no to skip ultralong ONT analysis, or rerun Verkko with nanopore reads." << std::endl;
		std::abort();
	}
	std::string outputPrefix = params["o"].as<std::string>();
	std::cerr << "do UL analysis: " << (doUL ? "yes" : "no") << std::endl;
	std::cerr << "output prefix: " << outputPrefix << std::endl;
	std::string ulTmpFolder = params["ul-tmp-folder"].as<std::string>();
	std::vector<std::vector<std::string>> tangleNodes;
	if (params.count("x") == 1)
	{
		if (params["x"].as<std::string>() == "human")
		{
			if (params.count("c") == 0 && params.count("guess-tangles-using-reference") == 0)
			{
				std::cerr << "guessing tangles" << std::endl;
				std::string refPath = std::string{RIBOTIN_TEMPLATE_PATH} + "/chm13_rDNAs.fa";
				std::cerr << "using reference from " << refPath << std::endl;
				tangleNodes = guessVerkkoRDNATangles(verkkoBasePath, std::vector<std::string>{refPath});
				std::cerr << "resulted in " << tangleNodes.size() << " tangles" << std::endl;
			}
		}
	}
	if (params.count("c") >= 1)
	{
		std::vector<std::string> tangleNodeFiles = params["c"].as<std::vector<std::string>>();
		std::cerr << "reading nodes per tangle" << std::endl;
		for (size_t i = 0; i < tangleNodeFiles.size(); i++)
		{
			tangleNodes.push_back(getNodesFromFile(tangleNodeFiles[i]));
		}
	}
	else if (params.count("guess-tangles-using-reference") >= 1)
	{
		std::cerr << "guessing tangles using reference:" << std::endl;
		std::vector<std::string> refPaths = params["guess-tangles-using-reference"].as<std::vector<std::string>>();
		for (auto path : refPaths)
		{
			std::cerr << path << std::endl;
		}
		tangleNodes = guessVerkkoRDNATangles(verkkoBasePath, refPaths);
		std::cerr << "resulted in " << tangleNodes.size() << " tangles" << std::endl;
	}
	size_t numTangles = tangleNodes.size();
	if (numTangles == 0)
	{
		std::cerr << "ERROR: No rDNA tangles found!" << std::endl;
		std::abort();
	}
	std::cerr << "assigning reads per tangle" << std::endl;
	auto reads = getReadNamesPerTangle(verkkoBasePath, tangleNodes);
	for (size_t i = 0; i < reads.size(); i++)
	{
		std::cerr << "tangle " << i << " has " << reads[i].size() << " hifi reads" << std::endl;
	}
	std::vector<std::string> hifiReadPaths = getRawReadFilenames(verkkoBasePath + "/verkko.yml", "HIFI_READS");
	for (size_t i = 0; i < numTangles; i++)
	{
		std::filesystem::create_directories(outputPrefix + std::to_string(i));
		writeNodes(outputPrefix + std::to_string(i) + "/nodes.txt", tangleNodes[i]);
	}
	ClusterParams clusterParams;
	commonParams.addToClusterOptions(clusterParams);
	std::vector<std::string> ulReadNames;
	if (doUL)
	{
		ulReadNames = getRawReadFilenames(verkkoBasePath + "/verkko.yml", "ONT_READS");
	}
	runMultipleTangles(clusterParams, outputPrefix, numTangles, reads, hifiReadPaths, doUL, ulTmpFolder, ulReadNames);
	std::cerr << "done" << std::endl;
}
