#include <cassert>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include <regex>
#include <cxxopts.hpp>
#include "fastqloader.h"
#include "ReadExtractor.h"
#include "ClusterHandler.h"
#include "HifiasmIntegration.h"
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
	std::cerr << "ribotin-hifiasm version " << VERSION << std::endl;
	cxxopts::Options options { "ribotin-hifiasm" };
	options.add_options()
		("h,help", "Print help")
		("v,version", "Print version")
		("a,assembly", "Input hifiasm assembly prefix (required)", cxxopts::value<std::string>())
		("i,hifi", "Input hifi reads (required)", cxxopts::value<std::vector<std::string>>())
		("nano", "Input ultralong ONT reads", cxxopts::value<std::vector<std::string>>())
		("o,out", "Output folder prefix", cxxopts::value<std::string>()->default_value("./result"))
		("c,tangles", "Input files for node tangles. Multiple files may be inputed with -c file1.txt -c file2.txt ... (required)", cxxopts::value<std::vector<std::string>>())
		("guess-tangles-using-reference", "Guess the rDNA tangles using k-mer matches to given reference sequence (required)", cxxopts::value<std::vector<std::string>>())
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
	if (params.count("assembly") == 0)
	{
		std::cerr << "Input hifiasm assembly is required" << std::endl;
		paramError = true;
	}
	if (params.count("hifi") == 0 || params["hifi"].as<std::vector<std::string>>().size() == 0)
	{
		std::cerr << "Input hifi reads are required" << std::endl;
		paramError = true;
	}
	if (!commonParams.hasGraphAligner() && params.count("nano") >= 1 && params["nano"].as<std::vector<std::string>>().size() >= 1)
	{
		std::cerr << "--graphaligner is required when using --nano" << std::endl;
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
	std::string hifiasmBasePath = params["assembly"].as<std::string>();
	bool doUL = false;
	if (params.count("nano") >= 1 && params["nano"].as<std::vector<std::string>>().size() >= 1) doUL = true;
	std::string outputPrefix = params["o"].as<std::string>();
	std::vector<std::string> hifiReadPaths = params["hifi"].as<std::vector<std::string>>();
	std::vector<std::string> ontReadPaths;
	if (params.count("nano") >= 1) ontReadPaths = params["nano"].as<std::vector<std::string>>();
	std::cerr << "assembly prefix: " << hifiasmBasePath << std::endl;
	std::cerr << "searching for hifiasm graph file" << std::endl;
	std::string hifiasmGraphFile = getHifiasmGraphFileName(hifiasmBasePath);
	std::cerr << "using hifiasm graph file: " << hifiasmGraphFile << std::endl;
	std::cerr << "hifi read paths:";
	for (auto path : hifiReadPaths)
	{
		std::cerr << " " << path;
	}
	std::cerr << std::endl;
	std::cerr << "do UL analysis: " << (doUL ? "yes" : "no") << std::endl;
	if (doUL)
	{
		std::cerr << "ultralong ONT read paths:";
		for (auto path : ontReadPaths)
		{
			std::cerr << " " << path;
		}
		std::cerr << std::endl;
	}
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
				tangleNodes = guessHifiasmRDNATangles(hifiasmGraphFile, std::vector<std::string>{refPath});
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
			std::cerr << "reading from " << tangleNodeFiles[i] << std::endl;
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
		tangleNodes = guessHifiasmRDNATangles(hifiasmGraphFile, refPaths);
		std::cerr << "resulted in " << tangleNodes.size() << " tangles" << std::endl;
	}
	size_t numTangles = tangleNodes.size();
	if (numTangles == 0)
	{
		std::cerr << "ERROR: No rDNA tangles found!" << std::endl;
		std::abort();
	}
	std::cerr << "assigning reads per tangle" << std::endl;
	auto reads = getHifiasmReadNamesPerTangle(hifiasmGraphFile, tangleNodes);
	for (size_t i = 0; i < reads.size(); i++)
	{
		std::cerr << "tangle " << i << " has " << reads[i].size() << " hifi reads" << std::endl;
	}
	for (size_t i = 0; i < numTangles; i++)
	{
		std::filesystem::create_directories(outputPrefix + std::to_string(i));
		writeNodes(outputPrefix + std::to_string(i) + "/nodes.txt", tangleNodes[i]);
	}
	ClusterParams clusterParams;
	commonParams.addToClusterOptions(clusterParams);
	runMultipleTangles(clusterParams, outputPrefix, numTangles, reads, hifiReadPaths, doUL, ulTmpFolder, ontReadPaths);
	std::cerr << "done" << std::endl;
}
