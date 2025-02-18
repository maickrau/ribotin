#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include <cxxopts.hpp>
#include "KmerMatcher.h"
#include "ClusterHandler.h"
#include "RibotinUtils.h"
#include "CommonParams.h"

void getKmers(std::string outputPrefix)
{
	std::ofstream file { outputPrefix + "/kmers.fa" };
	FastQ::streamFastqFromFile(outputPrefix + "/consensus.fa", false, [&file](FastQ& fastq)
	{
		file << ">consensus" << std::endl;
		file << fastq.sequence << std::endl;
	});
	std::ifstream variantgraph { outputPrefix + "/processed-graph.gfa" };
	while (variantgraph.good())
	{
		std::string line;
		getline(variantgraph, line);
		if (line.size() < 5 || line[0] != 'S') continue;
		std::stringstream sstr { line };
		std::string dummy, node, sequence;
		sstr >> dummy >> node >> sequence;
		file << ">node" << node << std::endl;
		file << sequence;
	}
}

int main(int argc, char** argv)
{
	std::cerr << "ribotin-ref version " << VERSION << std::endl;
	cxxopts::Options options { "ribotin-ref" };
	options.add_options()
		("h,help", "Print help")
		("v,version", "Print version")
		("i,in", "Input HiFi/duplex reads. Multiple files can be input with -i file1.fa -i file2.fa etc (required)", cxxopts::value<std::vector<std::string>>())
		("nano", "Input ultralong ONT reads. Multiple files can be input with --nano file1.fa --nano file2.fa etc", cxxopts::value<std::vector<std::string>>())
		("o,out", "Output folder", cxxopts::value<std::string>()->default_value("./result"))
		("r,reference", "Reference used for recruiting reads (required)", cxxopts::value<std::string>())
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
	if (!commonParams.hasGraphAligner() && params.count("nano") >= 1)
	{
		std::cerr << "--graphaligner is required when using ultralong ONT reads" << std::endl;
	}
	if (!commonParams.hasWinnowmap() && params.count("nano") >= 1)
	{
		std::cerr << "winnowmap not found, will not perform realignment of raw ONT loops to morph consensuses" << std::endl;
	}
	if (!commonParams.hasSamtools() && params.count("nano") >= 1)
	{
		std::cerr << "samtools not found, will not perform realignment of raw ONT loops to morph consensuses" << std::endl;
	}
	if (params.count("i") == 0)
	{
		std::cerr << "Input reads (-i) are required" << std::endl;
		paramError = true;
	}
	if (params.count("r") == 0 && params.count("x") == 0)
	{
		std::cerr << "Input reference (-r) is required" << std::endl;
		paramError = true;
	}
	if (paramError)
	{
		std::abort();
	}
	std::string refPath;
	std::vector<std::string> hifiReadPaths = params["i"].as<std::vector<std::string>>();
	std::vector<std::string> ontReadPaths;
	ClusterParams clusterParams;
	commonParams.addToClusterOptions(clusterParams);
	if (params.count("x") == 1)
	{
		if (params["x"].as<std::string>() == "human")
		{
			refPath = std::string{RIBOTIN_TEMPLATE_PATH} + "/chm13_rDNAs.fa";
		}
	}
	if (params.count("r") == 1) refPath = params["r"].as<std::string>();
	clusterParams.basePath = params["o"].as<std::string>();
	if (params.count("nano") >= 1) ontReadPaths = params["nano"].as<std::vector<std::string>>();
	std::cerr << "using reference from " << refPath << std::endl;
	std::cerr << "output folder: " << clusterParams.basePath << std::endl;
	std::filesystem::create_directories(clusterParams.basePath);
	// std::cerr << "extracting HiFi/duplex reads" << std::endl;
	// {
	// 	std::ofstream readsfile { clusterParams.basePath + "/hifi_reads.fa" };
	// 	iterateMatchingReads(refPath, hifiReadPaths, 101, 2000, [&readsfile](const FastQ& seq)
	// 	{
	// 		readsfile << ">" << seq.seq_id << std::endl;
	// 		readsfile << seq.sequence << std::endl;
	// 	});
	// }
	// std::cerr << "running" << std::endl;
	clusterParams.hifiReadPath = clusterParams.basePath + "/hifi_reads.fa";
	// HandleCluster(clusterParams);
	if (ontReadPaths.size() > 0)
	{
		// std::cerr << "extracting ultralong ONT reads" << std::endl;
		// std::ofstream readsfile { clusterParams.basePath + "/ont_reads.fa" };
		// getKmers(clusterParams.basePath);
		// std::vector<std::string> kmerFiles;
		// kmerFiles.push_back(refPath);
		// kmerFiles.push_back(clusterParams.basePath + "/kmers.fa");
		// size_t consensusLength = getSequenceLength(clusterParams.basePath + "/consensus.fa");
		// std::cerr << "consensus length " << consensusLength << ", using " << consensusLength/2 << " as minimum ONT match length" << std::endl;
		// iterateMatchingReads(kmerFiles, ontReadPaths, 21, consensusLength/2, [&readsfile](const FastQ& seq)
		// {
		// 	readsfile << ">" << seq.seq_id << std::endl;
		// 	readsfile << seq.sequence << std::endl;
		// });
		clusterParams.ontReadPath = clusterParams.basePath + "/ont_reads.fa";
	}
	if (ontReadPaths.size() > 0)
	{
		std::cerr << "start ultralong ONT analysis" << std::endl;
		// std::cerr << "aligning ultralong ONT reads to allele graph" << std::endl;
		// AlignONTReads(clusterParams.basePath, commonParams.GraphAlignerPath(), clusterParams.ontReadPath, clusterParams.basePath + "/processed-graph.gfa", clusterParams.basePath + "/ont-alns.gaf", clusterParams.numThreads);
		DoClusterONTAnalysis(clusterParams);
	}
	std::cerr << "done" << std::endl;
}
