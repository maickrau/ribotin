#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include <cxxopts.hpp>
#include "KmerMatcher.h"
#include "ClusterHandler.h"
#include "RibotinUtils.h"

void getKmers(std::string outputPrefix)
{
	std::ofstream file { outputPrefix + "/kmers.fa" };
	FastQ::streamFastqFromFile(outputPrefix + "/consensus.fa", false, [&file](FastQ& fastq)
	{
		file << ">consensus" << std::endl;
		file << fastq.sequence << std::endl;
	});
	std::ifstream variantgraph { outputPrefix + "/variant-graph.gfa" };
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
		("t", "Number of threads", cxxopts::value<size_t>()->default_value("1"))
		("x", "Preset parameters", cxxopts::value<std::string>())
		("sample-name", "Name of the sample added to all morph names", cxxopts::value<std::string>())
		("r,reference", "Reference used for recruiting reads (required)", cxxopts::value<std::string>())
		("orient-by-reference", "Rotate and possibly reverse complement the consensus to match the orientation of the given reference", cxxopts::value<std::string>())
		("approx-morphsize", "Approximate length of one morph")
		("k", "k-mer size", cxxopts::value<size_t>()->default_value("101"))
		("annotation-reference-fasta", "Lift over the annotations from given reference fasta+gff3 (requires liftoff)", cxxopts::value<std::string>())
		("annotation-gff3", "Lift over the annotations from given reference fasta+gff3 (requires liftoff)", cxxopts::value<std::string>())
		("morph-cluster-maxedit", "Maximum edit distance between two morphs to assign them into the same cluster", cxxopts::value<size_t>()->default_value("200"))
		("morph-recluster-minedit", "Minimum edit distance to recluster morphs", cxxopts::value<size_t>()->default_value("5"))
		("mbg", "MBG path", cxxopts::value<std::string>())
		("graphaligner", "GraphAligner path", cxxopts::value<std::string>())
	;
	std::string MBGPath;
	std::string GraphAlignerPath;
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
	if (params.count("x") == 1)
	{
		if (params["x"].as<std::string>() != "human")
		{
			std::cerr << "No preset \"" << params["x"].as<std::string>() << "\"" << std::endl;
			paramError = true;
		}
	}
	if (params.count("mbg") == 0)
	{
		std::cerr << "checking for MBG" << std::endl;
		int foundMBG = system("which MBG");
		if (foundMBG != 0)
		{
			std::cerr << "MBG not found" << std::endl;
			std::cerr << "MBG path (--mbg) is required" << std::endl;
			paramError = true;
		}
		else
		{
			MBGPath = "MBG";
		}
	}
	else
	{
		MBGPath = params["mbg"].as<std::string>();
	}
	if (params.count("graphaligner") == 1)
	{
		GraphAlignerPath = params["graphaligner"].as<std::string>();
	}
	if (params.count("graphaligner") == 0 && params.count("nano") >= 1)
	{
		std::cerr << "checking for GraphAligner" << std::endl;
		int foundGraphAligner = system("which GraphAligner");
		if (foundGraphAligner != 0)
		{
			std::cerr << "GraphAligner not found" << std::endl;
			std::cerr << "--graphaligner is required when using ultralong ONT reads" << std::endl;
			paramError = true;
		}
		else
		{
			GraphAlignerPath = "GraphAligner";
		}
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
	if (params.count("approx-morphsize") == 0 && params.count("x") == 0)
	{
		std::cerr << "Approximate size of one morph (--approx-morphsize) is required" << std::endl;
		paramError = true;
	}
	if (params.count("t") == 1 && params["t"].as<size_t>() == 0)
	{
		std::cerr << "number of threads can't be 0" << std::endl;
		paramError = true;
	}
	if (params.count("k") == 1 && params["k"].as<size_t>() < 31)
	{
		std::cerr << "k must be at least 31" << std::endl;
		paramError = true;
	}
	if (params.count("annotation-gff3") == 1 && params.count("annotation-reference-fasta") == 0)
	{
		std::cerr << "--annotation-reference-fasta is missing while --annotation-gff3 is used" << std::endl;
		paramError = true;
	}
	if (params.count("annotation-gff3") == 0 && params.count("annotation-reference-fasta") == 1)
	{
		std::cerr << "--annotation-gff3 is missing while --annotation-reference-fasta is used" << std::endl;
		paramError = true;
	}
	if (params.count("annotation-gff3") == 1 || params.count("annotation-reference-fasta") == 1)
	{
		std::cerr << "checking for liftoff" << std::endl;
		int foundMinimap2 = system("which liftoff");
		if (foundMinimap2 != 0)
		{
			std::cerr << "liftoff not found" << std::endl;
			paramError = true;
		}
	}
	if (paramError)
	{
		std::abort();
	}
	std::string refPath;
	std::vector<std::string> hifiReadPaths = params["i"].as<std::vector<std::string>>();
	std::vector<std::string> ontReadPaths;
	ClusterParams clusterParams;
	if (params.count("x") == 1)
	{
		if (params["x"].as<std::string>() == "human")
		{
			clusterParams.k = 101;
			clusterParams.maxClusterDifference = 200;
			clusterParams.minReclusterDistance = 5;
			clusterParams.maxResolveLength = 45000/5;
			clusterParams.orientReferencePath = std::string{RIBOTIN_TEMPLATE_PATH} + "/rDNA_one_unit.fasta";
			refPath = std::string{RIBOTIN_TEMPLATE_PATH} + "/chm13_rDNAs.fa";
		}
	}
	if (params.count("r") == 1) refPath = params["r"].as<std::string>();
	clusterParams.MBGPath = MBGPath;
	clusterParams.basePath = params["o"].as<std::string>();
	clusterParams.numThreads = 1;
	clusterParams.numThreads = params["t"].as<size_t>();
	clusterParams.maxClusterDifference = params["morph-cluster-maxedit"].as<size_t>();
	clusterParams.minReclusterDistance = params["morph-recluster-minedit"].as<size_t>();
	if (params.count("k") == 1) clusterParams.k = params["k"].as<size_t>();
	if (params.count("approx-morphsize") == 1) clusterParams.maxResolveLength = params["approx-morphsize"].as<size_t>()/5;
	if (params.count("nano") >= 1) ontReadPaths = params["nano"].as<std::vector<std::string>>();
	if (params.count("sample-name") == 1) clusterParams.namePrefix = params["sample-name"].as<std::string>();
	if (params.count("orient-by-reference") == 1) clusterParams.orientReferencePath = params["orient-by-reference"].as<std::string>();
	if (params.count("annotation-reference-fasta") == 1) clusterParams.annotationFasta = params["annotation-reference-fasta"].as<std::string>();
	if (params.count("annotation-gff3") == 1) clusterParams.annotationGff3 = params["annotation-gff3"].as<std::string>();
	clusterParams.GraphAlignerPath = GraphAlignerPath;
	std::cerr << "using reference from " << refPath << std::endl;
	std::cerr << "output folder: " << clusterParams.basePath << std::endl;
	std::filesystem::create_directories(clusterParams.basePath);
	std::cerr << "extracting HiFi/duplex reads" << std::endl;
	{
		std::ofstream readsfile { clusterParams.basePath + "/hifi_reads.fa" };
		iterateMatchingReads(refPath, hifiReadPaths, 101, 2000, [&readsfile](const FastQ& seq)
		{
			readsfile << ">" << seq.seq_id << std::endl;
			readsfile << seq.sequence << std::endl;
		});
	}
	std::cerr << "running" << std::endl;
	clusterParams.hifiReadPath = clusterParams.basePath + "/hifi_reads.fa";
	HandleCluster(clusterParams);
	if (ontReadPaths.size() > 0)
	{
		std::cerr << "extracting ultralong ONT reads" << std::endl;
		std::ofstream readsfile { clusterParams.basePath + "/ont_reads.fa" };
		getKmers(clusterParams.basePath);
		std::vector<std::string> kmerFiles;
		kmerFiles.push_back(refPath);
		kmerFiles.push_back(clusterParams.basePath + "/kmers.fa");
		size_t consensusLength = getSequenceLength(clusterParams.basePath + "/consensus.fa");
		std::cerr << "consensus length " << consensusLength << ", using " << consensusLength/2 << " as minimum ONT match length" << std::endl;
		iterateMatchingReads(kmerFiles, ontReadPaths, 21, consensusLength/2, [&readsfile](const FastQ& seq)
		{
			readsfile << ">" << seq.seq_id << std::endl;
			readsfile << seq.sequence << std::endl;
		});
		clusterParams.ontReadPath = clusterParams.basePath + "/ont_reads.fa";
	}
	if (ontReadPaths.size() > 0)
	{
		std::cerr << "start ultralong ONT analysis" << std::endl;
		std::cerr << "aligning ultralong ONT reads to allele graph" << std::endl;
		AlignONTReads(clusterParams.basePath, clusterParams.GraphAlignerPath, clusterParams.ontReadPath, clusterParams.basePath + "/allele-graph.gfa", clusterParams.basePath + "/ont-alns.gaf", clusterParams.numThreads);
		DoClusterONTAnalysis(clusterParams);
	}
}
