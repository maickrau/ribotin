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
	std::cerr << "ribotin-ref version " << VERSION << std::endl;
	cxxopts::Options options { "ribotin-ref" };
	options.add_options()
		("h,help", "Print help")
		("v,version", "Print version")
		("i,in", "Input reads. Multiple files can be input with -i file1.fa -i file2.fa etc (required)", cxxopts::value<std::vector<std::string>>())
		("o,out", "Output folder", cxxopts::value<std::string>()->default_value("./result"))
		("mbg", "MBG path (required)", cxxopts::value<std::string>())
		("r,reference", "Reference used for recruiting reads (required)", cxxopts::value<std::string>())
		("orient-by-reference", "Rotate and possibly reverse complement the consensus to match the orientation of the given reference", cxxopts::value<std::string>())
		("k", "k-mer size", cxxopts::value<size_t>()->default_value("101"))
		("annotation-reference-fasta", "Lift over the annotations from given reference fasta+gff3 (requires liftoff)", cxxopts::value<std::string>())
		("annotation-gff3", "Lift over the annotations from given reference fasta+gff3 (requires liftoff)", cxxopts::value<std::string>())
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
	std::string refPath = params["r"].as<std::string>();
	std::vector<std::string> readPaths = params["i"].as<std::vector<std::string>>();
	ClusterParams clusterParams;
	clusterParams.MBGPath = params["mbg"].as<std::string>();
	clusterParams.basePath = params["o"].as<std::string>();
	clusterParams.k = params["k"].as<size_t>();
	if (params.count("orient-by-reference") == 1) clusterParams.orientReferencePath = params["orient-by-reference"].as<std::string>();
	if (params.count("annotation-reference-fasta") == 1) clusterParams.annotationFasta = params["annotation-reference-fasta"].as<std::string>();
	if (params.count("annotation-gff3") == 1) clusterParams.annotationGff3 = params["annotation-gff3"].as<std::string>();
	std::cerr << "output folder: " << clusterParams.basePath << std::endl;
	std::filesystem::create_directories(clusterParams.basePath);
	std::cerr << "extracting reads" << std::endl;
	{
		std::ofstream readsfile { clusterParams.basePath + "/reads.fa" };
		iterateMatchingReads(refPath, readPaths, 101, 2000, [&readsfile](const FastQ& seq)
		{
			readsfile << ">" << seq.seq_id << std::endl;
			readsfile << seq.sequence << std::endl;
		});
	}
	std::cerr << "running" << std::endl;
	clusterParams.readPath = clusterParams.basePath + "/reads.fa";
	HandleCluster(clusterParams);
}
