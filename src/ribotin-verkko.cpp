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

void writeName(std::string filename, const std::string& name)
{
	std::ofstream file { filename };
	file << name;
}

void getKmers(std::string outputPrefix, size_t numTangles, std::string outputFile)
{
	std::ofstream file { outputFile };
	for (size_t i = 0; i < numTangles; i++)
	{
		FastQ::streamFastqFromFile(outputPrefix + std::to_string(i) + "/consensus.fa", false, [&file, i](FastQ& fastq)
		{
			file << ">consensus" << i << std::endl;
			file << fastq.sequence << std::endl;
		});
		std::ifstream variantgraph { outputPrefix + std::to_string(i) + "/variant-graph.gfa" };
		while (variantgraph.good())
		{
			std::string line;
			getline(variantgraph, line);
			if (line.size() < 5 || line[0] != 'S') continue;
			std::stringstream sstr { line };
			std::string dummy, node, sequence;
			sstr >> dummy >> node >> sequence;
			file << ">graph" << i << "node" << node << std::endl;
			file << sequence;
		}
	}
}

void mergeGraphs(std::string outputPrefix, size_t numTangles, std::string outputFile)
{
	std::ofstream file { outputFile };
	for (size_t i = 0; i < numTangles; i++)
	{
		std::ifstream variantgraph { outputPrefix + std::to_string(i) + "/allele-graph.gfa" };
		while (variantgraph.good())
		{
			std::string line;
			getline(variantgraph, line);
			std::stringstream sstr { line };
			std::string linetype;
			sstr >> linetype;
			if (linetype == "S")
			{
				std::string node, sequence;
				sstr >> node >> sequence;
				node = "graph" + std::to_string(i) + "node" + node;
				file << "S\t" << node << "\t" << sequence << std::endl;
			}
			else if (linetype == "L")
			{
				std::string fromnode, fromorient, tonode, toorient, overlap;
				sstr >> fromnode >> fromorient >> tonode >> toorient >> overlap;
				fromnode = "graph" + std::to_string(i) + "node" + fromnode;
				tonode = "graph" + std::to_string(i) + "node" + tonode;
				file << "L\t" << fromnode << "\t" << fromorient << "\t" << tonode << "\t" << toorient << "\t" << overlap << std::endl;
			}
		}
	}
}

size_t getTangleNum(const std::string& pathstr)
{
	size_t firstGraph = pathstr.find("graph");
	size_t firstNode = pathstr.find("node");
	assert(firstGraph != std::string::npos);
	assert(firstNode != std::string::npos);
	assert(firstNode > firstGraph+5);
	size_t result = std::stoull(pathstr.substr(firstGraph+5, firstNode-firstGraph-5));
	return result;
}

void splitAlignmentsPerTangle(std::string outputPrefix, size_t numTangles, std::string rawGafFile)
{
	std::unordered_map<std::string, size_t> uniqueMatch;
	{
		std::ifstream file { rawGafFile };
		while (file.good())
		{
			std::string line;
			getline(file, line);
			if (!file.good()) break;
			auto parts = split(line, '\t');
			std::string readname = parts[0];
			size_t readstart = std::stoull(parts[2]);
			size_t readend = std::stoull(parts[3]);
			std::string pathstr = parts[5];
			size_t mapq = std::stoull(parts[11]);
			if (mapq < 20) continue;
			assert(readend > readstart);
			if (readend - readstart < 20000) continue;
			size_t tangleNum = getTangleNum(pathstr);
			assert(tangleNum < numTangles);
			if (uniqueMatch.count(readname) == 0) uniqueMatch[readname] = tangleNum;
			if (uniqueMatch.count(readname) == 1 && uniqueMatch.at(readname) != tangleNum) uniqueMatch[readname] = std::numeric_limits<size_t>::max();
		}
	}
	std::vector<std::ofstream> perTangleFiles;
	for (size_t i = 0; i < numTangles; i++)
	{
		perTangleFiles.emplace_back(outputPrefix + std::to_string(i) + "/ont-alns.gaf");
	}
	std::ifstream file { rawGafFile };
	std::regex extraGraphRemover { "([<>])graph\\d+node" };
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (!file.good()) break;
		auto parts = split(line, '\t');
		std::string readname = parts[0];
		size_t readstart = std::stoull(parts[2]);
		size_t readend = std::stoull(parts[3]);
		std::string pathstr = parts[5];
		size_t mapq = std::stoull(parts[11]);
		if (mapq < 20) continue;
		assert(readend > readstart);
		if (readend - readstart < 20000) continue;
		if (uniqueMatch.at(readname) == std::numeric_limits<size_t>::max()) continue;
		size_t tangleNum = uniqueMatch.at(readname);
		assert(tangleNum < perTangleFiles.size());
		line = std::regex_replace(line, extraGraphRemover, "$1");
		perTangleFiles[tangleNum] << line << std::endl;
	}
}

size_t medianConsensusLength(const std::string& outputPrefix, size_t numTangles)
{
	std::vector<size_t> lengths;
	for (size_t i = 0; i < numTangles; i++)
	{
		size_t consensusLength = getSequenceLength(outputPrefix + std::to_string(i) + "/consensus.fa");
		lengths.push_back(consensusLength);
	}
	std::sort(lengths.begin(), lengths.end());
	return lengths[lengths.size()/2];
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
		("orient-by-reference", "Rotate and possibly reverse complement the consensus to match the orientation of the given reference", cxxopts::value<std::string>())
		("mbg", "MBG path", cxxopts::value<std::string>())
		("graphaligner", "GraphAligner path", cxxopts::value<std::string>())
		("do-ul", "Do ultralong ONT read analysis (requires GraphAligner)")
		("ul-tmp-folder", "Temporary folder for ultralong ONT read analysis", cxxopts::value<std::string>()->default_value("./tmp"))
		("k", "k-mer size", cxxopts::value<size_t>()->default_value("101"))
		("annotation-reference-fasta", "Lift over the annotations from given reference fasta+gff3 (requires liftoff)", cxxopts::value<std::string>())
		("annotation-gff3", "Lift over the annotations from given reference fasta+gff3 (requires liftoff)", cxxopts::value<std::string>())
		("morph-cluster-maxedit", "Maximum edit distance between two morphs to assign them into the same cluster", cxxopts::value<size_t>()->default_value("200"))
		("morph-recluster-minedit", "Minimum edit distance to recluster morphs", cxxopts::value<size_t>()->default_value("5"))
		("t", "Number of threads", cxxopts::value<size_t>()->default_value("1"))
		("approx-morphsize", "Approximate length of one morph", cxxopts::value<size_t>()->default_value("45000"))
		("sample-name", "Name of the sample added to all morph names", cxxopts::value<std::string>())
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
		std::exit(0);
	}
	bool paramError = false;
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
	if (params.count("do-ul") == 1 && params.count("graphaligner") == 0)
	{
		std::cerr << "checking for GraphAligner" << std::endl;
		int foundGraphAligner = system("which GraphAligner");
		if (foundGraphAligner != 0)
		{
			std::cerr << "GraphAligner not found" << std::endl;
			std::cerr << "--graphaligner is required when using --ul-ont" << std::endl;
			paramError = true;
		}
		else
		{
			GraphAlignerPath = "GraphAligner";
		}
	}
	if (params.count("i") == 0)
	{
		std::cerr << "Input verkko folder (-i) is required" << std::endl;
		paramError = true;
	}
	if (params.count("c") == 0 && params.count("guess-tangles-using-reference") == 0)
	{
		std::cerr << "Either node tangles (-c) or reference used for guessing (--guess-tangles-using-reference) are required" << std::endl;
		paramError = true;
	}
	if (params.count("c") == 1 && params.count("guess-tangles-using-reference") == 1)
	{
		std::cerr << "Only one of node tangles (-c) or reference used for guessing (--guess-tangles-using-reference) can be used" << std::endl;
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
		int foundLiftoff = system("which liftoff");
		if (foundLiftoff != 0)
		{
			std::cerr << "liftoff not found" << std::endl;
			paramError = true;
		}
	}
	if (paramError)
	{
		std::abort();
	}
	bool doUL = params.count("do-ul") == 1;
	std::string verkkoBasePath = params["i"].as<std::string>();
	if (doUL && !checkAssemblyHasNanopore(verkkoBasePath + "/verkko.yml"))
	{
		std::cerr << "Assembly did not use ultralong ONT reads, cannot do ultralong ONT analysis." << std::endl;
		std::cerr << "Try running without --do-ul to skip ultralong ONT analysis, or rerun Verkko with nanopore reads." << std::endl;
		std::abort();
	}
	std::string outputPrefix = params["o"].as<std::string>();
	size_t k = params["k"].as<size_t>();
	std::cerr << "output prefix: " << outputPrefix << std::endl;
	std::string orientReferencePath;
	std::string annotationFasta;
	std::string annotationGff3;
	std::string ulTmpFolder = params["ul-tmp-folder"].as<std::string>();
	std::string sampleName;
	size_t numThreads = params["t"].as<size_t>();
	size_t maxClusterDifference = params["morph-cluster-maxedit"].as<size_t>();
	size_t minReclusterDistance = params["morph-recluster-minedit"].as<size_t>();
	size_t maxResolveLength = params["approx-morphsize"].as<size_t>()/5;
	if (params.count("sample-name") == 1) sampleName = params["sample-name"].as<std::string>();
	if (params.count("orient-by-reference") == 1) orientReferencePath = params["orient-by-reference"].as<std::string>();
	if (params.count("annotation-reference-fasta") == 1) annotationFasta = params["annotation-reference-fasta"].as<std::string>();
	if (params.count("annotation-gff3") == 1) annotationGff3 = params["annotation-gff3"].as<std::string>();
	std::vector<std::vector<std::string>> tangleNodes;
	if (params.count("c") >= 1)
	{
		std::vector<std::string> tangleNodeFiles = params["c"].as<std::vector<std::string>>();
		std::cerr << "reading nodes per tangle" << std::endl;
		for (size_t i = 0; i < tangleNodeFiles.size(); i++)
		{
			tangleNodes.push_back(getNodesFromFile(tangleNodeFiles[i]));
		}
	}
	else
	{
		std::cerr << "guessing tangles" << std::endl;
		tangleNodes = guessVerkkoRDNATangles(verkkoBasePath, params["guess-tangles-using-reference"].as<std::vector<std::string>>());
		std::cerr << "resulted in " << tangleNodes.size() << " tangles" << std::endl;
	}
	size_t numTangles = tangleNodes.size();
	std::cerr << "assigning reads per tangle" << std::endl;
	auto reads = getReadNamesPerCluster(verkkoBasePath, tangleNodes);
	for (size_t i = 0; i < reads.size(); i++)
	{
		std::cerr << "tangle " << i << " has " << reads[i].size() << " hifi reads" << std::endl;
	}
	std::vector<std::string> readFileNames;
	for (size_t i = 0; i < numTangles; i++)
	{
		std::filesystem::create_directories(outputPrefix + std::to_string(i));
		readFileNames.push_back(outputPrefix + std::to_string(i) + "/hifi_reads.fa");
		writeNodes(outputPrefix + std::to_string(i) + "/nodes.txt", tangleNodes[i]);
	}
	std::cerr << "extracting HiFi/duplex reads per tangle" << std::endl;
	splitReads(getRawReadFilenames(verkkoBasePath + "/verkko.yml", "HIFI_READS"), reads, readFileNames);
	std::vector<size_t> tanglesWithoutReads;
	for (size_t i = 0; i < numTangles; i++)
	{
		if (reads[i].size() == 0)
		{
			std::cerr << "WARNING: tangle " << i << " has no HiFi/duplex reads, skipping" << std::endl;
			tanglesWithoutReads.push_back(i);
			continue;
		}
		ClusterParams clusterParams;
		clusterParams.maxClusterDifference = maxClusterDifference;
		clusterParams.minReclusterDistance = minReclusterDistance;
		clusterParams.namePrefix = "tangle" + std::to_string(i);
		if (sampleName != "") clusterParams.namePrefix = sampleName + "_" + clusterParams.namePrefix;
		clusterParams.basePath = outputPrefix + std::to_string(i);
		clusterParams.hifiReadPath = outputPrefix + std::to_string(i) + "/hifi_reads.fa";
		if (doUL)
		{
			clusterParams.ontReadPath = outputPrefix + std::to_string(i) + "/ont_reads.fa";
		}
		clusterParams.MBGPath = MBGPath;
		clusterParams.GraphAlignerPath = GraphAlignerPath;
		clusterParams.k = k;
		clusterParams.orientReferencePath = orientReferencePath;
		clusterParams.annotationFasta = annotationFasta;
		clusterParams.annotationGff3 = annotationGff3;
		clusterParams.maxResolveLength = maxResolveLength;
		std::cerr << "running tangle " << i << " in folder " << outputPrefix + std::to_string(i) << std::endl;
		HandleCluster(clusterParams);
	}
	if (doUL)
	{
		std::filesystem::create_directories(ulTmpFolder);
		std::cerr << "getting kmers from tangles" << std::endl;
		getKmers(outputPrefix, numTangles, ulTmpFolder + "/rdna_kmers.fa");
		std::cerr << "extracting ultralong ONT reads" << std::endl;
		auto fileNames = getRawReadFilenames(verkkoBasePath + "/verkko.yml", "ONT_READS");
		std::ofstream readsfile { ulTmpFolder + "/ont_reads.fa" };
		size_t consensusLength = medianConsensusLength(outputPrefix, numTangles);
		std::cerr << "median consensus length " << consensusLength << ", using " << consensusLength/2 << " as minimum ONT match length" << std::endl;
		iterateMatchingReads(ulTmpFolder + "/rdna_kmers.fa", fileNames, 21, consensusLength/2, [&readsfile](const FastQ& seq)
		{
			readsfile << ">" << seq.seq_id << std::endl;
			readsfile << seq.sequence << std::endl;
		});
		std::cerr << "merging allele graphs" << std::endl;
		mergeGraphs(outputPrefix, numTangles, ulTmpFolder + "/merged-allele-graph.gfa");
		std::cerr << "aligning ONT reads" << std::endl;
		AlignONTReads(ulTmpFolder, GraphAlignerPath, ulTmpFolder + "/ont_reads.fa", ulTmpFolder + "/merged-allele-graph.gfa", ulTmpFolder + "/ont-alns.gaf", numThreads);
		std::cerr << "splitting ONTs per tangle" << std::endl;
		splitAlignmentsPerTangle(outputPrefix, numTangles, ulTmpFolder + "/ont-alns.gaf");
		for (size_t i = 0; i < numTangles; i++)
		{
			if (reads[i].size() == 0)
			{
				continue;
			}
			std::cerr << "running tangle " << i << std::endl;
			ClusterParams clusterParams;
			clusterParams.maxClusterDifference = maxClusterDifference;
			clusterParams.minReclusterDistance = minReclusterDistance;
			clusterParams.namePrefix = "tangle" + std::to_string(i);
			if (sampleName != "") clusterParams.namePrefix = sampleName + "_" + clusterParams.namePrefix;
			clusterParams.basePath = outputPrefix + std::to_string(i);
			clusterParams.hifiReadPath = outputPrefix + std::to_string(i) + "/hifi_reads.fa";
			if (doUL)
			{
				clusterParams.ontReadPath = outputPrefix + std::to_string(i) + "/ont_reads.fa";
			}
			clusterParams.MBGPath = MBGPath;
			clusterParams.GraphAlignerPath = GraphAlignerPath;
			clusterParams.k = k;
			clusterParams.orientReferencePath = orientReferencePath;
			clusterParams.annotationFasta = annotationFasta;
			clusterParams.annotationGff3 = annotationGff3;
			DoClusterONTAnalysis(clusterParams);
		}
	}
	if (tanglesWithoutReads.size() > 0)
	{
		std::cerr << "WARNING: some tangles did not have any HiFi/duplex reads assigned:";
		for (auto tangle : tanglesWithoutReads) std::cerr << " " << tangle;
		std::cerr << ", something likely went wrong." << std::endl;
	}
}
