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
	if (params.count("c") == 1 && params.count("guess-tangles-using-reference") == 1)
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
	else if (params.count("guess-tangles-using-reference") == 1)
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
	std::vector<std::string> readFileNames;
	for (size_t i = 0; i < numTangles; i++)
	{
		std::filesystem::create_directories(outputPrefix + std::to_string(i));
		readFileNames.push_back(outputPrefix + std::to_string(i) + "/hifi_reads.fa");
		writeNodes(outputPrefix + std::to_string(i) + "/nodes.txt", tangleNodes[i]);
	}
	std::cerr << "extracting HiFi/duplex reads per tangle" << std::endl;
	bool allReadsFound = splitReads(hifiReadPaths, reads, readFileNames);
	if (!allReadsFound)
	{
		std::cerr << "WARNING: some HiFi reads were not found in the input read files. Double check that all HiFi reads given to hifiasm are inputed with -i" << std::endl;
	}
	std::vector<size_t> tanglesWithoutReads;
	std::string selectedONTPath = ulTmpFolder + "/ont_reads.fa";
	for (size_t i = 0; i < numTangles; i++)
	{
		if (reads[i].size() == 0)
		{
			std::cerr << "WARNING: tangle " << i << " has no HiFi/duplex reads, skipping" << std::endl;
			tanglesWithoutReads.push_back(i);
			continue;
		}
		ClusterParams clusterParams;
		commonParams.addToClusterOptions(clusterParams);
		clusterParams.namePrefix = clusterParams.namePrefix + (clusterParams.namePrefix.size() > 0 ? "_" : "") + "tangle" + std::to_string(i);
		clusterParams.basePath = outputPrefix + std::to_string(i);
		clusterParams.hifiReadPath = outputPrefix + std::to_string(i) + "/hifi_reads.fa";
		if (doUL)
		{
			clusterParams.ontReadPath = selectedONTPath;
		}
		std::cerr << "running tangle " << i << " in folder " << outputPrefix + std::to_string(i) << std::endl;
		HandleCluster(clusterParams);
	}
	if (doUL)
	{
		std::filesystem::create_directories(ulTmpFolder);
		std::cerr << "getting kmers from tangles" << std::endl;
		getKmers(outputPrefix, numTangles, ulTmpFolder + "/rdna_kmers.fa");
		std::cerr << "extracting ultralong ONT reads" << std::endl;
		std::ofstream readsfile { selectedONTPath };
		size_t consensusLength = medianConsensusLength(outputPrefix, numTangles);
		std::cerr << "median consensus length " << consensusLength << ", using " << consensusLength/2 << " as minimum ONT match length" << std::endl;
		iterateMatchingReads(ulTmpFolder + "/rdna_kmers.fa", ontReadPaths, 21, consensusLength/2, [&readsfile](const FastQ& seq)
		{
			readsfile << ">" << seq.seq_id << std::endl;
			readsfile << seq.sequence << std::endl;
		});
		std::cerr << "merging allele graphs" << std::endl;
		mergeGraphs(outputPrefix, numTangles, ulTmpFolder + "/merged-allele-graph.gfa");
		std::cerr << "aligning ONT reads" << std::endl;
		AlignONTReads(ulTmpFolder, commonParams.GraphAlignerPath(), ulTmpFolder + "/ont_reads.fa", ulTmpFolder + "/merged-allele-graph.gfa", ulTmpFolder + "/ont-alns.gaf", commonParams.numThreads());
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
			commonParams.addToClusterOptions(clusterParams);
			clusterParams.namePrefix = clusterParams.namePrefix + (clusterParams.namePrefix.size() > 0 ? "_" : "") + "tangle" + std::to_string(i);
			clusterParams.basePath = outputPrefix + std::to_string(i);
			clusterParams.hifiReadPath = outputPrefix + std::to_string(i) + "/hifi_reads.fa";
			if (doUL)
			{
				clusterParams.ontReadPath = selectedONTPath;
			}
			DoClusterONTAnalysis(clusterParams);
		}
	}
	std::cerr << "done" << std::endl;
	if (tanglesWithoutReads.size() > 0)
	{
		std::cerr << "WARNING: some tangles did not have any HiFi/duplex reads assigned:";
		for (auto tangle : tanglesWithoutReads) std::cerr << " " << tangle;
		std::cerr << ", something likely went wrong." << std::endl;
	}
}
