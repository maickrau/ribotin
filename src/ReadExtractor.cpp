#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <phmap.h>
#include "ReadExtractor.h"
#include "fastqloader.h"
#include "RibotinUtils.h"
#include "Logger.h"
#include "ClusterMisc.h"
#include "HeavyPath.h"
#include "KmerMatcher.h"

bool splitReads(std::vector<std::string> readFiles, const std::vector<std::vector<std::string>>& readsPerCluster, const std::vector<std::string>& outputFileNames)
{
	assert(readsPerCluster.size() == outputFileNames.size());
	std::unordered_map<std::string, std::unordered_set<size_t>> readBelongsToCluster;
	for (size_t i = 0; i < readsPerCluster.size(); i++)
	{
		for (const auto& read : readsPerCluster[i])
		{
			readBelongsToCluster[nameWithoutTags(read)].insert(i);
		}
	}
	std::vector<std::ofstream> outfiles;
	outfiles.reserve(outputFileNames.size());
	for (size_t i = 0; i < outputFileNames.size(); i++)
	{
		outfiles.emplace_back(outputFileNames[i]);
	}
	phmap::flat_hash_set<std::string> foundReads;
	for (const auto& file : readFiles)
	{
		Logger::Log.log(Logger::LogLevel::Always) << "extracting reads from " << file << std::endl;
		FastQ::streamFastqFromFile(std::string{ file }, false, [&readBelongsToCluster, &outfiles, &foundReads](FastQ& fastq)
		{
			if (readBelongsToCluster.count(nameWithoutTags(fastq.seq_id)) == 0) return;
			foundReads.insert(nameWithoutTags(fastq.seq_id));
			for (auto cluster : readBelongsToCluster.at(nameWithoutTags(fastq.seq_id)))
			{
				assert(cluster < outfiles.size());
				outfiles[cluster] << ">" << fastq.seq_id << std::endl;
				outfiles[cluster] << fastq.sequence << std::endl;
			}
		});
	}
	bool allReadsFound = true;
	for (size_t i = 0; i < readsPerCluster.size(); i++)
	{
		for (const auto& read : readsPerCluster[i])
		{
			if (foundReads.count(nameWithoutTags(read)) == 1) continue;
			Logger::Log.log(Logger::LogLevel::Always) << "WARNING: read " << nameWithoutTags(read) << " not found in input files" << std::endl;
			allReadsFound = false;
		}
	}
	return allReadsFound;
}

void extractMoreHifiReads(std::vector<std::string> wholeGenomeHifiFastas, const std::vector<std::string>& tempReadFileNames, const std::vector<std::string>& realReadFileNames, const std::string MBGPath, const std::string tmpFolder)
{
	const size_t k = 101;
	assert(tempReadFileNames.size() >= 1);
	assert(tempReadFileNames.size() == realReadFileNames.size());
	for (size_t tangle = 0; tangle < tempReadFileNames.size(); tangle++)
	{
		std::string graphFileName = tmpFolder + "/kmergraph" + std::to_string(tangle) + ".gfa";
		std::string mbgCommand = MBGPath + " -i " + tempReadFileNames[tangle] + " -o " + graphFileName + " -k " + std::to_string(k) + " -a 2 -u 5 --error-masking=msat 1> " + tmpFolder + "/mbgstdout" + std::to_string(tangle) + ".txt 2> " + tmpFolder + "/mbgstderr" + std::to_string(tangle) + ".txt";
		Logger::Log.log(Logger::LogLevel::Always) << "MBG command:" << std::endl;
		Logger::Log.log(Logger::LogLevel::Always) << mbgCommand << std::endl;
		int result = system(mbgCommand.c_str());
		if (result != 0)
		{
			Logger::Log.log(Logger::LogLevel::Always) << "MBG did not run successfully" << std::endl;
			std::abort();
		}
		GfaGraph graph;
		graph.loadFromFile(graphFileName);
		std::vector<bool> keptNodes;
		keptNodes.resize(graph.numNodes(), false);
		for (size_t node = 0; node < graph.numNodes(); node++)
		{
			if (getSelfDistance(graph, node) == std::numeric_limits<size_t>::max()) continue;
			keptNodes[node] = true;
		}
		std::ofstream kmerfile { tmpFolder + "/cycle_kmers" + std::to_string(tangle) + ".fa" };
		for (size_t node = 0; node < graph.numNodes(); node++)
		{
			if (!keptNodes[node]) continue;
			kmerfile << ">tangle" << tangle << "kmers" << node << std::endl;
			kmerfile << graph.nodeSeqs[node] << std::endl;
		}
	}
	phmap::flat_hash_map<std::string, size_t> alreadyAssignedReads;
	for (size_t tangle = 0; tangle < tempReadFileNames.size(); tangle++)
	{
		FastQ::streamFastqFromFile(tempReadFileNames[tangle], false, [&alreadyAssignedReads, tangle](const FastQ& read)
		{
			assert(alreadyAssignedReads.count(read.seq_id) == 0);
			alreadyAssignedReads[read.seq_id] = tangle;
		});
	}
	KmerMatcher matcher { k };
	std::vector<KmerMatcher> perTangleMatchers;
	for (size_t tangle = 0; tangle < tempReadFileNames.size(); tangle++)
	{
		perTangleMatchers.emplace_back(k);
		FastQ::streamFastqFromFile(tmpFolder + "/cycle_kmers" + std::to_string(tangle) + ".fa", false, [&matcher, &perTangleMatchers, tangle](FastQ& fastq)
		{
			matcher.addReferenceKmers(fastq.sequence);
			perTangleMatchers[tangle].addReferenceKmers(fastq.sequence);
		});
	}
	for (size_t tangle = 0; tangle < tempReadFileNames.size(); tangle++)
	{
		FastQ::streamFastqFromFile(tmpFolder + "/cycle_kmers" + std::to_string(tangle) + ".fa", false, [&matcher, &perTangleMatchers, tangle](FastQ& fastq)
		{
			for (size_t otherTangle = 0; otherTangle < perTangleMatchers.size(); otherTangle++)
			{
				if (otherTangle == tangle) continue;
				perTangleMatchers[otherTangle].removeReferenceKmers(fastq.sequence);
			}
		});
	}
	std::vector<std::ofstream> resultFiles;
	for (size_t tangle = 0; tangle < realReadFileNames.size(); tangle++)
	{
		resultFiles.emplace_back(realReadFileNames[tangle]);
	}
	std::vector<size_t> addCounts;
	addCounts.resize(tempReadFileNames.size(), 0);
	for (const auto& file : wholeGenomeHifiFastas)
	{
		FastQ::streamFastqFromFile(file, false, [&matcher, &alreadyAssignedReads, &resultFiles, &perTangleMatchers, &addCounts](FastQ& fastq)
		{
			if (alreadyAssignedReads.count(fastq.seq_id) == 1)
			{
				resultFiles[alreadyAssignedReads.at(fastq.seq_id)] << ">" << fastq.seq_id << std::endl;
				resultFiles[alreadyAssignedReads.at(fastq.seq_id)] << fastq.sequence << std::endl;
				return;
			}
			size_t matchLength = matcher.getMatchLength(fastq.sequence);
			if (matchLength < 2000) return;
			size_t assignedTangle = std::numeric_limits<size_t>::max();
			size_t bestMatchCount = 0;
			for (size_t tangle = 0; tangle < perTangleMatchers.size(); tangle++)
			{
				size_t matchCount = perTangleMatchers[tangle].getMatchLength(fastq.sequence);
				if (matchCount == bestMatchCount)
				{
					assignedTangle = std::numeric_limits<size_t>::max();
				}
				if (matchCount > bestMatchCount)
				{
					bestMatchCount = matchCount;
					assignedTangle = tangle;
				}
			}
			if (assignedTangle == std::numeric_limits<size_t>::max()) return;
			resultFiles[assignedTangle] << ">" << fastq.seq_id << std::endl;
			resultFiles[assignedTangle] << fastq.sequence << std::endl;
			addCounts[assignedTangle] += 1;
		});
	}
	Logger::Log.log(Logger::LogLevel::Always) << "new reads added per tangle:";
	for (size_t i = 0; i < addCounts.size(); i++)
	{
		Logger::Log.log(Logger::LogLevel::Always) << " " << addCounts[i];
	}
	Logger::Log.log(Logger::LogLevel::Always) << std::endl;
}
