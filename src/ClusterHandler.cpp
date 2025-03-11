#include <limits>
#include <queue>
#include <map>
#include <cstdlib>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <tuple>
#include <cassert>
#include <set>
#include <vector>
#include <algorithm>
#include <sstream>
#include <phmap.h>
#include <thread>
#include <chrono>
#include <filesystem>
#include "ReadExtractor.h"
#include "fastqloader.h"
#include "ClusterHandler.h"
#include "RibotinUtils.h"
#include "WfaHelper.h"
#include "TwobitString.h"
#include "edlib.h"
#include "ClusterMisc.h"
#include "HeavyPath.h"
#include "SequenceAligner.h"
#include "LoopPhasing.h"
#include "ConsensusHelper.h"
#include "FileHelper.h"

void runMBG(std::string basePath, std::string readPath, std::string MBGPath, const size_t k, const size_t maxResolveLength, const size_t numThreads)
{
	std::string mbgCommand;
	mbgCommand = MBGPath + " -o " + basePath + "/graph.gfa -t " + std::to_string(numThreads) + " -i " + readPath + " -k " + std::to_string(k) + " -w " + std::to_string(k-30) + " -a 2 -u 3 -r " + std::to_string(maxResolveLength) + " -R 4000 --error-masking=msat --output-sequence-paths " + basePath + "/paths.gaf --only-local-resolve --resolve-palindromes-global 1> " + basePath + "/tmp/mbg_stdout.txt 2> " + basePath + "/tmp/mbg_stderr.txt";
	std::cerr << "MBG command:" << std::endl;
	std::cerr << mbgCommand << std::endl;
	int result = system(mbgCommand.c_str());
	if (result != 0)
	{
		std::cerr << "MBG did not run successfully" << std::endl;
		std::abort();
	}
}

void liftoverAnnotationsToMorphs(const std::string& basePath, const std::vector<MorphConsensus>& morphConsensuses, const std::string& annotationFasta, const std::string& annotationGff3, const std::string& tmpPath, const std::string& liftoffPath)
{
	{
		std::ofstream typefile { tmpPath + "/liftoff_types.txt" };
		typefile << "rRNA" << std::endl;
		typefile << "misc_RNA" << std::endl;
		typefile << "repeat_region" << std::endl;
		typefile << "gene" << std::endl;
		typefile << "transcript" << std::endl;
		typefile << "exon" << std::endl;
		typefile << "pseudogene" << std::endl;
		typefile << "tandem_repeat" << std::endl;
	}
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		std::string tmpfilepath { tmpPath + "/tmp_seq.fa" };
		{
			std::ofstream tmpfile { tmpfilepath };
			tmpfile << ">" << morphConsensuses[i].name << std::endl;
			tmpfile << morphConsensuses[i].sequence << std::endl;
		}
		std::string command = liftoffPath + " -f " + tmpPath + "/liftoff_types.txt -g " + annotationGff3 + " -o " + tmpPath + "/tmp-morph-annotations-part" + std::to_string(i) + ".gff3 -u "+ tmpPath + "/morph-unmapped_features" + std::to_string(i) + ".txt -dir " + tmpPath + "/liftoff_intermediate_files/ " + tmpfilepath + " " + annotationFasta + " 1> " + tmpPath + "/liftoff_morphs_stdout" + std::to_string(i) + ".txt 2> " + tmpPath + "/liftoff_morphs_stderr" + std::to_string(i) + ".txt";
		std::cerr << "running liftoff with command:" << std::endl;
		std::cerr << command << std::endl;
		int result = system(command.c_str());
		if (result != 0)
		{
			std::cerr << "liftoff did not run successfully" << std::endl;
			std::abort();
		}
		command = "rm -r " + tmpPath + "/liftoff_intermediate_files/ " + tmpfilepath + ".fai " + tmpfilepath + ".mmi";
		std::cerr << command << std::endl;
		result = system(command.c_str());
	}
	std::cerr << "combining liftoff results" << std::endl;
	std::string outputFile = basePath + "/morph-annotations.gff3";
	std::string command = "echo \"##gff-version 3\" > " + outputFile + "\ncat " + tmpPath + "/tmp-morph-annotations-part*.gff3 | grep -v '#' >> " + outputFile;
	std::cerr << command << std::endl;
	int result = system(command.c_str());
	if (result != 0)
	{
		std::cerr << "failed to combine liftoff results" << std::endl;
		std::abort();
	}
}

void liftoverAnnotationsToConsensus(const std::string& basePath, const std::string& consensusPath, const std::string& annotationFasta, const std::string& annotationGff3, const std::string& tmpPath, const std::string& liftoffPath)
{
	{
		std::ofstream typefile { tmpPath + "/liftoff_types.txt" };
		typefile << "rRNA" << std::endl;
		typefile << "misc_RNA" << std::endl;
		typefile << "repeat_region" << std::endl;
		typefile << "gene" << std::endl;
		typefile << "transcript" << std::endl;
		typefile << "exon" << std::endl;
		typefile << "pseudogene" << std::endl;
		typefile << "tandem_repeat" << std::endl;
	}
	std::string command = liftoffPath + " -f " + tmpPath + "/liftoff_types.txt -g " + annotationGff3 + " -o " + basePath + "/consensus-annotation.gff3 -u " + tmpPath + "/consensus-unmapped_features.txt -dir " + tmpPath + "/liftoff_intermediate_files/ " + consensusPath + " " + annotationFasta + " 1> " + tmpPath + "/liftoff_consensus_stdout.txt 2> " + tmpPath + "/liftoff_consensus_stderr.txt";
	std::cerr << "running liftoff with command:" << std::endl;
	std::cerr << command << std::endl;
	int result = system(command.c_str());
	if (result != 0)
	{
		std::cerr << "liftoff did not run successfully" << std::endl;
		std::abort();
	}
}

void AlignONTReads(std::string basePath, std::string graphAlignerPath, std::string ontReadPath, std::string graphPath, std::string outputAlnPath, size_t numThreads)
{
	std::string graphalignerCommand;
	graphalignerCommand = graphAlignerPath + " -g " + graphPath + " -f " + ontReadPath + " -a " + outputAlnPath + " -t " + std::to_string(numThreads) + " --seeds-mxm-length 30 --seeds-mem-count 10000 --bandwidth 15 --multimap-score-fraction 0.99 --precise-clipping 0.85 --min-alignment-score 5000 --discard-cigar --clip-ambiguous-ends 100 --overlap-incompatible-cutoff 0.15 --mem-index-no-wavelet-tree --max-trace-count 5 1> " + basePath + "/tmp/graphaligner_stdout.txt 2> " + basePath + "/tmp/graphaligner_stderr.txt";
	std::cerr << "GraphAligner command:" << std::endl;
	std::cerr << graphalignerCommand << std::endl;
	int result = system(graphalignerCommand.c_str());
	if (result != 0)
	{
		std::cerr << "GraphAligner did not run successfully" << std::endl;
		std::abort();
	}
}

size_t getDPRowBacktracePos(const std::string& compareQuerySeq, const std::string& compareRefSeq, const size_t columnNum)
{
	std::vector<std::vector<size_t>> DPmatrix;
	DPmatrix.resize(compareQuerySeq.size()+1);
	for (size_t i = 0; i < compareQuerySeq.size()+1; i++)
	{
		DPmatrix[i].resize(compareRefSeq.size()+1, 0);
	}
	for (size_t i = 0; i < compareRefSeq.size()+1; i++)
	{
		DPmatrix[0][i] = i;
	}
	for (size_t i = 1; i < compareQuerySeq.size()+1; i++)
	{
		DPmatrix[i][0] = i;
		for (size_t j = 1; j < compareRefSeq.size()+1; j++)
		{
			bool match = (compareQuerySeq[i-1] == compareRefSeq[j-1]);
			DPmatrix[i][j] = std::min(DPmatrix[i][j-1]+1, DPmatrix[i-1][j]+1);
			DPmatrix[i][j] = std::min(DPmatrix[i][j], DPmatrix[i-1][j-1] + (match ? 0 : 1));
		}
	}
	size_t i = compareQuerySeq.size();
	size_t j = compareRefSeq.size();
	while (j > columnNum)
	{
		assert(j > 0);
		if (i == 0)
		{
			j -= 1;
			continue;
		}
		assert(i > 0);
		bool match = (compareQuerySeq[i-1] == compareRefSeq[j-1]);
		if (DPmatrix[i-1][j-1] + (match ? 0 : 1) == DPmatrix[i][j])
		{
			i -= 1;
			j -= 1;
		}
		else if (DPmatrix[i-1][j] + 1 == DPmatrix[i][j])
		{
			i -= 1;
		}
		else if (DPmatrix[i][j-1] + 1 == DPmatrix[i][j])
		{
			j -= 1;
		}
		else
		{
			assert(false);
		}
	}
	return i;
}

size_t getExactBreakPos(const std::string& nodeseq, const std::string& consensusSeq, const size_t approxPosition)
{
	const size_t flankSize = 10;
	assert(nodeseq.size() >= flankSize);
	assert(consensusSeq.size() >= flankSize);
	std::string compareQuerySeq;
	std::string compareRefSeq;
	size_t midPos = 0;
	if (approxPosition < flankSize)
	{
		compareQuerySeq += nodeseq.substr(0, approxPosition);
		compareRefSeq += consensusSeq.substr(consensusSeq.size() - approxPosition);
		midPos = approxPosition;
	}
	else
	{
		compareQuerySeq += nodeseq.substr(approxPosition-flankSize, flankSize);
		compareRefSeq += consensusSeq.substr(consensusSeq.size() - flankSize);
		midPos = flankSize;
	}
	assert(compareRefSeq.size() == compareQuerySeq.size());
	if (approxPosition+flankSize > nodeseq.size())
	{
		compareQuerySeq += nodeseq.substr(approxPosition, flankSize);
		compareRefSeq += consensusSeq.substr(0, nodeseq.size() - approxPosition);
	}
	else
	{
		compareQuerySeq += nodeseq.substr(approxPosition, flankSize);
		compareRefSeq += consensusSeq.substr(0, flankSize);
	}
	assert(compareRefSeq.size() == compareQuerySeq.size());
	size_t zeroPos = approxPosition;
	if (approxPosition >= flankSize)
	{
		zeroPos = approxPosition - flankSize;
	}
	else
	{
		zeroPos = 0;
	}
	size_t exactMatchPos = zeroPos + getDPRowBacktracePos(compareQuerySeq, compareRefSeq, midPos);
	return exactMatchPos;
}

std::tuple<std::unordered_set<size_t>, std::unordered_set<size_t>, std::unordered_map<Node, size_t>, std::unordered_map<Node, size_t>> getBorderNodes(const Path& heavyPath, const GfaGraph& graph)
{
	const size_t k = 101;
	const size_t borderSize = 200;
	std::string consensusSequence = heavyPath.getSequence(graph.nodeSeqs);
	std::unordered_map<std::string, size_t> startKmers;
	std::unordered_map<std::string, size_t> endKmers;
	for (size_t i = 0; i < borderSize; i++)
	{
		std::string subs = consensusSequence.substr(i, k);
		if (endKmers.count(subs) == 1 || startKmers.count(subs) == 1)
		{
			endKmers[subs] = std::numeric_limits<size_t>::max();
			startKmers[subs] = std::numeric_limits<size_t>::max();
		}
		else
		{
			startKmers[subs] = i;
		}
		subs = consensusSequence.substr(consensusSequence.size()-k-i, k);
		if (endKmers.count(subs) == 1 || startKmers.count(subs) == 1)
		{
			endKmers[subs] = std::numeric_limits<size_t>::max();
			startKmers[subs] = std::numeric_limits<size_t>::max();
		}
		else
		{
			endKmers[subs] = i;
		}
	}
	for (size_t i = borderSize; i < consensusSequence.size() - borderSize; i++)
	{
		std::string subs = consensusSequence.substr(i, k);
		if (endKmers.count(subs) == 1 || startKmers.count(subs) == 1)
		{
			endKmers[subs] = std::numeric_limits<size_t>::max();
			startKmers[subs] = std::numeric_limits<size_t>::max();
		}
	}
	std::unordered_set<size_t> borderNodes;
	std::unordered_set<size_t> anchorNodes;
	std::unordered_map<Node, size_t> pathStartClip;
	std::unordered_map<Node, size_t> pathEndClip;
	for (size_t nodeid = 0; nodeid < graph.nodeSeqs.size(); nodeid++)
	{
		if (graph.edges.count(Node { nodeid, true }) == 0)
		{
			pathStartClip[Node { nodeid, false }] = 0;
			pathEndClip[Node { nodeid, true }] = graph.nodeSeqs[nodeid].size();
			anchorNodes.emplace(nodeid);
			borderNodes.emplace(nodeid);
		}
		if (graph.edges.count(Node { nodeid, false }) == 0)
		{
			pathStartClip[Node { nodeid, true }] = 0;
			pathEndClip[Node { nodeid, false }] = graph.nodeSeqs[nodeid].size();
			anchorNodes.emplace(nodeid);
			borderNodes.emplace(nodeid);
		}
	}
	for (size_t nodeid = 0; nodeid < graph.nodeSeqs.size(); nodeid++)
	{
		std::vector<size_t> fwBreaks;
		for (size_t i = 0; i < graph.nodeSeqs[nodeid].size()-k; i++)
		{
			std::string subs = graph.nodeSeqs[nodeid].substr(i, k);
			if (startKmers.count(subs) == 1)
			{
				size_t pos = startKmers.at(subs);
				if (pos != std::numeric_limits<size_t>::max() && pos <= i)
				{
					fwBreaks.push_back(i - pos);
				}
			}
			if (endKmers.count(subs) == 1)
			{
				size_t pos = endKmers.at(subs);
				if (pos != std::numeric_limits<size_t>::max() && i+pos+k < graph.nodeSeqs[nodeid].size())
				{
					fwBreaks.push_back(i+pos+k);
				}
			}
		}
		std::vector<size_t> bwBreaks;
		std::string reverseSeq = revcomp(graph.nodeSeqs[nodeid]);
		for (size_t i = 0; i < reverseSeq.size()-k; i++)
		{
			std::string subs = reverseSeq.substr(i, k);
			if (startKmers.count(subs) == 1)
			{
				size_t pos = startKmers.at(subs);
				if (pos != std::numeric_limits<size_t>::max() && pos <= i)
				{
					bwBreaks.push_back(i - pos);
				}
			}
			if (endKmers.count(subs) == 1)
			{
				size_t pos = endKmers.at(subs);
				if (pos != std::numeric_limits<size_t>::max() && i+pos+k < reverseSeq.size())
				{
					bwBreaks.push_back(i+pos+k);
				}
			}
		}
		if (fwBreaks.size() >= 3)
		{
			size_t breakPos = getExactBreakPos(graph.nodeSeqs[nodeid], consensusSequence, fwBreaks[fwBreaks.size()/2]);
			// size_t breakPos = fwBreaks[fwBreaks.size()/2];
			pathStartClip[Node { nodeid, true }] = breakPos;
			pathEndClip[Node { nodeid, true }] = graph.nodeSeqs[nodeid].size() - breakPos;
			breakPos = graph.nodeSeqs[nodeid].size() - 1 - breakPos;
			pathStartClip[Node { nodeid, false }] = breakPos;
			pathEndClip[Node { nodeid, false }] = graph.nodeSeqs[nodeid].size() - breakPos;
			borderNodes.insert(nodeid);
		}
		if (bwBreaks.size() >= 3)
		{
			size_t breakPos = getExactBreakPos(revcomp(graph.nodeSeqs[nodeid]), consensusSequence, bwBreaks[fwBreaks.size()/2]);
			// size_t breakPos = bwBreaks[bwBreaks.size()/2];
			pathStartClip[Node { nodeid, false }] = breakPos;
			pathEndClip[Node { nodeid, false }] = graph.nodeSeqs[nodeid].size() - breakPos;
			breakPos = graph.nodeSeqs[nodeid].size() - 1 - breakPos;
			pathStartClip[Node { nodeid, true }] = breakPos;
			pathEndClip[Node { nodeid, true }] = graph.nodeSeqs[nodeid].size() - breakPos;
			borderNodes.insert(nodeid);
		}
	}
	return std::make_tuple(borderNodes, anchorNodes, pathStartClip, pathEndClip);
}

std::vector<OntLoop> extractLoopSequences(const std::vector<ReadPath>& correctedPaths, const Path& heavyPath, const size_t minLength, const GfaGraph& graph, const std::unordered_set<size_t>& borderNodes, const std::unordered_set<size_t>& anchorNodes, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
{
	std::vector<OntLoop> result;
	for (const auto& read : correctedPaths)
	{
		size_t lastBreak = std::numeric_limits<size_t>::max();
		std::vector<size_t> approxReadStartPoses;
		std::vector<size_t> approxReadEndPoses;
		size_t pathLength = getPathLength(read.path, graph.nodeSeqs, graph.edges);
		assert(pathLength > read.pathEndClip);
		size_t currentPos = 0;
		for (size_t i = 0; i < read.path.size(); i++)
		{
			approxReadStartPoses.push_back(currentPos);
			approxReadEndPoses.push_back(currentPos + graph.nodeSeqs.at(read.path[i].id()).size());
			currentPos += graph.nodeSeqs.at(read.path[i].id()).size();
			if (i+1 < read.path.size()) currentPos -= getOverlap(read.path[i], read.path[i+1], graph.edges);
		}
		size_t minPos = pathLength;
		size_t maxPos = 0;
		for (size_t i = 0; i < read.path.size(); i++)
		{
			if (approxReadStartPoses[i] < read.pathStartClip) approxReadStartPoses[i] = read.pathStartClip;
			if (approxReadEndPoses[i] < read.pathStartClip) approxReadEndPoses[i] = read.pathStartClip;
			if (approxReadEndPoses[i] > pathLength - read.pathEndClip) approxReadEndPoses[i] = pathLength - read.pathEndClip;
			if (approxReadStartPoses[i] > pathLength - read.pathEndClip) approxReadStartPoses[i] = pathLength - read.pathEndClip;
			assert(approxReadStartPoses[i] <= approxReadEndPoses[i]);
			minPos = std::min(minPos, approxReadStartPoses[i]);
			maxPos = std::max(maxPos, approxReadEndPoses[i]);
		}
		double readScaleMultiplier = (double)(read.readEnd - read.readStart) / (double)(maxPos - minPos);
		for (size_t i = 0; i < read.path.size(); i++)
		{
			approxReadStartPoses[i] -= minPos;
			approxReadStartPoses[i] *= readScaleMultiplier;
			approxReadStartPoses[i] += read.readStart;
			approxReadEndPoses[i] -= minPos;
			approxReadEndPoses[i] *= readScaleMultiplier;
			approxReadEndPoses[i] += read.readStart;
		}
		for (size_t i = 0; i < read.path.size(); i++)
		{
			if (borderNodes.count(read.path[i].id()) == 0) continue;
			if (lastBreak == std::numeric_limits<size_t>::max())
			{
				lastBreak = i;
				continue;
			}
			std::vector<Node> looppath { read.path.begin()+lastBreak, read.path.begin()+i+1 };
			size_t len = getPathLength(looppath, graph.nodeSeqs, graph.edges);
			size_t clip = pathStartClip.at(read.path[lastBreak]) + pathEndClip.at(read.path[i]);
			bool anchored = anchorNodes.count(read.path[lastBreak].id()) == 1 || anchorNodes.count(read.path[i].id()) == 1;
			if (len > clip && (len - clip >= minLength || anchored))
			{
				result.emplace_back();
				result.back().originalReadLength = read.readLength;
				result.back().path = looppath;
				result.back().readName = nameWithoutTags(read.readName);
				result.back().approxStart = approxReadStartPoses[lastBreak] + pathStartClip.at(read.path[lastBreak]);
				result.back().approxEnd = approxReadEndPoses[i] - pathEndClip.at(read.path[i]);
				if (read.reverse)
				{
					result.back().approxStart = read.readEnd - (result.back().approxStart - read.readStart);
					result.back().approxEnd = read.readEnd - (result.back().approxEnd - read.readStart);
				}
			}
			lastBreak = i;
		}
	}
	return result;
}

size_t countNeedsAligning(const std::vector<size_t>& loopLengths, size_t maxEdits)
{
	size_t result = 0;
	size_t j = 0;
	for (size_t i = 1; i < loopLengths.size(); i++)
	{
		while (loopLengths[j]+maxEdits < loopLengths[i]) j += 1;
		assert(j <= i);
		result += i-j;
	}
	return result;
}

std::vector<std::vector<OntLoop>> roughClusterLoopSequences(const std::vector<OntLoop>& loops, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::unordered_set<size_t>& coreNodes, const size_t maxEdits)
{
	std::vector<size_t> parent;
	parent.resize(loops.size());
	for (size_t i = 0; i < loops.size(); i++)
	{
		parent[i] = i;
	}
	std::vector<size_t> loopLengths;
	loopLengths.reserve(loops.size());
	for (size_t i = 0; i < loops.size(); i++)
	{
		size_t len = getPathLength(loops[i].path, graph.nodeSeqs, graph.edges);
		len -= pathStartClip.at(loops[i].path[0]);
		len -= pathEndClip.at(loops[i].path.back());
		loopLengths.emplace_back(len);
	}
	std::vector<phmap::flat_hash_map<Node, size_t>> nodeCountIndex;
	std::vector<phmap::flat_hash_map<Node, size_t>> nodePosIndex;
	nodeCountIndex.resize(loops.size());
	nodePosIndex.resize(loops.size());
	for (size_t i = 0; i < loops.size(); i++)
	{
		for (size_t j = 0; j < loops[i].path.size(); j++)
		{
			nodeCountIndex[i][loops[i].path[j]] += 1;
			nodePosIndex[i][loops[i].path[j]] = j;
		}
	}
	std::unordered_map<std::pair<std::vector<Node>, std::vector<Node>>, size_t> memoizedEditDistances;
	size_t sumAligned = 0;
	size_t needsAligning = countNeedsAligning(loopLengths, maxEdits);
	for (size_t i = 0; i < loops.size(); i++)
	{
		while (parent[i] != parent[parent[i]]) parent[i] = parent[parent[i]];
		for (size_t j = i-1; j < i; j--)
		{
			if (sumAligned % 1000000 == 0) std::cerr << "aligning morph path pair " << sumAligned << " / " << needsAligning << std::endl;
			sumAligned += 1;
			assert(loopLengths[j] <= loopLengths[i]);
			if (loopLengths[j]+maxEdits < loopLengths[i]) break;
			while (parent[j] != parent[parent[j]]) parent[j] = parent[parent[j]];
			if (parent[i] == parent[j]) continue;
			size_t edits = getEditDistance(loops[i].path, i, loops[j].path, j, graph, pathStartClip, pathEndClip, maxEdits, coreNodes, nodeCountIndex, nodePosIndex, memoizedEditDistances);
			assert(edits == getEditDistance(loops[j].path, j, loops[i].path, i, graph, pathStartClip, pathEndClip, maxEdits, coreNodes, nodeCountIndex, nodePosIndex, memoizedEditDistances));
			if (edits > maxEdits) continue;
			parent[parent[j]] = parent[i];
		}
	}
	std::vector<size_t> parentToCluster;
	parentToCluster.resize(parent.size(), std::numeric_limits<size_t>::max());
	std::vector<std::vector<OntLoop>> result;
	for (size_t i = 0; i < parent.size(); i++)
	{
		while (parent[i] != parent[parent[i]]) parent[i] = parent[parent[i]];
		if (parent[i] != i) continue;
		parentToCluster[i] = result.size();
		result.emplace_back();
	}
	for (size_t i = 0; i < parent.size(); i++)
	{
		result[parentToCluster[parent[i]]].emplace_back(loops[i]);
	}
	return result;
}

std::vector<std::vector<OntLoop>> clusterByDbscan(const std::vector<OntLoop>& cluster, const size_t epsilon, const size_t minPoints, const std::vector<std::vector<uint16_t>>& editDistanceMatrix)
{
	// O(n^2) is fast enough for the cluster sizes here
	std::vector<bool> corePoint;
	corePoint.resize(cluster.size(), false);
	for (size_t i = 0; i < cluster.size(); i++)
	{
		size_t nearPoints = 0;
		for (size_t j = 0; j < cluster.size(); j++)
		{
			if (i == j) continue;
			size_t canoni = std::max(i, j);
			size_t canonj = std::min(i, j);
			if (editDistanceMatrix[canoni][canonj] > epsilon) continue;
			nearPoints += 1;
			if (nearPoints == minPoints)
			{
				corePoint[i] = true;
				break;
			}
		}
	}
	std::vector<size_t> parent;
	parent.resize(cluster.size(), 0);
	for (size_t i = 0; i < cluster.size(); i++)
	{
		parent[i] = i;
	}
	for (size_t i = 0; i < cluster.size(); i++)
	{
		if (!corePoint[i]) continue;
		for (size_t j = 0; j < cluster.size(); j++)
		{
			if (i == j) continue;
			if (!corePoint[j]) continue;
			size_t canoni = std::max(i, j);
			size_t canonj = std::min(i, j);
			if (editDistanceMatrix[canoni][canonj] > epsilon) continue;
			while (parent[parent[i]] != parent[i]) parent[i] = parent[parent[i]];
			while (parent[parent[j]] != parent[j]) parent[j] = parent[parent[j]];
			parent[parent[j]] = parent[i];
		}
	}
	for (size_t j = 0; j < parent.size(); j++)
	{
		while (parent[j] != parent[parent[j]]) parent[j] = parent[parent[j]];
	}
	std::vector<size_t> uniqueClusterForEdgePoint;
	uniqueClusterForEdgePoint.resize(cluster.size(), std::numeric_limits<size_t>::max()-1);
	for (size_t i = 0; i < cluster.size(); i++)
	{
		if (corePoint[i]) continue;
		assert(parent[i] == i);
		for (size_t j = 0; j < cluster.size(); j++)
		{
			if (i == j) continue;
			if (!corePoint[j]) continue;
			size_t canoni = std::max(i, j);
			size_t canonj = std::min(i, j);
			if (editDistanceMatrix[canoni][canonj] > epsilon) continue;
			if (uniqueClusterForEdgePoint[i] == std::numeric_limits<size_t>::max()-1)
			{
				assert(parent[j] == parent[parent[j]]);
				uniqueClusterForEdgePoint[i] = parent[j];
				parent[i] = parent[j];
			}
			else if (uniqueClusterForEdgePoint[i] != parent[j])
			{
				uniqueClusterForEdgePoint[i] = std::numeric_limits<size_t>::max();
			}
		}
	}
	std::vector<size_t> clusterMapping;
	std::vector<std::vector<OntLoop>> result;
	clusterMapping.resize(cluster.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < cluster.size(); i++)
	{
		if (!corePoint[i] && (uniqueClusterForEdgePoint[i] == std::numeric_limits<size_t>::max() || uniqueClusterForEdgePoint[i] == std::numeric_limits<size_t>::max()-1)) continue;
		if (clusterMapping[parent[i]] == std::numeric_limits<size_t>::max())
		{
			clusterMapping[parent[i]] = result.size();
			result.emplace_back();
		}
		result[clusterMapping[parent[i]]].emplace_back(cluster[i]);
	}
	return result;
}

std::vector<std::vector<OntLoop>> densityClusterLoops(const std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::unordered_set<size_t>& coreNodes, const size_t maxEdits, const size_t minPoints, const size_t minEdits)
{
	std::vector<size_t> editHistogram;
	editHistogram.resize(maxEdits, 0);
	std::vector<std::vector<std::vector<uint16_t>>> editDistanceMatrices;
	editDistanceMatrices.resize(clusters.size());
	for (size_t clusteri = 0; clusteri < clusters.size(); clusteri++)
	{
		editDistanceMatrices[clusteri].resize(clusters[clusteri].size());
		std::vector<phmap::flat_hash_map<Node, size_t>> nodeCountIndex;
		std::vector<phmap::flat_hash_map<Node, size_t>> nodePosIndex;
		nodeCountIndex.resize(clusters[clusteri].size());
		nodePosIndex.resize(clusters[clusteri].size());
		for (size_t i = 0; i < clusters[clusteri].size(); i++)
		{
			for (size_t j = 0; j < clusters[clusteri][i].path.size(); j++)
			{
				nodeCountIndex[i][clusters[clusteri][i].path[j]] += 1;
				nodePosIndex[i][clusters[clusteri][i].path[j]] = j;
			}
		}
		std::unordered_map<std::pair<std::vector<Node>, std::vector<Node>>, size_t> memoizedEditDistances;
		for (size_t i = 0; i < clusters[clusteri].size(); i++)
		{
			editDistanceMatrices[clusteri][i].resize(i);
			for (size_t j = 0; j < i; j++)
			{
				size_t edits = getEditDistance(clusters[clusteri][i].path, i, clusters[clusteri][j].path, j, graph, pathStartClip, pathEndClip, maxEdits, coreNodes, nodeCountIndex, nodePosIndex, memoizedEditDistances);
				assert(edits == getEditDistance(clusters[clusteri][j].path, j, clusters[clusteri][i].path, i, graph, pathStartClip, pathEndClip, maxEdits, coreNodes, nodeCountIndex, nodePosIndex, memoizedEditDistances));
				if (edits >= maxEdits) edits = maxEdits-1;
				editDistanceMatrices[clusteri][i][j] = edits;
				editHistogram[edits] += 1;
			}
		}
	}
	size_t histogramPeak = 0;
	for (size_t i = 0; i < editHistogram.size()-1; i++)
	{
		if (editHistogram[i] >= editHistogram[histogramPeak]) histogramPeak = i;
	}
	std::cerr << "edit distance peak at " << histogramPeak << std::endl;
	size_t newEditDistance = histogramPeak;
	if (newEditDistance >= maxEdits)
	{
		return clusters;
	}
	if (newEditDistance < minEdits) newEditDistance = minEdits;
	std::cerr << "recluster with max edit distance " << newEditDistance << ", min points " << minPoints << std::endl;
	std::vector<std::vector<OntLoop>> result;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		auto partialResult = clusterByDbscan(clusters[i], newEditDistance, minPoints, editDistanceMatrices[i]);
		std::cerr << "cluster " << i << " with " << clusters[i].size() << " reads reclustered to " << partialResult.size() << " clusters, sizes:";
		size_t countInClusters = 0;
		for (size_t j = 0; j < partialResult.size(); j++)
		{
			std::cerr << " " << partialResult[j].size();
			countInClusters += partialResult[j].size();
		}
		if (partialResult.size() < 2 || countInClusters < clusters[i].size()*0.9)
		{
			std::cerr << ", don't split" << std::endl;
			result.emplace_back(clusters[i]);
			continue;
		}
		std::cerr << std::endl;
		while (partialResult.size() > 0)
		{
			result.emplace_back();
			std::swap(result.back(), partialResult.back());
			partialResult.pop_back();
		}
	}
	return result;
}

std::vector<MorphConsensus> getMorphConsensuses(const std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::string& namePrefix)
{
	std::vector<MorphConsensus> result;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		result.emplace_back();
		result.back().path = getConsensusPath(clusters[i], graph);
		result.back().preCorrectionSequence = getConsensusSequence(result.back().path, graph, pathStartClip, pathEndClip);
		result.back().sequence = result.back().preCorrectionSequence;
		result.back().coverage = clusters[i].size();
		result.back().ontLoops = clusters[i];
		result.back().name = "morphconsensus" + std::to_string(i) + "_coverage" + std::to_string(result.back().coverage);
		if (namePrefix != "") result.back().name = namePrefix + "_" + result.back().name;
	}
	return result;
}

void polishMorphConsensuses(std::vector<MorphConsensus>& morphConsensuses, const std::vector<std::vector<OntLoop>>& ontLoopSequences, const std::string MBGPath, const std::string tmpPath, const size_t numThreads)
{
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		if (ontLoopSequences[i].size() < 3)
		{
			std::cerr << "skip polishing morph " << i << " with coverage " << ontLoopSequences[i].size() << std::endl;
			continue;
		}
		std::vector<std::string> seqs;
		for (size_t k = 0; k < ontLoopSequences[i].size(); k++)
		{
			seqs.emplace_back(ontLoopSequences[i][k].rawSequence);
		}
		size_t sizeBeforePolish = morphConsensuses[i].sequence.size();
		morphConsensuses[i].sequence = polishConsensus(morphConsensuses[i].sequence, seqs, numThreads);
		std::cerr << "polished from size " << sizeBeforePolish << " to " << morphConsensuses[i].sequence.size() << std::endl;
	}
}

void writeMorphConsensuses(std::string outFile, std::string preconsensusOutputFile, const std::vector<MorphConsensus>& morphConsensuses)
{
	std::ofstream file { outFile };
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		file << ">" << morphConsensuses[i].name << std::endl;
		file << morphConsensuses[i].sequence << std::endl;
	}
	std::ofstream filePreCorrection { preconsensusOutputFile };
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		filePreCorrection << ">" << morphConsensuses[i].name << std::endl;
		filePreCorrection << morphConsensuses[i].preCorrectionSequence << std::endl;
	}
}

void orderLoopsByLength(std::vector<OntLoop>& loops, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
{
	std::vector<std::pair<size_t, size_t>> indexAndLength;
	indexAndLength.reserve(loops.size());
	for (size_t i = 0; i < loops.size(); i++)
	{
		size_t len = getPathLength(loops[i].path, graph.nodeSeqs, graph.edges);
		size_t leftclip = pathStartClip.at(loops[i].path[0]);
		size_t rightclip = pathEndClip.at(loops[i].path.back());
		assert(len > leftclip+rightclip);
		indexAndLength.emplace_back(i, len-leftclip-rightclip);
	}
	std::sort(indexAndLength.begin(), indexAndLength.end(), [](auto left, auto right) { return left.second < right.second; });
	std::vector<OntLoop> result;
	result.resize(loops.size());
	for (size_t i = 0; i < loops.size(); i++)
	{
		std::swap(result[i], loops[indexAndLength[i].first]);
	}
	std::swap(result, loops);
}

void processGraphAndWrite(const Path& heavyPath, const GfaGraph& graph, const std::string outputGraphName)
{
	const size_t minCoverage = 5;
	const size_t minAnchorLength = 10000;
	phmap::flat_hash_set<size_t> coveredNodes;
	std::vector<Node> checkStack;
	phmap::flat_hash_set<Node> reachableNodes;
	for (size_t i = 0; i < heavyPath.nodes.size(); i++)
	{
		coveredNodes.insert(heavyPath.nodes[i].id());
		checkStack.emplace_back(Node { heavyPath.nodes[i].id(), true });
		checkStack.emplace_back(Node { heavyPath.nodes[i].id(), false });
	}
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (graph.nodeCoverages[i] < minCoverage) continue;
		coveredNodes.insert(i);
	}
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (reachableNodes.count(top) == 1) continue;
		reachableNodes.insert(top);
		if (graph.edges.count(top) == 1)
		{
			for (auto edge : graph.edges.at(top))
			{
				auto target = std::get<0>(edge);
				if (coveredNodes.count(target.id()) == 0) continue;
				checkStack.emplace_back(target);
			}
		}
	}
	phmap::flat_hash_set<size_t> keptNodesInCycles;
	for (size_t node : coveredNodes)
	{
		if (reachableNodes.count(Node { node, true }) == 0) continue;
		if (reachableNodes.count(Node { node, false }) == 0) continue;
		keptNodesInCycles.insert(node);
	}
	phmap::flat_hash_map<Node, size_t> shortestDistanceToCyclicComponent;
	std::vector<std::pair<Node, size_t>> checkStack2;
	for (size_t node : keptNodesInCycles)
	{
		checkStack2.emplace_back(Node { node, true }, 0);
		checkStack2.emplace_back(Node { node, false }, 0);
	}
	while (checkStack2.size() >= 1)
	{
		auto top = checkStack2.back();
		checkStack2.pop_back();
		if (shortestDistanceToCyclicComponent.count(top.first) == 1)
		{
			assert(shortestDistanceToCyclicComponent.at(top.first) <= top.second);
			continue;
		}
		shortestDistanceToCyclicComponent[top.first] = top.second;
		bool addedAny = false;
		if (graph.edges.count(top.first) == 1)
		{
			for (auto edge : graph.edges.at(top.first))
			{
				if (keptNodesInCycles.count(std::get<0>(edge).id()) == 1) continue;
				checkStack2.emplace_back(std::get<0>(edge), top.second + graph.nodeSeqs[std::get<0>(edge).id()].size() - std::get<1>(edge));
				assert(checkStack2.back().second > top.second);
				addedAny = true;
			}
		}
		if (addedAny) std::sort(checkStack2.begin(), checkStack2.end(), [](auto left, auto right) { return left.second > right.second; });
	}
	phmap::flat_hash_set<size_t> keptNodesInAnchors;
	std::vector<Node> checkStack3;
	phmap::flat_hash_set<Node> nodesWithEnoughDistance;
	for (auto pair : shortestDistanceToCyclicComponent)
	{
		if (keptNodesInCycles.count(pair.first.id()) == 1) continue;
		if (shortestDistanceToCyclicComponent.count(reverse(pair.first)) == 1) continue;
		if (pair.second < minAnchorLength) continue;
		checkStack3.emplace_back(reverse(pair.first));
		nodesWithEnoughDistance.emplace(reverse(pair.first));
	}
	while (checkStack3.size() >= 1)
	{
		auto top = checkStack3.back();
		checkStack3.pop_back();
		if (keptNodesInAnchors.count(top.id()) == 1) continue;
		keptNodesInAnchors.insert(top.id());
		assert(graph.edges.count(top) == 1);
		for (auto edge : graph.edges.at(top))
		{
			if (keptNodesInCycles.count(std::get<0>(edge).id()) == 1) continue;
			if (shortestDistanceToCyclicComponent.count(reverse(std::get<0>(edge))) == 0) continue;
			checkStack3.emplace_back(std::get<0>(edge));
		}
	}
	for (Node node : nodesWithEnoughDistance)
	{
		assert(keptNodesInAnchors.count(node.id()) == 1);
		assert(graph.edges.count(node) == 1);
		bool allNeighborsHaveEnoughDistance = true;
		for (auto edge : graph.edges.at(node))
		{
			if (nodesWithEnoughDistance.count(std::get<0>(edge)) == 0)
			{
				allNeighborsHaveEnoughDistance = false;
				break;
			}
		}
		if (allNeighborsHaveEnoughDistance)
		{
			keptNodesInAnchors.erase(node.id());
		}
	}
	phmap::flat_hash_set<size_t> finalKeptNodes;
	finalKeptNodes.insert(keptNodesInCycles.begin(), keptNodesInCycles.end());
	finalKeptNodes.insert(keptNodesInAnchors.begin(), keptNodesInAnchors.end());
	std::ofstream file { outputGraphName };
	for (size_t node : finalKeptNodes)
	{
		file << "S\t" << graph.nodeNames[node] << "\t" << graph.nodeSeqs[node] << "\tll:f:" << graph.nodeCoverages[node] << "\tFC:i:" << (graph.nodeCoverages[node] * graph.nodeSeqs[node].size()) << std::endl;
	}
	for (const auto& pair : graph.edges)
	{
		if (finalKeptNodes.count(pair.first.id()) == 0) continue;
		for (const auto& target : pair.second)
		{
			if (finalKeptNodes.count(std::get<0>(target).id()) == 0) continue;
			file << "L\t" << graph.nodeNames[pair.first.id()] << "\t" << (pair.first.forward() ? "+" : "-") << "\t" << graph.nodeNames[std::get<0>(target).id()] << "\t" << (std::get<0>(target).forward() ? "+" : "-") << "\t" << std::get<1>(target) << "M\tec:i:" << std::get<2>(target) << std::endl;
		}
	}
}

void HandleCluster(const ClusterParams& params)
{
	std::filesystem::create_directories(params.basePath + "/tmp");
	std::cerr << "running MBG" << std::endl;
	runMBG(params.basePath, params.hifiReadPath, params.MBGPath, params.k, params.maxResolveLength, params.numThreads);
	std::cerr << "reading graph" << std::endl;
	GfaGraph graph;
	graph.loadFromFile(params.basePath + "/graph.gfa");
	std::cerr << "reading read paths" << std::endl;
	std::vector<ReadPath> readPaths = loadReadPaths(params.basePath + "/paths.gaf", graph);
	std::cerr << "getting consensus" << std::endl;
	Path heavyPath = getHeavyPath(graph, readPaths, params.maxResolveLength);
	std::cerr << "consensus length " << heavyPath.getSequence(graph.nodeSeqs).size() << "bp" << std::endl;
	if (params.orientReferencePath.size() > 0)
	{
		std::cerr << "orienting consensus" << std::endl;
		heavyPath = orientPath(graph, heavyPath, params.orientReferencePath, 101);
	}
	std::cerr << "writing consensus" << std::endl;
	writePathSequence(heavyPath, graph, params.basePath + "/consensus.fa", params.namePrefix);
	writePathGaf(heavyPath, graph, params.basePath + "/consensus_path.gaf");
	std::cerr << "process graph" << std::endl;
	processGraphAndWrite(heavyPath, graph, params.basePath + "/processed-graph.gfa");
	if (params.annotationFasta.size() > 0)
	{
		std::cerr << "lifting over annotations to consensus" << std::endl;
		liftoverAnnotationsToConsensus(params.basePath, params.basePath + "/consensus.fa", params.annotationFasta, params.annotationGff3, params.basePath + "/tmp", params.liftoffPath);
	}
}

void writeMorphPaths(const std::string& outputFile, const std::vector<MorphConsensus>& morphConsensuses, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip)
{
	std::ofstream file { outputFile };
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		file << morphConsensuses[i].name << "\t" << morphConsensuses[i].sequence.size() << "\t" << 0 << "\t" << morphConsensuses[i].sequence.size() << "\t+\t";
		for (auto node : morphConsensuses[i].path)
		{
			file << (node.forward() ? ">" : "<") << graph.nodeNames.at(node.id());
		}
		size_t size = getSequence(morphConsensuses[i].path, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges).size();
		file << "\t" << size << "\t" << pathStartClip.at(morphConsensuses[i].path[0]) << "\t" << (size - pathEndClip.at(morphConsensuses[i].path.back())) << "\t" << morphConsensuses[i].sequence.size() << "\t" << morphConsensuses[i].sequence.size() << "\t60" << std::endl;
	}
}

void writeMorphGraphAndReadPaths(const std::string& graphFile, const std::string& pathsFile, const std::vector<MorphConsensus>& morphConsensuses)
{
	std::unordered_map<std::string, size_t> originalReadLength;
	std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t, size_t>>> readFwMatches;
	std::unordered_map<std::string, std::vector<std::tuple<size_t, size_t, size_t>>> readBwMatches;
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		for (const auto& loop : morphConsensuses[i].ontLoops)
		{
			assert(loop.originalReadLength != 0);
			originalReadLength[loop.readName] = loop.originalReadLength;
			if (loop.approxEnd < loop.approxStart)
			{
				readBwMatches[loop.readName].emplace_back(loop.approxStart, loop.approxEnd, i);
			}
			else
			{
				readFwMatches[loop.readName].emplace_back(loop.approxStart, loop.approxEnd, i);
			}
		}
	}
	std::vector<ReadPath> readPaths;
	std::vector<size_t> pathLength;
	for (auto& pair : readFwMatches)
	{
		assert(pair.second.size() >= 1);
		std::sort(pair.second.begin(), pair.second.end(), [](auto left, auto right){ return std::get<0>(left) < std::get<0>(right); });
		for (size_t i = 0; i < pair.second.size(); i++)
		{
			if (i == 0 || std::get<0>(pair.second[i]) + 1000 < std::get<1>(pair.second[i-1]) || std::get<0>(pair.second[i]) > std::get<1>(pair.second[i-1]) + 1000)
			{
				readPaths.emplace_back();
				pathLength.emplace_back(0);
				readPaths.back().readLength = originalReadLength.at(pair.first);
				readPaths.back().readName = pair.first;
				readPaths.back().pathStartClip = 0;
				readPaths.back().pathEndClip = 0;
				readPaths.back().reverse = false;
				readPaths.back().readStart = std::get<0>(pair.second[i]);
				readPaths.back().readEnd = std::get<1>(pair.second[i]);
				readPaths.back().path.emplace_back(std::get<2>(pair.second[i]), true);
				pathLength.back() = morphConsensuses[std::get<2>(pair.second[i])].sequence.size();
			}
			else
			{
				readPaths.back().readEnd = std::get<1>(pair.second[i]);
				readPaths.back().path.emplace_back(std::get<2>(pair.second[i]), true);
				pathLength.back() += morphConsensuses[std::get<2>(pair.second[i])].sequence.size();
			}
		}
	}
	for (auto& pair : readBwMatches)
	{
		assert(pair.second.size() >= 1);
		std::sort(pair.second.begin(), pair.second.end(), [](auto left, auto right){ return std::get<0>(left) < std::get<0>(right); });
		for (size_t i = 0; i < pair.second.size(); i++)
		{
			if (i == 0 || std::get<1>(pair.second[i]) + 1000 < std::get<0>(pair.second[i-1]) || std::get<1>(pair.second[i]) > std::get<0>(pair.second[i-1]) + 1000)
			{
				readPaths.emplace_back();
				pathLength.emplace_back(0);
				readPaths.back().readLength = originalReadLength.at(pair.first);
				readPaths.back().readName = pair.first;
				readPaths.back().pathStartClip = 0;
				readPaths.back().pathEndClip = 0;
				readPaths.back().reverse = false;
				readPaths.back().readStart = std::get<1>(pair.second[i]);
				readPaths.back().readEnd = std::get<0>(pair.second[i]);
				readPaths.back().path.emplace_back(std::get<2>(pair.second[i]), false);
				pathLength.back() = morphConsensuses[std::get<2>(pair.second[i])].sequence.size();
			}
			else
			{
				readPaths.back().readEnd = std::get<0>(pair.second[i]);
				readPaths.back().path.emplace_back(std::get<2>(pair.second[i]), false);
				pathLength.back() += morphConsensuses[std::get<2>(pair.second[i])].sequence.size();
			}
		}
	}
	std::map<std::pair<Node, Node>, size_t> edgeCoverage;
	for (const auto& path : readPaths)
	{
		for (size_t i = 1; i < path.path.size(); i++)
		{
			auto key = canon(path.path[i-1], path.path[i]);
			edgeCoverage[key] += 1;
		}
	}
	{
		std::ofstream file { graphFile };
		for (size_t i = 0; i < morphConsensuses.size(); i++)
		{
			file << "S\t" << morphConsensuses[i].name << "\t*\tLN:i:" << morphConsensuses[i].sequence.size() << "\tll:f:" << morphConsensuses[i].coverage << "\tFC:i:" << (morphConsensuses[i].coverage * morphConsensuses[i].sequence.size()) << std::endl;
		}
		for (auto pair : edgeCoverage)
		{
			file << "L\t" << morphConsensuses[pair.first.first.id()].name << "\t" << (pair.first.first.forward() ? "+" : "-") << "\t" << morphConsensuses[pair.first.second.id()].name << "\t" << (pair.first.second.forward() ? "+" : "-") << "\t0M\tec:i:" << pair.second << std::endl;
		}
	}
	{
		std::ofstream file { pathsFile };
		for (size_t i = 0; i < readPaths.size(); i++)
		{
			file << readPaths[i].readName << "\t" << readPaths[i].readLength << " \t" << readPaths[i].readStart << "\t" << readPaths[i].readEnd << "\t+\t";
			for (auto node : readPaths[i].path)
			{
				file << (node.forward() ? ">" : "<") << morphConsensuses[node.id()].name;
			}
			file << "\t" << pathLength[i] << "\t" << 0 << "\t" << pathLength[i] << "\t" << pathLength[i] << "\t" << pathLength[i] << "\t" << 60 << std::endl;
		}
	}
}

std::pair<size_t, size_t> refineStartEndPoses(const std::string& readSequence, const std::string& consensusSequence, const phmap::flat_hash_map<uint64_t, size_t>& consensusKmerPositions, const size_t refiningK, const size_t initialStartPos, const size_t initialEndPos)
{
	const uint64_t mask = (1ull << (2ull*refiningK)) - 1ull;
	size_t checkAreaStart = initialStartPos;
	if (checkAreaStart > 5000)
	{
		checkAreaStart -= 5000;
	}
	else
	{
		checkAreaStart = 0;
	}
	size_t checkAreaEnd = std::min(initialEndPos + 5000, readSequence.size());
	uint64_t kmer = 0;
	for (size_t i = 0; i < refiningK; i++)
	{
		kmer <<= 2;
		switch(readSequence[checkAreaStart+i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
	}
	std::vector<std::pair<size_t, size_t>> matches;
	if (consensusKmerPositions.count(kmer) == 1)
	{
		matches.emplace_back(checkAreaStart, consensusKmerPositions.at(kmer));
	}
	for (size_t i = refiningK; checkAreaStart+i < checkAreaEnd; i++)
	{
		kmer <<= 2;
		kmer &= mask;
		switch(readSequence[checkAreaStart+i])
		{
			case 'A':
				kmer += 0;
				break;
			case 'C':
				kmer += 1;
				break;
			case 'G':
				kmer += 2;
				break;
			case 'T':
				kmer += 3;
				break;
			default:
				assert(false);
				break;
		}
		if (consensusKmerPositions.count(kmer) == 1)
		{
			matches.emplace_back(checkAreaStart+i-refiningK+1, consensusKmerPositions.at(kmer));
		}
	}
	removeNonDiagonalKmers(matches);
	removeRepeatKmersEitherRepeat(matches);
	matches = getIncreasingChain(matches);
	if (matches.size() < 100)
	{
		return std::make_pair(initialStartPos, initialEndPos);
	}
	if (matches[0].second > 200)
	{
		return std::make_pair(initialStartPos, initialEndPos);
	}
	if (matches.back().second+200 < consensusSequence.size())
	{
		return std::make_pair(initialStartPos, initialEndPos);
	}
	size_t approxStart = matches[0].first;
	if (matches[0].second < approxStart)
	{
		approxStart -= matches[0].second;
	}
	else
	{
		approxStart = 0;
	}
	size_t approxEnd = std::min(matches.back().first + (consensusSequence.size() - matches.back().second), readSequence.size());
	approxStart = getExactBreakPos(readSequence, consensusSequence, approxStart);
	approxEnd = getExactBreakPos(readSequence, consensusSequence, approxEnd);
	return std::make_pair(approxStart, approxEnd);
}

std::vector<std::vector<std::string>> getRawLoopSequences(const std::vector<std::vector<OntLoop>>& clusters, const std::string& rawOntPath, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t numThreads)
{
	const size_t refiningK = 15;
	phmap::flat_hash_map<std::string, std::vector<std::tuple<size_t, size_t, size_t, size_t>>> loopPositionsInReads;
	std::vector<phmap::flat_hash_map<uint64_t, size_t>> clusterConsensusKmerPositions;
	std::vector<std::string> clusterConsensuses;
	for (size_t i = 0; i < clusters.size(); i++)
	{
		auto path = getConsensusPath(clusters[i], graph);
		clusterConsensuses.emplace_back(getConsensusSequence(path, graph, pathStartClip, pathEndClip));
		clusterConsensusKmerPositions.emplace_back(getRefKmers(clusterConsensuses.back(), refiningK));
	}
	std::vector<std::vector<std::string>> result;
	result.resize(clusters.size());
	for (size_t i = 0; i < clusters.size(); i++)
	{
		result[i].resize(clusters[i].size(), "");
		for (size_t j = 0; j < clusters[i].size(); j++)
		{
			loopPositionsInReads[clusters[i][j].readName].emplace_back(i, j, clusters[i][j].approxStart, clusters[i][j].approxEnd);
		}
	}
	std::vector<std::thread> threads;
	std::atomic<bool> readDone;
	readDone = false;
	std::mutex stackMutex;
	std::vector<FastQ> readStack;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&readDone, &stackMutex, &readStack, &result, &clusterConsensusKmerPositions, &clusterConsensuses, &loopPositionsInReads]()
		{
			while (true)
			{
				FastQ fastq;
				{
					std::lock_guard<std::mutex> lock { stackMutex };
					if (readStack.size() == 0)
					{
						if (readDone)
						{
							if (readStack.size() == 0)
							{
								break;
							}
						}
						continue;
					}
					std::swap(readStack.back(), fastq);
					readStack.pop_back();
				}
				std::string readname = nameWithoutTags(fastq.seq_id);
				if (loopPositionsInReads.count(readname) == 0) continue;
				std::string revcompRead = revcomp(fastq.sequence);
				for (auto t : loopPositionsInReads.at(readname))
				{
					size_t startPos = std::get<2>(t);
					size_t endPos = std::get<3>(t);
					bool reverse = false;
					if (startPos > endPos)
					{
						reverse = true;
					}
					if (reverse)
					{
						startPos = fastq.sequence.size() - startPos;
						endPos = fastq.sequence.size() - endPos;
						assert(endPos > startPos);
						std::tie(startPos, endPos) = refineStartEndPoses(revcompRead, clusterConsensuses[std::get<0>(t)], clusterConsensusKmerPositions[std::get<0>(t)], refiningK, startPos, endPos);
						startPos = fastq.sequence.size() - startPos;
						endPos = fastq.sequence.size() - endPos;
						std::swap(startPos, endPos);
						assert(endPos > startPos);
					}
					else
					{
						std::tie(startPos, endPos) = refineStartEndPoses(fastq.sequence, clusterConsensuses[std::get<0>(t)], clusterConsensusKmerPositions[std::get<0>(t)], refiningK, startPos, endPos);
					}
					if (startPos > 500)
					{
						startPos -= 500;
					}
					else
					{
						startPos = 0;
					}
					if (endPos+500 < fastq.sequence.size())
					{
						endPos += 500;
					}
					else
					{
						endPos = fastq.sequence.size();
					}
					std::string str = fastq.sequence.substr(startPos, endPos-startPos);
					if (reverse)
					{
						str = revcomp(str);
					}
					assert(result[std::get<0>(t)][std::get<1>(t)] == "");
					assert(str.size() >= 1);
					result[std::get<0>(t)][std::get<1>(t)] = str;
				}
			}
		});
	}
	FastQ::streamFastqFromFile(rawOntPath, false, [&readStack, &stackMutex](FastQ& fastq)
	{
		std::lock_guard<std::mutex> lock { stackMutex };
		readStack.emplace_back();
		std::swap(readStack.back(), fastq);
	});
	readDone = true;
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t j = 0; j < result[i].size(); j++)
		{
			assert(result[i][j] != "");
		}
	}
	return result;
}

void writeSelfCorrectedOntLoopSequences(const std::string outputFile, const std::vector<std::vector<OntLoop>>& loopSequences)
{
	std::ofstream file { outputFile };
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		for (size_t j = 0; j < loopSequences[i].size(); j++)
		{
			file << ">" << loopSequences[i][j].rawLoopName << std::endl;
			file << loopSequences[i][j].selfCorrectedSequence << std::endl;
		}
	}
}

void writeRawOntLoopSequences(const std::string outputFile, const std::vector<std::vector<OntLoop>>& loopSequences)
{
	std::ofstream file { outputFile };
	for (size_t i = 0; i < loopSequences.size(); i++)
	{
		for (size_t j = 0; j < loopSequences[i].size(); j++)
		{
			file << ">" << loopSequences[i][j].rawLoopName << std::endl;
			file << loopSequences[i][j].rawSequence << std::endl;
		}
	}
}

void alignSelfCorrectedOntLoopsToMorphConsensuses(const std::vector<std::vector<OntLoop>>& loopSequences, std::string outputFile, std::string tmppath, const size_t numThreads, const std::vector<MorphConsensus>& morphConsensuses, const std::string& winnowmapPath, const std::string& samtoolsPath)
{
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		{
			std::ofstream tmpReadsFile { tmppath + "/tmpreads.fa" };
			std::ofstream tmpConsensusFile { tmppath + "/tmpconsensus.fa" };
			for (size_t j = 0; j < loopSequences[i].size(); j++)
			{
				tmpReadsFile << ">" << loopSequences[i][j].rawLoopName << std::endl;
				tmpReadsFile << loopSequences[i][j].selfCorrectedSequence << std::endl;
			}
			tmpConsensusFile << ">" << morphConsensuses[i].name << std::endl;
			tmpConsensusFile << morphConsensuses[i].sequence << std::endl;
		}
		std::string winnowmapCommand = winnowmapPath + " -x map-ont -a -t " + std::to_string(numThreads) + " " + tmppath + "/tmpconsensus.fa " + tmppath + "/tmpreads.fa | " + samtoolsPath + " view -b | " + samtoolsPath + " sort > " + tmppath + "/tmpalns_" + std::to_string(i) + ".bam";
		std::cerr << "winnowmap command:" << std::endl;
		std::cerr << winnowmapCommand << std::endl;
		int result = system(winnowmapCommand.c_str());
		if (result != 0)
		{
			std::cerr << "alignment did not run successfully" << std::endl;
			std::abort();
		}
	}
	std::string mergeCommand = samtoolsPath + " merge -o " + outputFile;
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		mergeCommand += " " + tmppath + "/tmpalns_" + std::to_string(i) + ".bam";
	}
	std::cerr << "samtools command:" << std::endl;
	std::cerr << mergeCommand << std::endl;
	int result = system(mergeCommand.c_str());
	if (result != 0)
	{
		std::cerr << "samtools merge did not run successfully" << std::endl;
		std::abort();
	}
	std::string indexCommand = samtoolsPath + " index -b " + outputFile;
	std::cerr << "samtools command:" << std::endl;
	std::cerr << indexCommand << std::endl;
	result = system(indexCommand.c_str());
	if (result != 0)
	{
		std::cerr << "samtools index did not run successfully" << std::endl;
		std::abort();
	}
}

void alignRawOntLoopsToMorphConsensuses(const std::vector<std::vector<OntLoop>>& loopSequences, std::string outputFile, std::string tmppath, const size_t numThreads, const std::vector<MorphConsensus>& morphConsensuses, const std::string& winnowmapPath, const std::string& samtoolsPath)
{
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		{
			std::ofstream tmpReadsFile { tmppath + "/tmpreads.fa" };
			std::ofstream tmpConsensusFile { tmppath + "/tmpconsensus.fa" };
			for (size_t j = 0; j < loopSequences[i].size(); j++)
			{
				tmpReadsFile << ">" << loopSequences[i][j].rawLoopName << std::endl;
				tmpReadsFile << loopSequences[i][j].rawSequence << std::endl;
			}
			tmpConsensusFile << ">" << morphConsensuses[i].name << std::endl;
			tmpConsensusFile << morphConsensuses[i].sequence << std::endl;
		}
		std::string winnowmapCommand = winnowmapPath + " -x map-ont -a -t " + std::to_string(numThreads) + " " + tmppath + "/tmpconsensus.fa " + tmppath + "/tmpreads.fa | " + samtoolsPath + " view -b | " + samtoolsPath + " sort > " + tmppath + "/tmpalns_" + std::to_string(i) + ".bam";
		std::cerr << "winnowmap command:" << std::endl;
		std::cerr << winnowmapCommand << std::endl;
		int result = system(winnowmapCommand.c_str());
		if (result != 0)
		{
			std::cerr << "alignment did not run successfully" << std::endl;
			std::abort();
		}
	}
	std::string mergeCommand = samtoolsPath + " merge -o " + outputFile;
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		mergeCommand += " " + tmppath + "/tmpalns_" + std::to_string(i) + ".bam";
	}
	std::cerr << "samtools command:" << std::endl;
	std::cerr << mergeCommand << std::endl;
	int result = system(mergeCommand.c_str());
	if (result != 0)
	{
		std::cerr << "samtools merge did not run successfully" << std::endl;
		std::abort();
	}
	std::string indexCommand = samtoolsPath + " index -b " + outputFile;
	std::cerr << "samtools command:" << std::endl;
	std::cerr << indexCommand << std::endl;
	result = system(indexCommand.c_str());
	if (result != 0)
	{
		std::cerr << "samtools index did not run successfully" << std::endl;
		std::abort();
	}
}

void addRawSequencesToLoops(std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::string ontReadPath, const size_t numThreads)
{
	std::vector<std::vector<std::string>> sequences = getRawLoopSequences(clusters, ontReadPath, graph, pathStartClip, pathEndClip, numThreads);
	assert(sequences.size() == clusters.size());
	for (size_t i = 0; i < clusters.size(); i++)
	{
		assert(sequences[i].size() == clusters[i].size());
		for (size_t j = 0; j < clusters[i].size(); j++)
		{
			clusters[i][j].rawSequence = sequences[i][j];
		}
	}
}

void addRawSequenceNamesToLoops(std::vector<std::vector<OntLoop>>& clusters)
{
	for (size_t i = 0; i < clusters.size(); i++)
	{
		for (size_t j = 0; j < clusters[i].size(); j++)
		{
			clusters[i][j].rawLoopName = "morph_" + std::to_string(i) + "_rawloop_" + std::to_string(j) + "_" + clusters[i][j].readName + "_" + std::to_string(clusters[i][j].approxStart) + "_" + std::to_string(clusters[i][j].approxEnd);
		}
	}
}

void addSelfCorrectedOntLoopSequences(std::vector<std::vector<OntLoop>>& clusters)
{
	for (size_t i = 0; i < clusters.size(); i++)
	{
		std::cerr << "self-correct cluster " << i << " with size " << clusters[i].size() << std::endl;
		auto selfCorrectedLoops = getSelfCorrectedLoops(clusters[i]);
		assert(selfCorrectedLoops.size() == clusters[i].size());
		for (size_t j = 0; j < clusters[i].size(); j++)
		{
			clusters[i][j].selfCorrectedSequence = selfCorrectedLoops[j];
			std::cerr << "corrected read " << j << " from size " << clusters[i][j].rawSequence.size() << " to " << clusters[i][j].selfCorrectedSequence.size() << std::endl;
		}
	}
}

void DoClusterONTAnalysis(const ClusterParams& params)
{
	std::cerr << "reading allele graph" << std::endl;
	GfaGraph graph;
	graph.loadFromFile(params.basePath + "/processed-graph.gfa");
	std::cerr << "reading consensus" << std::endl;
	Path heavyPath = readHeavyPath(graph, params.basePath + "/consensus_path.gaf");
	std::cerr << "extract corrected ultralong paths" << std::endl;
	size_t heavyPathLength = heavyPath.getSequence(graph.nodeSeqs).size();
	size_t minLength = heavyPathLength * 0.5;
	std::cerr << "consensus path length " << heavyPathLength << ", using " << minLength << " as minimum morph length" << std::endl;
	auto ontPaths = extractCorrectedONTPaths(params.basePath + "/ont-alns.gaf", heavyPath, minLength, graph);
	std::cerr << ontPaths.size() << " corrected paths" << std::endl;
	std::cerr << "extract loops from ONTs" << std::endl;
	std::unordered_map<Node, size_t> pathStartClip;
	std::unordered_map<Node, size_t> pathEndClip;
	std::unordered_set<size_t> borderNodes;
	std::unordered_set<size_t> anchorNodes;
	std::tie(borderNodes, anchorNodes, pathStartClip, pathEndClip) = getBorderNodes(heavyPath, graph);
	assert(borderNodes.size() > 0);
	auto loopSequences = extractLoopSequences(ontPaths, heavyPath, minLength, graph, borderNodes, anchorNodes, pathStartClip, pathEndClip);
	auto coreNodes = getCoreNodes(loopSequences);
	std::cerr << loopSequences.size() << " loops in ONTs" << std::endl;
	{
		std::ofstream file { params.basePath + "/tmp/loops.fa" };
		for (size_t i = 0; i < loopSequences.size(); i++)
		{
			file << ">loop_" << i << "_read_" << loopSequences[i].readName << "_start_" << loopSequences[i].approxStart << "_end_" << loopSequences[i].approxEnd << std::endl;
			std::string loopSeq = getSequence(loopSequences[i].path, graph.nodeSeqs, graph.revCompNodeSeqs, graph.edges);
			size_t leftClipBp = pathStartClip.at(loopSequences[i].path[0]);
			size_t rightClipBp = pathEndClip.at(loopSequences[i].path.back());
			loopSeq = loopSeq.substr(leftClipBp, loopSeq.size() - leftClipBp - rightClipBp);
			file << loopSeq << std::endl;
		}
	}
	std::cerr << "cluster loops roughly" << std::endl;
	orderLoopsByLength(loopSequences, graph, pathStartClip, pathEndClip);
	std::cerr << "max clustering edit distance " << params.maxClusterDifference << std::endl;
	auto clusters = roughClusterLoopSequences(loopSequences, graph, pathStartClip, pathEndClip, coreNodes, params.maxClusterDifference);
	std::cerr << clusters.size() << " rough clusters" << std::endl;
	std::cerr << "getting exact locations of raw loop sequences" << std::endl;
	addRawSequencesToLoops(clusters, graph, pathStartClip, pathEndClip, params.ontReadPath, params.numThreads);
	if (params.extraPhasing)
	{
		std::cerr << "phase clusters by raw sequences" << std::endl;
		clusters = editSplitClusters(clusters, graph, pathStartClip, pathEndClip, params.numThreads, params.MBGPath, params.basePath + "/tmp");
	}
	else
	{
		std::cerr << "cluster loops by density" << std::endl;
		clusters = densityClusterLoops(clusters, graph, pathStartClip, pathEndClip, coreNodes, params.maxClusterDifference, 5, params.minReclusterDistance);
	}
	std::cerr << clusters.size() << " clusters" << std::endl;
	std::sort(clusters.begin(), clusters.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
	addRawSequenceNamesToLoops(clusters);
	std::cerr << "write raw ONT loop sequences" << std::endl;
	writeRawOntLoopSequences(params.basePath + "/raw_loops.fa", clusters);
	std::cerr << "get self-corrected ONT loop sequences" << std::endl;
	addSelfCorrectedOntLoopSequences(clusters);
	std::cerr << "write self-corrected ONT loop sequences" << std::endl;
	writeSelfCorrectedOntLoopSequences(params.basePath + "/tmp/loops_selfcorrected.fa", clusters);
	std::cerr << "getting morph consensuses" << std::endl;
	auto morphConsensuses = getMorphConsensuses(clusters, graph, pathStartClip, pathEndClip, params.namePrefix);
	std::cerr << "polishing morph consensuses" << std::endl;
	polishMorphConsensuses(morphConsensuses, clusters, params.MBGPath, params.basePath + "/tmp", params.numThreads);
	std::cerr << "write morph consensuses" << std::endl;
	writeMorphConsensuses(params.basePath + "/morphs.fa", params.basePath + "/tmp/morphs_preconsensus.fa", morphConsensuses);
	std::cerr << "write morph paths" << std::endl;
	writeMorphPaths(params.basePath + "/morphs.gaf", morphConsensuses, graph, pathStartClip, pathEndClip);
	std::cerr << "write morph graph and read paths" << std::endl;
	writeMorphGraphAndReadPaths(params.basePath + "/morphgraph.gfa", params.basePath + "/readpaths-morphgraph.gaf", morphConsensuses);
	if (params.winnowmapPath != "" && params.samtoolsPath != "")
	{
		std::cerr << "realign raw ONT loop sequences to morph consensuses" << std::endl;
		alignRawOntLoopsToMorphConsensuses(clusters, params.basePath + "/raw_loop_to_morphs_alignments.bam", params.basePath + "/tmp", params.numThreads, morphConsensuses, params.winnowmapPath, params.samtoolsPath);
		std::cerr << "realign self-corrected ONT loop sequences to morph consensuses" << std::endl;
		alignSelfCorrectedOntLoopsToMorphConsensuses(clusters, params.basePath + "/tmp/selfcorrected_loop_to_morphs_alignments.bam", params.basePath + "/tmp", params.numThreads, morphConsensuses, params.winnowmapPath, params.samtoolsPath);
	}
	else
	{
		if (params.winnowmapPath == "")
		{
			std::cerr << "winnowmap not found, ";
		}
		if (params.samtoolsPath != "")
		{
			std::cerr << "bamtools not found, ";
		}
		std::cerr << "skipping alignment of raw ONT loop sequences to morph consensuses" << std::endl;
	}
	if (params.annotationFasta.size() > 0)
	{
		std::cerr << "lifting over annotations to morphs" << std::endl;
		liftoverAnnotationsToMorphs(params.basePath, morphConsensuses, params.annotationFasta, params.annotationGff3, params.basePath + "/tmp", params.liftoffPath);
	}
}
