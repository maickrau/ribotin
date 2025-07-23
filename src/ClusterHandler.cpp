#include <limits>
#include <queue>
#include <map>
#include <cstdlib>
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
#include "LoopAnalysis.h"
#include "Logger.h"
#include "KmerMatcher.h"

void runMBG(std::string basePath, std::string readPath, std::string MBGPath, const size_t k, const size_t maxResolveLength, const size_t numThreads)
{
	std::string mbgCommand;
	mbgCommand = MBGPath + " -o " + basePath + "/graph.gfa -t " + std::to_string(numThreads) + " -i " + readPath + " -k " + std::to_string(k) + " -w " + std::to_string(k-30) + " -a 2 -u 3 -r " + std::to_string(maxResolveLength) + " -R 4000 --error-masking=msat --output-sequence-paths " + basePath + "/paths.gaf --only-local-resolve --resolve-palindromes-global 1> " + basePath + "/tmp/mbg_stdout.txt 2> " + basePath + "/tmp/mbg_stderr.txt";
	Logger::Log.log(Logger::LogLevel::Always) << "MBG command:" << std::endl;
	Logger::Log.log(Logger::LogLevel::Always) << mbgCommand << std::endl;
	int result = system(mbgCommand.c_str());
	if (result != 0)
	{
		Logger::Log.log(Logger::LogLevel::Always) << "MBG did not run successfully" << std::endl;
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
	Logger::Log.log(Logger::LogLevel::Always) << "running liftoff with command:" << std::endl;
	Logger::Log.log(Logger::LogLevel::Always) << command << std::endl;
	int result = system(command.c_str());
	if (result != 0)
	{
		Logger::Log.log(Logger::LogLevel::Always) << "liftoff did not run successfully" << std::endl;
		std::abort();
	}
}

void AlignONTReads(std::string basePath, std::string graphAlignerPath, std::string ontReadPath, std::string graphPath, std::string outputAlnPath, size_t numThreads)
{
	std::string graphalignerCommand;
	graphalignerCommand = graphAlignerPath + " -g " + graphPath + " -f " + ontReadPath + " -a " + outputAlnPath + " -t " + std::to_string(numThreads) + " --seeds-mxm-length 30 --seeds-mem-count 10000 --bandwidth 15 --multimap-score-fraction 0.99 --precise-clipping 0.85 --min-alignment-score 5000 --discard-cigar --clip-ambiguous-ends 100 --overlap-incompatible-cutoff 0.15 --mem-index-no-wavelet-tree --max-trace-count 5 1> " + basePath + "/tmp/graphaligner_stdout.txt 2> " + basePath + "/tmp/graphaligner_stderr.txt";
	Logger::Log.log(Logger::LogLevel::Always) << "GraphAligner command:" << std::endl;
	Logger::Log.log(Logger::LogLevel::Always) << graphalignerCommand << std::endl;
	int result = system(graphalignerCommand.c_str());
	if (result != 0)
	{
		Logger::Log.log(Logger::LogLevel::Always) << "GraphAligner did not run successfully" << std::endl;
		std::abort();
	}
}

void processGraphAndWrite(const Path& heavyPath, const GfaGraph& graph, const std::vector<ReadPath>& readPaths, const std::string outputGraphName)
{
	const size_t minCoverage = 5;
	const size_t minAnchorLength = 5000;
	phmap::flat_hash_set<size_t> keptNodesInCycles = getNodesInMainCycle(graph, heavyPath, minCoverage);
	phmap::flat_hash_map<Node, size_t> shortestDistanceToCyclicComponent = getNodeDistancesToMainCycle(graph, keptNodesInCycles);
	phmap::flat_hash_set<std::vector<Node>> distinctNoncyclePaths;
	phmap::flat_hash_map<std::vector<Node>, size_t> pathCoverage;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		size_t firstCyclic = std::numeric_limits<size_t>::max();
		size_t lastCyclic = std::numeric_limits<size_t>::max();
		for (size_t j = 0; j < readPaths[i].path.size(); j++)
		{
			if (keptNodesInCycles.count(readPaths[i].path[j].id()) == 1)
			{
				if (firstCyclic == std::numeric_limits<size_t>::max()) firstCyclic = j;
				lastCyclic = j;
			}
			if (shortestDistanceToCyclicComponent.count(readPaths[i].path[j]) == 1 && shortestDistanceToCyclicComponent.count(reverse(readPaths[i].path[j])) == 1)
			{
				if (firstCyclic == std::numeric_limits<size_t>::max()) firstCyclic = j;
				lastCyclic = j;
			}
		}
		if (firstCyclic == std::numeric_limits<size_t>::max() && lastCyclic == std::numeric_limits<size_t>::max()) continue;
		if (firstCyclic > 0)
		{
			std::vector<Node> path;
			path.insert(path.end(), readPaths[i].path.begin(), readPaths[i].path.begin()+firstCyclic);
			std::reverse(path.begin(), path.end());
			for (size_t j = 0; j < path.size(); j++)
			{
				path[j] = reverse(path[j]);
			}
			if (getPathLength(path, graph.nodeSeqs, graph.edges) > minAnchorLength)
			{
				path.insert(path.begin(), reverse(readPaths[i].path[firstCyclic]));
				pathCoverage[path] += 1;
				distinctNoncyclePaths.emplace(path);
			}
		}
		if (lastCyclic+1 < readPaths[i].path.size())
		{
			std::vector<Node> path;
			path.insert(path.end(), readPaths[i].path.begin()+lastCyclic+1, readPaths[i].path.end());
			if (getPathLength(path, graph.nodeSeqs, graph.edges) > minAnchorLength)
			{
				path.insert(path.begin(), readPaths[i].path[lastCyclic]);
				pathCoverage[path] += 1;
				distinctNoncyclePaths.emplace(path);
			}
		}
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << distinctNoncyclePaths.size() << " distinct noncycle paths" << std::endl;
	phmap::flat_hash_set<std::vector<Node>> distinctNoncyclePathsLengthFilter;
	for (const auto& path : distinctNoncyclePaths)
	{
		size_t clipCount = 0;
		std::vector<Node> pathSubstr;
		pathSubstr.insert(pathSubstr.end(), path.begin()+1, path.end());
		assert(pathSubstr.size() >= 1);
		size_t pathLength = getPathLength(pathSubstr, graph.nodeSeqs, graph.edges);
		while (pathLength > minAnchorLength + graph.nodeSeqs[pathSubstr.back().id()].size())
		{
			pathSubstr.pop_back();
			pathLength = getPathLength(pathSubstr, graph.nodeSeqs, graph.edges);
			assert(pathSubstr.size() >= 1);
			clipCount += 1;
		}
		pathSubstr.insert(pathSubstr.begin(), path[0]);
		distinctNoncyclePathsLengthFilter.insert(pathSubstr);
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << distinctNoncyclePathsLengthFilter.size() << " distinct noncycle paths length filtered" << std::endl;
	phmap::flat_hash_set<std::vector<Node>> containedPaths;
	for (const auto& path : distinctNoncyclePathsLengthFilter)
	{
		for (const auto& path2 : distinctNoncyclePathsLengthFilter)
		{
			if (path == path2) continue;
			if (path2.size() >= path.size()) continue;
			bool contained = true;
			for (size_t i = 0; i < path2.size(); i++)
			{
				if (path2[i] != path[i]) contained = false;
			}
			if (contained) containedPaths.emplace(path2);
		}
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << containedPaths.size() << " contained paths" << std::endl;
	std::vector<std::vector<Node>> borderPaths;
	for (const auto& path : distinctNoncyclePathsLengthFilter)
	{
		if (containedPaths.count(path) == 1) continue;
		borderPaths.emplace_back(path);
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "path ";
		for (size_t i = 0; i < path.size(); i++)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << (path[i].forward() ? ">" : "<") << graph.nodeNames[path[i].id()];
		}
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "path (revcomp) ";
		for (size_t i = 0; i < path.size(); i++)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << (path[path.size()-1-i].forward() ? "<" : ">") << graph.nodeNames[path[path.size()-1-i].id()];
		}
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
	}
	std::vector<size_t> borderPathCoverages;
	for (size_t i = 0; i < borderPaths.size(); i++)
	{
		borderPathCoverages.emplace_back(0);
		for (const auto& pair : pathCoverage)
		{
			if (pair.first.size() < borderPaths[i].size()) continue;
			assert(pair.first.size() >= borderPaths[i].size());
			bool contained = true;
			for (size_t j = 0; j < borderPaths[i].size(); j++)
			{
				if (pair.first[j] != borderPaths[i][j]) contained = false;
			}
			if (contained) borderPathCoverages.back() += pair.second;
		}
	}
	std::vector<std::string> borderPathSequences;
	borderPathSequences.resize(borderPaths.size());
	for (size_t i = 0; i < borderPaths.size(); i++)
	{
		assert(borderPaths[i].size() >= 2);
		for (size_t j = 1; j < borderPaths[i].size(); j++)
		{
			std::string add = graph.nodeSeqs[borderPaths[i][j].id()];
			if (!borderPaths[i][j].forward())
			{
				add = revcomp(add);
			}
			if (j >= 2)
			{
				add = add.substr(getOverlap(borderPaths[i][j-1], borderPaths[i][j], graph.edges));
			}
			borderPathSequences[i] += add;
		}
	}
	size_t shortestBorderPath = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < borderPathSequences.size(); i++)
	{
		shortestBorderPath = std::min(shortestBorderPath, borderPathSequences[i].size());
	}
	Logger::Log.log(Logger::LogLevel::DebugInfo) << "shortest border path " << shortestBorderPath << std::endl;
	assert(shortestBorderPath > 0);
	for (size_t i = 0; i < borderPathSequences.size(); i++)
	{
		assert(borderPathSequences[i].size() >= shortestBorderPath);
		borderPathSequences[i] = borderPathSequences[i].substr(0, shortestBorderPath);
	}
	phmap::flat_hash_set<size_t> finalKeptNodes;
	finalKeptNodes.insert(keptNodesInCycles.begin(), keptNodesInCycles.end());
	std::ofstream file { outputGraphName };
	for (size_t i = 0; i < borderPaths.size(); i++)
	{
		file << "S\tborder" << i << "\t" << borderPathSequences[i] << "\tll:f:" << borderPathCoverages[i] << "\tFC:i:" << (borderPathSequences[i].size() * borderPathCoverages[i]) << std::endl;
		file << "L\t" << graph.nodeNames[borderPaths[i][0].id()] << "\t" << (borderPaths[i][0].forward() ? "+" : "-") << "\tborder" << i << "\t+\t" << getOverlap(borderPaths[i][0], borderPaths[i][1], graph.edges) << "M\tec:i:" << borderPathCoverages[i] << std::endl;
	}
	for (size_t node : finalKeptNodes)
	{
		std::string sequence = graph.nodeSeqs[node];
		file << "S\t" << graph.nodeNames[node] << "\t" << sequence << "\tll:f:" << graph.nodeCoverages[node] << "\tFC:i:" << (graph.nodeCoverages[node] * sequence.size()) << std::endl;
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
	Logger::Log.log(Logger::LogLevel::Always) << "running MBG" << std::endl;
	runMBG(params.basePath, params.hifiReadPath, params.MBGPath, params.k, params.maxResolveLength, params.numThreads);
	Logger::Log.log(Logger::LogLevel::Always) << "reading graph" << std::endl;
	GfaGraph graph;
	graph.loadFromFile(params.basePath + "/graph.gfa");
	Logger::Log.log(Logger::LogLevel::Always) << "reading read paths" << std::endl;
	std::vector<ReadPath> readPaths = loadReadPaths(params.basePath + "/paths.gaf", graph);
	Logger::Log.log(Logger::LogLevel::Always) << "getting consensus" << std::endl;
	Path heavyPath = getHeavyPath(graph, readPaths, params.maxResolveLength);
	Logger::Log.log(Logger::LogLevel::Always) << "consensus length " << heavyPath.getSequence(graph.nodeSeqs).size() << "bp" << std::endl;
	if (params.orientReferencePath.size() > 0)
	{
		Logger::Log.log(Logger::LogLevel::Always) << "orienting consensus" << std::endl;
		heavyPath = orientPath(graph, heavyPath, params.orientReferencePath, 101);
	}
	Logger::Log.log(Logger::LogLevel::Always) << "writing consensus" << std::endl;
	writePathSequence(heavyPath, graph, params.basePath + "/consensus.fa", params.namePrefix);
	writePathGaf(heavyPath, graph, params.basePath + "/consensus_path.gaf");
	Logger::Log.log(Logger::LogLevel::Always) << "process graph" << std::endl;
	processGraphAndWrite(heavyPath, graph, readPaths, params.basePath + "/processed-graph.gfa");
	if (params.annotationFasta.size() > 0)
	{
		Logger::Log.log(Logger::LogLevel::Always) << "lifting over annotations to consensus" << std::endl;
		liftoverAnnotationsToConsensus(params.basePath, params.basePath + "/consensus.fa", params.annotationFasta, params.annotationGff3, params.basePath + "/tmp", params.liftoffPath);
	}
}

std::vector<MorphConsensus> GetONTClusters(const ClusterParams& params)
{
	Logger::Log.log(Logger::LogLevel::Always) << "reading allele graph" << std::endl;
	GfaGraph graph;
	graph.loadFromFile(params.basePath + "/processed-graph.gfa");
	Logger::Log.log(Logger::LogLevel::Always) << "reading consensus" << std::endl;
	Path heavyPath = readHeavyPath(graph, params.basePath + "/consensus_path.gaf");
	Logger::Log.log(Logger::LogLevel::Always) << "extract corrected ultralong paths" << std::endl;
	size_t heavyPathLength = heavyPath.getSequence(graph.nodeSeqs).size();
	size_t minLength = heavyPathLength * 0.5;
	Logger::Log.log(Logger::LogLevel::Always) << "consensus path length " << heavyPathLength << ", using " << minLength << " as minimum morph length" << std::endl;
	auto ontPaths = extractCorrectedONTPaths(params.basePath + "/ont-alns.gaf", heavyPath, minLength, graph);
	Logger::Log.log(Logger::LogLevel::Always) << ontPaths.size() << " corrected paths" << std::endl;
	Logger::Log.log(Logger::LogLevel::Always) << "augment broken paths" << std::endl;
	augmentGraphWithFakeBreaks(ontPaths, graph, params.ontReadPath);
	Logger::Log.log(Logger::LogLevel::Always) << ontPaths.size() << " merged paths" << std::endl;
	Logger::Log.log(Logger::LogLevel::Always) << "extract loops from ONTs" << std::endl;
	std::unordered_map<Node, size_t> pathStartClip;
	std::unordered_map<Node, size_t> pathEndClip;
	std::unordered_set<size_t> borderNodes;
	std::unordered_set<size_t> anchorNodes;
	std::tie(borderNodes, anchorNodes, pathStartClip, pathEndClip) = getBorderNodes(heavyPath, graph);
	assert(borderNodes.size() > 0);
	auto loopSequences = extractLoopSequences(ontPaths, heavyPath, minLength, graph, borderNodes, anchorNodes, pathStartClip, pathEndClip);
	auto coreNodes = getCoreNodes(loopSequences);
	Logger::Log.log(Logger::LogLevel::Always) << loopSequences.size() << " loops in ONTs" << std::endl;
	writeLoopGraphSequences(graph, params.basePath + "/tmp/loops.fa", loopSequences, pathStartClip, pathEndClip);
	Logger::Log.log(Logger::LogLevel::Always) << "cluster loops roughly" << std::endl;
	orderLoopsByLength(loopSequences, graph, pathStartClip, pathEndClip);
	Logger::Log.log(Logger::LogLevel::Always) << "max clustering edit distance " << params.maxClusterDifference << std::endl;
	auto clusters = roughClusterLoopSequences(loopSequences, graph, pathStartClip, pathEndClip, coreNodes, params.maxClusterDifference, params.numThreads);
	Logger::Log.log(Logger::LogLevel::Always) << clusters.size() << " rough clusters" << std::endl;
	Logger::Log.log(Logger::LogLevel::Always) << "getting exact locations of raw loop sequences" << std::endl;
	addRawSequencesToLoops(clusters, graph, pathStartClip, pathEndClip, params.ontReadPath, params.numThreads);
	if (params.extraPhasing)
	{
		Logger::Log.log(Logger::LogLevel::Always) << "phase clusters by raw sequences" << std::endl;
		clusters = editSplitClusters(clusters, graph, pathStartClip, pathEndClip, params.numThreads, params.MBGPath, params.basePath + "/tmp");
	}
	else
	{
		Logger::Log.log(Logger::LogLevel::Always) << "cluster loops by density" << std::endl;
		clusters = densityClusterLoops(clusters, graph, pathStartClip, pathEndClip, coreNodes, params.maxClusterDifference, 5, params.minReclusterDistance, params.numThreads);
	}
	Logger::Log.log(Logger::LogLevel::Always) << clusters.size() << " clusters" << std::endl;
	std::sort(clusters.begin(), clusters.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
	Logger::Log.log(Logger::LogLevel::Always) << "get possible missing ONT loop sequences" << std::endl;
	auto missingLoops = getMissingLoopSequences(clusters, graph, borderNodes, anchorNodes, pathStartClip, pathEndClip, ontPaths, minLength);
	Logger::Log.log(Logger::LogLevel::Always) << missingLoops.size() << " possible missing ONT loop sequences" << std::endl;
	Logger::Log.log(Logger::LogLevel::Always) << "write possible missing ONT loop sequences" << std::endl;
	writeMissingOntLoopSequences(params.basePath + "/missing_loops.fa", missingLoops);
	Logger::Log.log(Logger::LogLevel::Always) << "getting morph consensuses" << std::endl;
	auto morphConsensuses = getMorphConsensuses(clusters, graph, pathStartClip, pathEndClip);
//	Logger::Log.log(Logger::LogLevel::Always) << "polishing morph consensuses" << std::endl;
//	polishMorphConsensuses(morphConsensuses, params.MBGPath, params.basePath + "/tmp", params.numThreads);
	addMorphTypes(morphConsensuses, graph, borderNodes, anchorNodes);
	nameMorphConsensuses(morphConsensuses, params.namePrefix);
	Logger::Log.log(Logger::LogLevel::Always) << "write morph paths" << std::endl;
	writeMorphPaths(params.basePath + "/morphs.gaf", morphConsensuses, graph, pathStartClip, pathEndClip);
	return morphConsensuses;
}

void PostprocessONTClusters(std::vector<MorphConsensus>& morphConsensuses, const ClusterParams& params)
{
	Logger::Log.log(Logger::LogLevel::Always) << "polishing morph consensuses" << std::endl;
	polishMorphConsensuses(morphConsensuses, params.MBGPath, params.basePath + "/tmp", params.numThreads);
	addRawSequenceNamesToLoops(morphConsensuses);
	Logger::Log.log(Logger::LogLevel::Always) << "get self-corrected ONT loop sequences" << std::endl;
	addSelfCorrectedOntLoopSequences(morphConsensuses);
	nameMorphConsensuses(morphConsensuses, params.namePrefix);
}

void WriteONTClusters(const ClusterParams& params, const std::vector<MorphConsensus>& morphConsensuses)
{
	Logger::Log.log(Logger::LogLevel::Always) << "write raw ONT loop sequences" << std::endl;
	writeRawOntLoopSequences(params.basePath + "/raw_loops.fa", morphConsensuses);
	Logger::Log.log(Logger::LogLevel::Always) << "write self-corrected ONT loop sequences" << std::endl;
	writeSelfCorrectedOntLoopSequences(params.basePath + "/tmp/loops_selfcorrected.fa", morphConsensuses);
	Logger::Log.log(Logger::LogLevel::Always) << "write morph consensuses" << std::endl;
	writeMorphConsensuses(params.basePath + "/morphs.fa", params.basePath + "/tmp/morphs_preconsensus.fa", morphConsensuses);
	Logger::Log.log(Logger::LogLevel::Always) << "write morph graph and read paths" << std::endl;
	auto morphGraph = getMorphGraph(morphConsensuses);
	writeMorphGraphAndReadPaths(params.basePath + "/morphgraph.gfa", params.basePath + "/readpaths-morphgraph.gaf", morphGraph);
	if (params.winnowmapPath != "" && params.samtoolsPath != "")
	{
		Logger::Log.log(Logger::LogLevel::Always) << "realign raw ONT loop sequences to morph consensuses" << std::endl;
		alignRawOntLoopsToMorphConsensuses(params.basePath + "/raw_loop_to_morphs_alignments.bam", params.basePath + "/tmp", params.numThreads, morphConsensuses, params.winnowmapPath, params.samtoolsPath);
		Logger::Log.log(Logger::LogLevel::Always) << "realign self-corrected ONT loop sequences to morph consensuses" << std::endl;
		alignSelfCorrectedOntLoopsToMorphConsensuses(params.basePath + "/tmp/selfcorrected_loop_to_morphs_alignments.bam", params.basePath + "/tmp", params.numThreads, morphConsensuses, params.winnowmapPath, params.samtoolsPath);
	}
	else
	{
		if (params.winnowmapPath == "")
		{
			Logger::Log.log(Logger::LogLevel::Always) << "winnowmap not found, ";
		}
		if (params.samtoolsPath != "")
		{
			Logger::Log.log(Logger::LogLevel::Always) << "bamtools not found, ";
		}
		Logger::Log.log(Logger::LogLevel::Always) << "skipping alignment of raw ONT loop sequences to morph consensuses" << std::endl;
	}
	if (params.annotationFasta.size() > 0)
	{
		Logger::Log.log(Logger::LogLevel::Always) << "lifting over annotations to morphs" << std::endl;
		liftoverAnnotationsToMorphs(params.basePath, morphConsensuses, params.annotationFasta, params.annotationGff3, params.basePath + "/tmp", params.liftoffPath);
	}
}

void expandCircularComponent(std::vector<bool>& partOfCircularComponent, const GfaGraph& graph, const size_t startNode)
{
	phmap::flat_hash_set<Node> reachable;
	std::vector<Node> checkStack;
	checkStack.emplace_back(startNode, true);
	checkStack.emplace_back(startNode, false);
	while (checkStack.size() >= 1)
	{
		auto top = checkStack.back();
		checkStack.pop_back();
		if (reachable.count(top) == 1) continue;
		reachable.insert(top);
		if (graph.edges.count(top) == 1)
		{
			for (auto edge : graph.edges.at(top))
			{
				checkStack.emplace_back(std::get<0>(edge));
			}
		}
	}
	for (Node node : reachable)
	{
		if (reachable.count(reverse(node)) == 0) continue;
		partOfCircularComponent[node.id()] = true;
	}
}

phmap::flat_hash_set<size_t> getContigsInCircularComponent(const GfaGraph& graph)
{
	std::vector<bool> partOfCircularComponent;
	partOfCircularComponent.resize(graph.numNodes(), false);
	for (size_t i = 0; i < graph.numNodes(); i++)
	{
		if (partOfCircularComponent[i]) continue;
		if (getSelfDistance(graph, i) == std::numeric_limits<size_t>::max()) continue;
		partOfCircularComponent[i] = true;
		expandCircularComponent(partOfCircularComponent, graph, i);
	}
	phmap::flat_hash_set<size_t> result;
	for (size_t i = 0; i < partOfCircularComponent.size(); i++)
	{
		if (partOfCircularComponent[i]) result.emplace(i);
	}
	return result;
}

void createMergedProcessedGraph(const std::string& outputPrefix, const size_t numTangles, const std::string& resultGraphPath)
{
	std::ofstream file { resultGraphPath };
	size_t shortestBorderNode = std::numeric_limits<size_t>::max();
	size_t longestBorderOverlap = 0;
	for (size_t i = 0; i < numTangles; i++)
	{
		GfaGraph graph;
		graph.loadFromFile(outputPrefix + std::to_string(i) + "/processed-graph.gfa");
		phmap::flat_hash_set<size_t> skipBorders;
		for (size_t j = 0; j < graph.nodeNames.size(); j++)
		{
			if (graph.nodeNames[j].substr(0, 6) == "border")
			{
				skipBorders.insert(j);
				shortestBorderNode = std::min(shortestBorderNode, graph.nodeSeqs[j].size());
			}
			graph.nodeNames[j] = "tangle" + std::to_string(i) + "_" + graph.nodeNames[j];
		}
		for (size_t node = 0; node < graph.numNodes(); node++)
		{
			if (skipBorders.count(node) == 1) continue;
			file << "S\t" << graph.nodeNames[node] << "\t" << graph.nodeSeqs[node] << "\tll:f:" << graph.nodeCoverages[node] << "\tFC:i:" << (graph.nodeCoverages[node] * graph.nodeSeqs[node].size()) << std::endl;
		}
		for (const auto& pair : graph.edges)
		{
			for (const auto& target : pair.second)
			{
				if (skipBorders.count(pair.first.id()) == 1 || skipBorders.count(std::get<0>(target).id()) == 1)
				{
					longestBorderOverlap = std::max(longestBorderOverlap, std::get<1>(target));
				}
				file << "L\t" << graph.nodeNames[pair.first.id()] << "\t" << (pair.first.forward() ? "+" : "-") << "\t" << graph.nodeNames[std::get<0>(target).id()] << "\t" << (std::get<0>(target).forward() ? "+" : "-") << "\t" << std::get<1>(target) << "M\tec:i:" << std::get<2>(target) << std::endl;
			}
		}
	}
	size_t clipToLength = shortestBorderNode;
	if (longestBorderOverlap >= clipToLength) clipToLength = longestBorderOverlap+1;
	for (size_t i = 0; i < numTangles; i++)
	{
		GfaGraph graph;
		graph.loadFromFile(outputPrefix + std::to_string(i) + "/processed-graph.gfa");
		phmap::flat_hash_set<size_t> keepBorders;
		for (size_t j = 0; j < graph.nodeNames.size(); j++)
		{
			if (graph.nodeNames[j].substr(0, 6) == "border")
			{
				keepBorders.insert(j);
			}
			graph.nodeNames[j] = "tangle" + std::to_string(i) + "_" + graph.nodeNames[j];
		}
		for (auto node : keepBorders)
		{
			size_t keepAmount = std::min(graph.nodeSeqs[node].size(), clipToLength);
			file << "S\t" << graph.nodeNames[node] << "\t" << graph.nodeSeqs[node].substr(0, keepAmount) << "\tll:f:" << graph.nodeCoverages[node] << "\tFC:i:" << (graph.nodeCoverages[node] * graph.nodeSeqs[node].size()) << std::endl;
		}
	}
}

size_t parseTangleNumber(const std::string& nodeName)
{
	assert(nodeName.size() >= 8);
	assert(nodeName.substr(0, 6) == "tangle");
	size_t firstUnderscore = 7;
	while (nodeName[firstUnderscore] != '_')
	{
		firstUnderscore += 1;
		assert(firstUnderscore < nodeName.size());
	}
	size_t result = std::stoull(nodeName.substr(6, firstUnderscore-6));
	return result;
}

void splitONTReadsPerTangle(const std::string& ontBasePath, const std::string& outputPrefix, const size_t numTangles, const std::string& mergedProcessedGraphPath, const std::string& ontTangleAlignmentPath, const std::string& ontAlignmentPath)
{
	phmap::flat_hash_map<std::string, phmap::flat_hash_set<size_t>> bestTanglesPerRead;
	phmap::flat_hash_map<std::string, phmap::flat_hash_set<size_t>> secondBestTanglesPerRead;
	phmap::flat_hash_map<std::string, phmap::flat_hash_set<size_t>> anyTanglesPerRead;
	phmap::flat_hash_set<std::string> allReads;
	{
		GfaGraph graph;
		graph.loadFromFile(mergedProcessedGraphPath);
		std::vector<ReadPath> paths = loadReadPaths(ontTangleAlignmentPath, graph);
		for (const auto& path : paths)
		{
			allReads.insert(path.readName);
			size_t tangleNum = parseTangleNumber(graph.nodeNames[path.path[0].id()]);
			assert(tangleNum < numTangles);
			if (path.readStart < 1000 && path.readEnd + 1000 > path.readLength && path.mapq >= 20)
			{
				bestTanglesPerRead[path.readName].insert(tangleNum);
			}
			else if (path.mapq >= 20)
			{
				secondBestTanglesPerRead[path.readName].insert(tangleNum);
			}
			else
			{
				anyTanglesPerRead[path.readName].insert(tangleNum);
			}
		}
	}
	phmap::flat_hash_map<std::string, size_t> readAssignment;
	for (const std::string& read : allReads)
	{
		if (bestTanglesPerRead.count(read) == 1)
		{
			if (bestTanglesPerRead.at(read).size() == 1)
			{
				readAssignment[read] = *bestTanglesPerRead.at(read).begin();
			}
			else
			{
				readAssignment[read] = std::numeric_limits<size_t>::max();
			}
		}
		else if (secondBestTanglesPerRead.count(read) == 1)
		{
			if (secondBestTanglesPerRead.at(read).size() == 1)
			{
				readAssignment[read] = *secondBestTanglesPerRead.at(read).begin();
			}
			else
			{
				readAssignment[read] = std::numeric_limits<size_t>::max();
			}
		}
		else if (anyTanglesPerRead.count(read) == 1)
		{
			if (anyTanglesPerRead.at(read).size() == 1)
			{
				readAssignment[read] = *anyTanglesPerRead.at(read).begin();
			}
			else
			{
				readAssignment[read] = std::numeric_limits<size_t>::max();
			}
		}
	}
	{
		std::vector<std::unique_ptr<std::ofstream>> readFilesPerTangle;
		std::vector<std::unique_ptr<std::ofstream>> alnFilesPerTangle;
		for (size_t i = 0; i < numTangles; i++)
		{
			readFilesPerTangle.emplace_back(new std::ofstream { outputPrefix + std::to_string(i) + "/ont_reads.fa" });
			alnFilesPerTangle.emplace_back(new std::ofstream { outputPrefix + std::to_string(i) + "/ont-alns.gaf" });
		}
		FastQ::streamFastqFromFile(ontBasePath, false, [&readFilesPerTangle, &readAssignment](const FastQ& read)
		{
			if (readAssignment.count(read.seq_id) == 0) return;
			if (readAssignment.at(read.seq_id) == std::numeric_limits<size_t>::max()) return;
			size_t tangle = readAssignment.at(read.seq_id);
			assert(tangle < readFilesPerTangle.size());
			(*readFilesPerTangle[tangle]) << ">" << read.seq_id << std::endl;
			(*readFilesPerTangle[tangle]) << read.sequence << std::endl;
		});
		{
			std::ifstream file { ontAlignmentPath };
			while (file.good())
			{
				std::string line;
				getline(file, line);
				if (!file.good()) break;
				auto parts = split(line, '\t');
				std::string readname = parts[0];
				if (readAssignment.count(readname) == 0) continue;
				if (readAssignment.at(readname) == std::numeric_limits<size_t>::max()) continue;
				size_t tangle = readAssignment.at(readname);
				assert(tangle < alnFilesPerTangle.size());
				(*alnFilesPerTangle[tangle]) << line << std::endl;
			}
		}
	}
}

void copyHifiFiles(const std::string& fromPath, const std::string& toPath)
{
	std::string command = "cp " + fromPath + "/processed-graph.gfa " + toPath + "/processed-graph.gfa";
	system(command.c_str());
	command = "cp " + fromPath + "/consensus_path.gaf " + toPath + "/consensus_path.gaf";
	system(command.c_str());
}

void copyTangleHifiFiles(const std::string& basePath)
{
	std::string command = "cp " + basePath + "/processed-graph.gfa " + basePath + "/tangle-processed-graph.gfa";
	system(command.c_str());
	command = "cp " + basePath + "/consensus_path.gaf " + basePath + "/tangle_consensus_path.gaf";
	system(command.c_str());
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
		std::ifstream variantgraph { outputPrefix + std::to_string(i) + "/processed-graph.gfa" };
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

void runMultipleTangles(const ClusterParams& baseParams, const std::string& outputPrefix, const size_t numTangles, const std::vector<std::vector<std::string>>& hifiReadsPerTangle, const std::vector<std::string>& hifiReadFiles, const bool doUL, const std::string& ulTmpFolder, const std::vector<std::string>& ontReadFiles)
{
	std::vector<std::string> readFileNames;
	for (size_t i = 0; i < numTangles; i++)
	{
		readFileNames.push_back(outputPrefix + std::to_string(i) + "/hifi_reads.fa");
	}
	std::cerr << "extracting HiFi/duplex reads per tangle" << std::endl;
	bool allReadsFound = splitReads(hifiReadFiles, hifiReadsPerTangle, readFileNames);
	if (!allReadsFound)
	{
		std::cerr << "WARNING: some HiFi reads were not found in the input read files. Double check that the HiFi read files have not been deleted or moved after running verkko." << std::endl;
	}
	std::vector<size_t> tanglesWithoutReads;
	std::vector<std::string> tangleBasePaths;
	for (size_t i = 0; i < numTangles; i++)
	{
		if (hifiReadsPerTangle[i].size() == 0)
		{
			std::cerr << "WARNING: tangle " << i << " has no HiFi/duplex reads, skipping" << std::endl;
			tanglesWithoutReads.push_back(i);
			continue;
		}
		ClusterParams clusterParams = baseParams;
		clusterParams.namePrefix = clusterParams.namePrefix + (clusterParams.namePrefix.size() > 0 ? "_" : "") + "tangle" + std::to_string(i);
		clusterParams.basePath = outputPrefix + std::to_string(i);
		clusterParams.hifiReadPath = outputPrefix + std::to_string(i) + "/hifi_reads.fa";
		std::cerr << "running tangle " << i << " in folder " << outputPrefix + std::to_string(i) << std::endl;
		HandleCluster(clusterParams);
		copyTangleHifiFiles(clusterParams.basePath);
	}
	if (doUL)
	{
		std::string selectedONTPath = ulTmpFolder + "/ont_reads.fa";
		ClusterParams clusterParams = baseParams;
		clusterParams.basePath = ulTmpFolder;
		clusterParams.hifiReadPath = clusterParams.basePath + "/hifi_reads.fa";
		std::filesystem::create_directories(ulTmpFolder);
		std::filesystem::create_directories(ulTmpFolder + "/tmp");
		std::cerr << "getting kmers from tangles" << std::endl;
		std::string tmpKmerFile = ulTmpFolder + "/rdna_kmers.fa";
		getKmers(outputPrefix, numTangles, tmpKmerFile);
		std::cerr << "extracting HiFi/duplex reads" << std::endl;
		{
			std::ofstream readsfile { clusterParams.basePath + "/hifi_reads.fa" };
			iterateMatchingReads(tmpKmerFile, hifiReadFiles, 101, 2000, [&readsfile](const FastQ& seq)
			{
				readsfile << ">" << seq.seq_id << std::endl;
				readsfile << seq.sequence << std::endl;
			});
		}
		std::cerr << "running" << std::endl;
		HandleCluster(clusterParams);
		std::cerr << "extracting ultralong ONT reads" << std::endl;
		size_t consensusLength = medianConsensusLength(outputPrefix, numTangles);
		std::cerr << "median consensus length " << consensusLength << ", using " << consensusLength/2 << " as minimum ONT match length" << std::endl;
		{
			std::ofstream readsfile { selectedONTPath };
			iterateMatchingReads(ulTmpFolder + "/rdna_kmers.fa", ontReadFiles, 21, consensusLength/2, [&readsfile](const FastQ& seq)
			{
				readsfile << ">" << seq.seq_id << std::endl;
				readsfile << seq.sequence << std::endl;
			});
		}
		clusterParams.ontReadPath = selectedONTPath;
		std::cerr << "creating merged tanglegraph" << std::endl;
		createMergedProcessedGraph(outputPrefix, numTangles, clusterParams.basePath + "/merged_processed_graph.gfa");
		std::cerr << "aligning ONT reads to merged tanglegraph" << std::endl;
		AlignONTReads(clusterParams.basePath, clusterParams.GraphAlignerPath, clusterParams.ontReadPath, clusterParams.basePath + "/merged_processed_graph.gfa", clusterParams.basePath + "/ont-alns-pertangle.gaf", clusterParams.numThreads);
		std::cerr << "aligning ONT reads to k-mer graph" << std::endl;
		AlignONTReads(clusterParams.basePath, clusterParams.GraphAlignerPath, clusterParams.ontReadPath, clusterParams.basePath + "/processed-graph.gfa", clusterParams.basePath + "/ont-alns.gaf", clusterParams.numThreads);
		std::cerr << "split ONT reads by tangle" << std::endl;
		splitONTReadsPerTangle(clusterParams.ontReadPath, outputPrefix, numTangles, clusterParams.basePath + "/merged_processed_graph.gfa", clusterParams.basePath + "/ont-alns-pertangle.gaf", clusterParams.basePath + "/ont-alns.gaf");
		for (size_t i = 0; i < numTangles; i++)
		{
			std::cerr << "running tangle " << i << std::endl;
			ClusterParams tangleClusterParams = clusterParams;
			tangleClusterParams.namePrefix = tangleClusterParams.namePrefix + (tangleClusterParams.namePrefix.size() > 0 ? "_" : "") + "tangle" + std::to_string(i);
			tangleClusterParams.basePath = outputPrefix + std::to_string(i);
			copyHifiFiles(clusterParams.basePath, tangleClusterParams.basePath);
			auto clusters = GetONTClusters(tangleClusterParams);
			PostprocessONTClusters(clusters, tangleClusterParams);
			WriteONTClusters(tangleClusterParams, clusters);
		}
	}
	if (tanglesWithoutReads.size() > 0)
	{
		std::cerr << "WARNING: some tangles did not have any HiFi/duplex reads assigned:";
		for (auto tangle : tanglesWithoutReads) std::cerr << " " << tangle;
		std::cerr << ", something likely went wrong." << std::endl;
	}
}
