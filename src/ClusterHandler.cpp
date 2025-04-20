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
	processGraphAndWrite(heavyPath, graph, params.basePath + "/processed-graph.gfa");
	if (params.annotationFasta.size() > 0)
	{
		Logger::Log.log(Logger::LogLevel::Always) << "lifting over annotations to consensus" << std::endl;
		liftoverAnnotationsToConsensus(params.basePath, params.basePath + "/consensus.fa", params.annotationFasta, params.annotationGff3, params.basePath + "/tmp", params.liftoffPath);
	}
}

void DoClusterONTAnalysis(const ClusterParams& params)
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
	auto clusters = roughClusterLoopSequences(loopSequences, graph, pathStartClip, pathEndClip, coreNodes, params.maxClusterDifference);
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
		clusters = densityClusterLoops(clusters, graph, pathStartClip, pathEndClip, coreNodes, params.maxClusterDifference, 5, params.minReclusterDistance);
	}
	Logger::Log.log(Logger::LogLevel::Always) << clusters.size() << " clusters" << std::endl;
	std::sort(clusters.begin(), clusters.end(), [](const auto& left, const auto& right) { return left.size() > right.size(); });
	addRawSequenceNamesToLoops(clusters);
	Logger::Log.log(Logger::LogLevel::Always) << "write raw ONT loop sequences" << std::endl;
	writeRawOntLoopSequences(params.basePath + "/raw_loops.fa", clusters);
	Logger::Log.log(Logger::LogLevel::Always) << "get possible missing ONT loop sequences" << std::endl;
	auto missingLoops = getMissingLoopSequences(clusters, graph, borderNodes, anchorNodes, pathStartClip, pathEndClip, ontPaths, minLength);
	std::cerr << missingLoops.size() << " possible missing ONT loop sequences" << std::endl;
	Logger::Log.log(Logger::LogLevel::Always) << "write possible missing ONT loop sequences" << std::endl;
	writeMissingOntLoopSequences(params.basePath + "/missing_loops.fa", missingLoops);
	Logger::Log.log(Logger::LogLevel::Always) << "get self-corrected ONT loop sequences" << std::endl;
	addSelfCorrectedOntLoopSequences(clusters);
	Logger::Log.log(Logger::LogLevel::Always) << "write self-corrected ONT loop sequences" << std::endl;
	writeSelfCorrectedOntLoopSequences(params.basePath + "/tmp/loops_selfcorrected.fa", clusters);
	Logger::Log.log(Logger::LogLevel::Always) << "getting morph consensuses" << std::endl;
	auto morphConsensuses = getMorphConsensuses(clusters, graph, pathStartClip, pathEndClip);
	Logger::Log.log(Logger::LogLevel::Always) << "polishing morph consensuses" << std::endl;
	polishMorphConsensuses(morphConsensuses, clusters, params.MBGPath, params.basePath + "/tmp", params.numThreads);
	Logger::Log.log(Logger::LogLevel::Always) << "write morph consensuses" << std::endl;
	nameMorphConsensuses(morphConsensuses, graph, borderNodes, anchorNodes, params.namePrefix);
	writeMorphConsensuses(params.basePath + "/morphs.fa", params.basePath + "/tmp/morphs_preconsensus.fa", morphConsensuses);
	Logger::Log.log(Logger::LogLevel::Always) << "write morph paths" << std::endl;
	writeMorphPaths(params.basePath + "/morphs.gaf", morphConsensuses, graph, pathStartClip, pathEndClip);
	Logger::Log.log(Logger::LogLevel::Always) << "write morph graph and read paths" << std::endl;
	writeMorphGraphAndReadPaths(params.basePath + "/morphgraph.gfa", params.basePath + "/readpaths-morphgraph.gaf", morphConsensuses);
	if (params.winnowmapPath != "" && params.samtoolsPath != "")
	{
		Logger::Log.log(Logger::LogLevel::Always) << "realign raw ONT loop sequences to morph consensuses" << std::endl;
		alignRawOntLoopsToMorphConsensuses(clusters, params.basePath + "/raw_loop_to_morphs_alignments.bam", params.basePath + "/tmp", params.numThreads, morphConsensuses, params.winnowmapPath, params.samtoolsPath);
		Logger::Log.log(Logger::LogLevel::Always) << "realign self-corrected ONT loop sequences to morph consensuses" << std::endl;
		alignSelfCorrectedOntLoopsToMorphConsensuses(clusters, params.basePath + "/tmp/selfcorrected_loop_to_morphs_alignments.bam", params.basePath + "/tmp", params.numThreads, morphConsensuses, params.winnowmapPath, params.samtoolsPath);
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
