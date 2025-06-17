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

void processGraphAndWrite(const Path& heavyPath, const GfaGraph& graph, const std::string outputGraphName)
{
	const size_t minCoverage = 5;
	const size_t minAnchorLength = 10000; // this should be a bit larger than anchorDistanceFromMainLoop in LoopAnalysis.cpp:getBorderNodes
	phmap::flat_hash_set<size_t> keptNodesInCycles = getNodesInMainCycle(graph, heavyPath, minCoverage);
	phmap::flat_hash_map<Node, size_t> shortestDistanceToCyclicComponent = getNodeDistancesToMainCycle(graph, keptNodesInCycles);
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
	phmap::flat_hash_map<Node, size_t> nodeTrim;
	for (Node node : nodesWithEnoughDistance)
	{
		assert(shortestDistanceToCyclicComponent.count(node) == 0);
		assert(shortestDistanceToCyclicComponent.count(reverse(node)) == 1);
		assert(shortestDistanceToCyclicComponent.at(reverse(node)) >= minAnchorLength);
		if (shortestDistanceToCyclicComponent.at(reverse(node)) == minAnchorLength) continue;
		nodeTrim[node] = shortestDistanceToCyclicComponent.at(reverse(node)) - minAnchorLength;
		if (nodeTrim[node] >= graph.nodeSeqs[node.id()].size())
		{
			if (keptNodesInAnchors.count(node.id()) == 1)
			{
				keptNodesInAnchors.erase(node.id());
			}
		}
	}
	phmap::flat_hash_set<size_t> finalKeptNodes;
	finalKeptNodes.insert(keptNodesInCycles.begin(), keptNodesInCycles.end());
	finalKeptNodes.insert(keptNodesInAnchors.begin(), keptNodesInAnchors.end());
	std::ofstream file { outputGraphName };
	for (size_t node : finalKeptNodes)
	{
		std::string sequence = graph.nodeSeqs[node];
		assert(nodeTrim.count(Node { node, true }) == 0 || nodeTrim.count(Node { node, false }) == 0);
		if (nodeTrim.count(Node { node, true }) == 1)
		{
			assert(nodeTrim.at(Node { node, true }) < sequence.size());
			sequence = sequence.substr(nodeTrim.at(Node { node, true }));
		}
		if (nodeTrim.count(Node { node, false }) == 1)
		{
			assert(nodeTrim.at(Node { node, false }) < sequence.size());
			sequence = sequence.substr(0, sequence.size()-nodeTrim.at(Node { node, false }));
		}
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
	processGraphAndWrite(heavyPath, graph, params.basePath + "/processed-graph.gfa");
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
	Logger::Log.log(Logger::LogLevel::Always) << "get possible missing ONT loop sequences" << std::endl;
	auto missingLoops = getMissingLoopSequences(clusters, graph, borderNodes, anchorNodes, pathStartClip, pathEndClip, ontPaths, minLength);
	std::cerr << missingLoops.size() << " possible missing ONT loop sequences" << std::endl;
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

std::vector<std::vector<MorphConsensus>> splitClustersByTangle(const std::vector<MorphConsensus>& morphConsensuses, const std::string outputPrefix, const size_t numTangles)
{
	assert(numTangles >= 1);
	if (numTangles == 1)
	{
		std::vector<std::vector<MorphConsensus>> result;
		result.emplace_back(morphConsensuses);
		return result;
	}
	assert(numTangles >= 2);
	const size_t k = 101;
	std::vector<KmerMatcher> perTangleMatchers;
	for (size_t tangle = 0; tangle < numTangles; tangle++)
	{
		perTangleMatchers.emplace_back(k, k);
		GfaGraph graph;
		graph.loadFromFile(outputPrefix + std::to_string(tangle) + "/processed-graph.gfa");
		phmap::flat_hash_set<size_t> contigsInCircularComponent = getContigsInCircularComponent(graph);
		for (size_t node : contigsInCircularComponent)
		{
			perTangleMatchers[tangle].addReferenceKmers(graph.nodeSeqs[node]);
			auto rc = revcomp(graph.nodeSeqs[node]);
			perTangleMatchers[tangle].addReferenceKmers(rc);
		}
	}
	for (size_t tangle = 0; tangle < numTangles; tangle++)
	{
		GfaGraph graph;
		graph.loadFromFile(outputPrefix + std::to_string(tangle) + "/processed-graph.gfa");
		for (size_t node = 0; node < graph.numNodes(); node++)
		{
			for (size_t otherTangle = 0; otherTangle < numTangles; otherTangle++)
			{
				if (otherTangle == tangle) continue;
				perTangleMatchers[otherTangle].removeReferenceKmers(graph.nodeSeqs[node]);
				auto rc = revcomp(graph.nodeSeqs[node]);
				perTangleMatchers[otherTangle].removeReferenceKmers(rc);
			}
		}
	}
	std::vector<std::vector<size_t>> matchCount;
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		matchCount.emplace_back();
		matchCount.back().resize(numTangles);
		for (size_t j = 0; j < numTangles; j++)
		{
			matchCount.back()[j] = perTangleMatchers[j].getMatchKmerCount(morphConsensuses[i].sequence);
		}
	}
	MorphGraph morphgraph = getMorphGraph(morphConsensuses);
	std::vector<size_t> morphClusterAssignment;
	morphClusterAssignment.resize(morphConsensuses.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		if (morphConsensuses[i].type != MorphConsensus::MorphType::Inner) continue;
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "morph " << i << " tangle match counts";
		for (size_t j = 0; j < numTangles; j++)
		{
			Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << " " << matchCount[i][j];
		}
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << std::endl;
		size_t maxTangle = std::numeric_limits<size_t>::max();
		size_t maxTangleMatches = 0;
		for (size_t j = 0; j < numTangles; j++)
		{
			if (matchCount[i][j] == maxTangleMatches)
			{
				maxTangle = std::numeric_limits<size_t>::max();
			}
			if (matchCount[i][j] > maxTangleMatches)
			{
				maxTangle = j;
				maxTangleMatches = matchCount[i][j];
			}
		}
		morphClusterAssignment[i] = maxTangle;
	}
	std::vector<size_t> morphUniqueTangleFw;
	std::vector<size_t> morphUniqueTangleBw;
	morphUniqueTangleBw.resize(morphConsensuses.size(), std::numeric_limits<size_t>::max());
	morphUniqueTangleFw.resize(morphConsensuses.size(), std::numeric_limits<size_t>::max());
	for (auto pair : morphgraph.edgeCoverage)
	{
		Node from = pair.first.first;
		Node to = pair.first.second;
		if (morphClusterAssignment[to.id()] != std::numeric_limits<size_t>::max())
		{
			if (from.forward())
			{
				if (morphUniqueTangleFw[from.id()] == std::numeric_limits<size_t>::max())
				{
					morphUniqueTangleFw[from.id()] = morphClusterAssignment[to.id()];
				}
				else if (morphUniqueTangleFw[from.id()] != morphClusterAssignment[to.id()])
				{
					morphUniqueTangleFw[from.id()] = std::numeric_limits<size_t>::max()-1;
				}
			}
			else
			{
				if (morphUniqueTangleBw[from.id()] == std::numeric_limits<size_t>::max())
				{
					morphUniqueTangleBw[from.id()] = morphClusterAssignment[to.id()];
				}
				else if (morphUniqueTangleBw[from.id()] != morphClusterAssignment[to.id()])
				{
					morphUniqueTangleBw[from.id()] = std::numeric_limits<size_t>::max()-1;
				}
			}
		}
		if (morphClusterAssignment[from.id()] != std::numeric_limits<size_t>::max())
		{
			if (to.forward())
			{
				if (morphUniqueTangleFw[to.id()] == std::numeric_limits<size_t>::max())
				{
					morphUniqueTangleFw[to.id()] = morphClusterAssignment[from.id()];
				}
				else if (morphUniqueTangleFw[to.id()] != morphClusterAssignment[from.id()])
				{
					morphUniqueTangleFw[to.id()] = std::numeric_limits<size_t>::max()-1;
				}
			}
			else
			{
				if (morphUniqueTangleBw[to.id()] == std::numeric_limits<size_t>::max())
				{
					morphUniqueTangleBw[to.id()] = morphClusterAssignment[from.id()];
				}
				else if (morphUniqueTangleBw[to.id()] != morphClusterAssignment[from.id()])
				{
					morphUniqueTangleBw[to.id()] = std::numeric_limits<size_t>::max()-1;
				}
			}
		}
	}
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		if (morphUniqueTangleBw[i] != morphUniqueTangleFw[i]) continue;
		if (morphUniqueTangleBw[i] == morphClusterAssignment[i]) continue;
		if (morphConsensuses[i].type == MorphConsensus::MorphType::Inner)
		{
			if (morphClusterAssignment[i] != std::numeric_limits<size_t>::max())
			{
				if (matchCount[i][morphUniqueTangleBw[i]] < matchCount[i][morphClusterAssignment[i]] * 0.5)
				{
					continue;
				}
			}
		}
		morphClusterAssignment[i] = morphUniqueTangleBw[i];
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "topology based reassign morph " << i << " to tangle " << morphClusterAssignment[i] << std::endl;
	}
	phmap::flat_hash_map<std::string, size_t> readTangleAssignment;
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		if (morphClusterAssignment[i] == std::numeric_limits<size_t>::max()) continue;
		for (size_t j = 0; j < morphConsensuses[i].ontLoops.size(); j++)
		{
			if (readTangleAssignment.count(morphConsensuses[i].ontLoops[j].readName) == 0)
			{
				readTangleAssignment[morphConsensuses[i].ontLoops[j].readName] = morphClusterAssignment[i];
			}
			else if (readTangleAssignment.at(morphConsensuses[i].ontLoops[j].readName) != morphClusterAssignment[i])
			{
				readTangleAssignment[morphConsensuses[i].ontLoops[j].readName] = std::numeric_limits<size_t>::max();
			}
		}
	}
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		if (morphClusterAssignment[i] != std::numeric_limits<size_t>::max()) continue;
		std::vector<size_t> readMatchesPerCluster;
		readMatchesPerCluster.resize(numTangles, 0);
		for (size_t j = 0; j < morphConsensuses[i].ontLoops.size(); j++)
		{
			if (readTangleAssignment.count(morphConsensuses[i].ontLoops[j].readName) == 0) continue;
			if (readTangleAssignment.at(morphConsensuses[i].ontLoops[j].readName) == std::numeric_limits<size_t>::max()) continue;
			readMatchesPerCluster[readTangleAssignment.at(morphConsensuses[i].ontLoops[j].readName)] += 1;
		}
		size_t maxTangle = std::numeric_limits<size_t>::max();
		size_t maxTangleMatches = 0;
		size_t countWithMatches = 0;
		for (size_t j = 0; j < numTangles; j++)
		{
			if (readMatchesPerCluster[j] > 0) countWithMatches += 1;
			if (readMatchesPerCluster[j] == maxTangleMatches) maxTangle = std::numeric_limits<size_t>::max();
			if (readMatchesPerCluster[j] > maxTangleMatches)
			{
				maxTangleMatches = readMatchesPerCluster[j];
				maxTangle = j;
			}
		}
		morphClusterAssignment[i] = maxTangle;
		if (countWithMatches >= 2 && morphConsensuses[i].type != MorphConsensus::MorphType::Inner)
		{
			morphClusterAssignment[i] = std::numeric_limits<size_t>::max()-1;
		}
		Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "read based reassign morph " << i << " to tangle " << morphClusterAssignment[i] << std::endl;
	}
	std::vector<std::vector<MorphConsensus>> result;
	result.resize(numTangles+1);
	for (size_t i = 0; i < morphConsensuses.size(); i++)
	{
		if (morphClusterAssignment[i] == std::numeric_limits<size_t>::max())
		{
			result.back().emplace_back(morphConsensuses[i]);
			continue;
		}
		if (morphClusterAssignment[i] == std::numeric_limits<size_t>::max()-1)
		{
			std::vector<std::vector<size_t>> readsPerTangle;
			readsPerTangle.resize(numTangles);
			for (size_t j = 0; j < morphConsensuses[i].ontLoops.size(); j++)
			{
				if (readTangleAssignment.count(morphConsensuses[i].ontLoops[j].readName) == 0) continue;
				if (readTangleAssignment.at(morphConsensuses[i].ontLoops[j].readName) == std::numeric_limits<size_t>::max()) continue;
				readsPerTangle[readTangleAssignment.at(morphConsensuses[i].ontLoops[j].readName)].emplace_back(j);
			}
			for (size_t j = 0; j < numTangles; j++)
			{
				if (readsPerTangle[j].size() == 0) continue;
				Logger::Log.log(Logger::LogLevel::DetailedDebugInfo) << "split morph " << i << " to tangle " << j << " read count " << readsPerTangle[j].size() << std::endl;
				result[j].emplace_back();
				result[j].back() = morphConsensuses[i];
				result[j].back().ontLoops.clear();
				for (size_t k : readsPerTangle[j])
				{
					result[j].back().ontLoops.emplace_back(morphConsensuses[i].ontLoops[k]);
				}
				result[j].back().coverage = result[j].back().ontLoops.size();
			}
			continue;
		}
		result[morphClusterAssignment[i]].emplace_back(morphConsensuses[i]);
	}
	{
		std::ofstream file { "tangle_match_counts.csv" };
		file << "node\tmatch_counts\tassignment" << std::endl;
		for (size_t i = 0; i < morphConsensuses.size(); i++)
		{
			file << morphConsensuses[i].name << "\t";
			for (size_t j = 0; j < matchCount[i].size(); j++)
			{
				if (j != 0) file << "_";
				file << matchCount[i][j];
			}
			file << "\t" << morphClusterAssignment[i] << std::endl;
		}
	}
	return result;
}
