#ifndef LoopAnalysis_h
#define LoopAnalysis_h

#include <vector>
#include <string>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include "ClusterMisc.h"

class MorphGraph
{
public:
	std::vector<ReadPath> readPaths;
	std::vector<size_t> pathLength;
	std::map<std::pair<Node, Node>, size_t> edgeCoverage;
	std::vector<MorphConsensus> morphConsensuses;
};

std::tuple<std::unordered_set<size_t>, std::unordered_set<size_t>, std::unordered_map<Node, size_t>, std::unordered_map<Node, size_t>> getBorderNodes(const Path& heavyPath, const GfaGraph& graph);
std::vector<OntLoop> extractLoopSequences(const std::vector<ReadPath>& correctedPaths, const Path& heavyPath, const size_t minLength, const GfaGraph& graph, const std::unordered_set<size_t>& borderNodes, const std::unordered_set<size_t>& anchorNodes, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip);
std::vector<std::vector<OntLoop>> roughClusterLoopSequences(const std::vector<OntLoop>& loops, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::unordered_set<size_t>& coreNodes, const size_t maxEdits, const size_t numThreads);
std::vector<std::vector<OntLoop>> densityClusterLoops(const std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::unordered_set<size_t>& coreNodes, const size_t maxEdits, const size_t minPoints, const size_t minEdits, const size_t numThreads);
std::vector<MorphConsensus> getMorphConsensuses(const std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip);
void polishMorphConsensus(MorphConsensus& morphConsensus, const std::string MBGPath, const std::string tmpPath, const size_t numThreads);
void polishMorphConsensuses(std::vector<MorphConsensus>& morphConsensuses, const std::string MBGPath, const std::string tmpPath, const size_t numThreads);
void writeMorphConsensuses(std::string outFile, std::string preconsensusOutputFile, const std::vector<MorphConsensus>& morphConsensuses);
void writeMorphPaths(const std::string& outputFile, const std::vector<MorphConsensus>& morphConsensuses, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip);
void addSelfCorrectedOntLoopSequences(std::vector<MorphConsensus>& clusters);
void addRawSequenceNamesToLoops(std::vector<MorphConsensus>& clusters);
void addRawSequencesToLoops(std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::string ontReadPath, const size_t numThreads);
void alignRawOntLoopsToMorphConsensuses(std::string outputFile, std::string tmppath, const size_t numThreads, const std::vector<MorphConsensus>& morphConsensuses, const std::string& winnowmapPath, const std::string& samtoolsPath);
void alignSelfCorrectedOntLoopsToMorphConsensuses(std::string outputFile, std::string tmppath, const size_t numThreads, const std::vector<MorphConsensus>& morphConsensuses, const std::string& winnowmapPath, const std::string& samtoolsPath);
void writeRawOntLoopSequences(const std::string outputFile, const std::vector<MorphConsensus>& loopSequences);
void writeSelfCorrectedOntLoopSequences(const std::string outputFile, const std::vector<MorphConsensus>& loopSequences);
void orderLoopsByLength(std::vector<OntLoop>& loops, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip);
void liftoverAnnotationsToMorphs(const std::string& basePath, const std::vector<MorphConsensus>& morphConsensuses, const std::string& annotationFasta, const std::string& annotationGff3, const std::string& tmpPath, const std::string& liftoffPath);
void writeLoopGraphSequences(const GfaGraph& graph, const std::string& outputFile, const std::vector<OntLoop>& loopSequences, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip);
void addMorphTypes(std::vector<MorphConsensus>& morphConsensuses, const GfaGraph& graph, const std::unordered_set<size_t>& borderNodes, const std::unordered_set<size_t>& anchorNodes);
void nameMorphConsensuses(std::vector<MorphConsensus>& morphConsensuses, const std::string& namePrefix);
std::vector<OntLoop> getMissingLoopSequences(const std::vector<std::vector<OntLoop>>& presentClusters, const GfaGraph& graph, const std::unordered_set<size_t>& borderNodes, const std::unordered_set<size_t>& anchorNodes, const std::unordered_map<Node, size_t>&  pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const std::vector<ReadPath>& correctedPaths, const size_t minLength);
void writeMissingOntLoopSequences(const std::string outputFile, const std::vector<OntLoop>& loopSequences);
MorphGraph getMorphGraph(const std::vector<MorphConsensus>& morphConsensuses);
void writeMorphGraphAndReadPaths(const std::string& graphFile, const std::string& pathsFile, const MorphGraph& morphGraph);

#endif
