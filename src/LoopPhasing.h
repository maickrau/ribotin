#ifndef LoopPhasing_h
#define LoopPhasing_h

#include <cstddef>
#include <string>
#include <vector>
#include "ClusterMisc.h"

std::vector<std::vector<OntLoop>> editSplitClusters(const std::vector<std::vector<OntLoop>>& clusters, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t numThreads, const std::string MBGPath, const std::string tmpPath);

#endif
