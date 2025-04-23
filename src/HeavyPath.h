#ifndef HeavyPath_h
#define HeavyPath_h

#include <cstddef>
#include <vector>
#include "ClusterMisc.h"

Path getHeavyPath(const GfaGraph& graph, const std::vector<ReadPath>& readPaths, const size_t localResolveLength);
size_t getShortestPathDistance(const GfaGraph& graph, const Node startNode, const Node endNode);
size_t getSelfDistance(const GfaGraph& graph, const size_t startNode);

#endif
