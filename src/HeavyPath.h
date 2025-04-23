#ifndef HeavyPath_h
#define HeavyPath_h

#include <cstddef>
#include <vector>
#include "ClusterMisc.h"

size_t getSelfDistance(const GfaGraph& graph, const size_t startNode);
Path getHeavyPath(const GfaGraph& graph, const std::vector<ReadPath>& readPaths, const size_t localResolveLength);

#endif
