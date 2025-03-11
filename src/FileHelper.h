#ifndef FileHelper_h
#define FileHelper_h

#include <string>
#include <vector>
#include "ClusterMisc.h"

std::string getPathGaf(const Path& path, const GfaGraph& graph);
Path readHeavyPath(const GfaGraph& graph, const std::string& heavyPathGafFile);
void writePathGaf(const Path& path, const GfaGraph& graph, std::string outputFile);
std::vector<ReadPath> loadReadPaths(const std::string& filename, const GfaGraph& graph);
std::vector<ReadPath> extractCorrectedONTPaths(std::string gafFile, const Path& heavyPath, const size_t minLength, const GfaGraph& graph);
void writePathSequence(const Path& path, const GfaGraph& graph, std::string outputFile, const std::string& namePrefix);

#endif
