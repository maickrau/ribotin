#ifndef ClusterHandler_h
#define ClusterHandler_h

#include <string>
#include "ClusterMisc.h"

void HandleCluster(const ClusterParams& params);
void DoClusterONTAnalysis(const ClusterParams& params);
void AlignONTReads(std::string basePath, std::string graphAlignerPath, std::string ontReadPath, std::string graphPath, std::string outputAlnPath, size_t numThreads);

#endif
