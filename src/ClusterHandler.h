#ifndef ClusterHandler_h
#define ClusterHandler_h

#include <string>
#include "ClusterMisc.h"

class ClusterParams
{
public:
	ClusterParams() = default;
	std::string basePath;
	std::string hifiReadPath;
	std::string ontReadPath;
	std::string MBGPath;
	std::string GraphAlignerPath;
	std::string winnowmapPath;
	std::string samtoolsPath;
	std::string liftoffPath;
	size_t k;
	std::string orientReferencePath;
	std::string annotationFasta;
	std::string annotationGff3;
	std::string namePrefix;
	size_t numThreads;
	size_t maxClusterDifference;
	size_t maxResolveLength;
	size_t minReclusterDistance;
	bool extraPhasing;
};

void HandleCluster(const ClusterParams& params);
std::vector<MorphConsensus> GetONTClusters(const ClusterParams& params);
void WriteONTClusters(const ClusterParams& params, const std::vector<MorphConsensus>& morphConsensuses);
void AlignONTReads(std::string basePath, std::string graphAlignerPath, std::string ontReadPath, std::string graphPath, std::string outputAlnPath, size_t numThreads);
void PostprocessONTClusters(std::vector<MorphConsensus>& morphConsensuses, const ClusterParams& params);
void runMultipleTangles(const ClusterParams& baseParams, const std::string& outputPrefix, const size_t numTangles, const std::vector<std::vector<std::string>>& hifiReadsPerTangle, const std::vector<std::string>& hifiReadFiles, const bool doUL, const std::string& ulTmpFolder, const std::vector<std::string>& ontReadFiles);

#endif
