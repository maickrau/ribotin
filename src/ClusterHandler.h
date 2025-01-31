#ifndef ClusterHandler_h
#define ClusterHandler_h

#include <string>

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
	size_t k;
	std::string orientReferencePath;
	std::string annotationFasta;
	std::string annotationGff3;
	std::string namePrefix;
	size_t numThreads;
	size_t maxClusterDifference;
	size_t maxResolveLength;
	size_t minReclusterDistance;
};

void HandleCluster(const ClusterParams& params);
void DoClusterONTAnalysis(const ClusterParams& params);
void AlignONTReads(std::string basePath, std::string graphAlignerPath, std::string ontReadPath, std::string graphPath, std::string outputAlnPath, size_t numThreads);

#endif
