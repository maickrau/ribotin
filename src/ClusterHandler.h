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
	size_t k;
	std::string orientReferencePath;
	std::string annotationFasta;
	std::string annotationGff3;
	size_t numThreads;
};

void HandleCluster(const ClusterParams& params);

#endif
