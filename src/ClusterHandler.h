#ifndef ClusterHandler_h
#define ClusterHandler_h

#include <string>

class ClusterParams
{
public:
	ClusterParams() = default;
	std::string basePath;
	std::string readPath;
	std::string MBGPath;
	size_t k;
	std::string orientReferencePath;
	std::string annotationFasta;
	std::string annotationGff3;
};

void HandleCluster(const ClusterParams& params);

#endif
