#ifndef ClusterHandler_h
#define ClusterHandler_h

#include <string>

class ClusterParams
{
public:
	std::string basePath;
	std::string readPath;
	std::string MBGPath;
	size_t k;
	std::string orientReferencePath;
};

void HandleCluster(const ClusterParams& params);

#endif
