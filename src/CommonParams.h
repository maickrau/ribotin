#ifndef CommonParams_h
#define CommonParams_h

#include <string>
#include <tuple>
#include <cstddef>
#include <cxxopts.hpp>
#include "ClusterHandler.h"

class CommonParams
{
public:
	CommonParams();
	void addParamsToOptions(cxxopts::Options& options);
	bool parseParamsAndPrintErrors(const cxxopts::ParseResult& params);
	void addToClusterOptions(ClusterParams& clusterParams);
	bool hasGraphAligner() const;
	bool hasWinnowmap() const;
	bool hasSamtools() const;
	std::string GraphAlignerPath() const;
	size_t numThreads() const;
private:
	std::string GraphAlignerPathpriv;
	std::string MBGPath;
	std::string winnowmapPath;
	std::string samtoolsPath;
	std::string liftoffPath;
	std::string sampleName;
	std::string orientReferencePath;
	size_t maxResolveLength;
	size_t morphClusterMaxDistance;
	size_t morphReclusterMinDistance;
	size_t numThreadspriv;
	size_t k;
	bool extraPhasing;
	std::pair<std::string, std::string> annotationReferenceFastaAndGff3;
};

#endif
