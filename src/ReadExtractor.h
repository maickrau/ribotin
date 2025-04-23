#ifndef ReadExtractor_h
#define ReadExtractor_h

#include <vector>
#include <string>

bool splitReads(std::vector<std::string> readFiles, const std::vector<std::vector<std::string>>& readsPerCluster, const std::vector<std::string>& outputFileNames);
void extractMoreHifiReads(std::vector<std::string> wholeGenomeHifiFastas, const std::vector<std::string>& tempReadFileNames, const std::vector<std::string>& realReadFileNames, const std::string MBGPath, const std::string tmpFolder);

#endif
