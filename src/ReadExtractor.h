#ifndef ReadExtractor_h
#define ReadExtractor_h

#include <vector>
#include <string>

std::string nameWithoutTags(const std::string& rawname);
bool splitReads(std::vector<std::string> readFiles, const std::vector<std::vector<std::string>>& readsPerCluster, const std::vector<std::string>& outputFileNames);

#endif
