#ifndef HifiasmIntegration_h
#define HifiasmIntegration_h

#include <vector>
#include <string>

std::vector<std::vector<std::string>> getHifiasmReadNamesPerTangle(std::string hifiasmBasePath, const std::vector<std::vector<std::string>>& nodesPerTangle);
std::vector<std::vector<std::string>> guessHifiasmRDNATangles(std::string hifiasmBasePath, const std::vector<std::string>& referencePath);

#endif
