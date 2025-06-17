#ifndef HifiasmIntegration_h
#define HifiasmIntegration_h

#include <vector>
#include <string>

std::string getHifiasmGraphFileName(std::string hifiasmBasePath);
std::vector<std::vector<std::string>> getHifiasmReadNamesPerTangle(std::string hifiasmGraphFile, const std::vector<std::vector<std::string>>& nodesPerTangle);
std::vector<std::vector<std::string>> guessHifiasmRDNATangles(std::string hifiasmGraphFile, const std::vector<std::string>& referencePath);

#endif
