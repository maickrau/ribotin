#ifndef VerkkoClusterGuesser_h
#define VerkkoClusterGuesser_h

#include <vector>
#include <string>

std::vector<std::vector<std::string>> guessVerkkoRDNAClusters(std::string verkkoBasePath, const std::vector<std::string>& referencePath);

#endif
