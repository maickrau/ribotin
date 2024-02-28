#ifndef TangleGuesser_h
#define TangleGuesser_h

#include <vector>
#include <string>
#include "KmerMatcher.h"

std::vector<std::vector<std::string>> guessTangles(const KmerMatcher& matcher, const std::string& gfaPath);

#endif
