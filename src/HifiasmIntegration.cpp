#include <phmap.h>
#include "KmerMatcher.h"
#include "HifiasmIntegration.h"
#include "TangleGuesser.h"

bool fileExists(const std::string& filename)
{
	std::ifstream file { filename };
	return file.good();
}

std::string getHifiasmGraphFileName(std::string hifiasmBasePath)
{
	std::cerr << "checking for " << hifiasmBasePath << ".dip.r_utg.gfa" << std::endl;
	if (fileExists(hifiasmBasePath + ".dip.r_utg.gfa"))
	{
		return hifiasmBasePath + ".dip.r_utg.gfa";
	}
	std::cerr << "checking for " << hifiasmBasePath << ".bp.r_utg.gfa" << std::endl;
	if (fileExists(hifiasmBasePath + ".bp.r_utg.gfa"))
	{
		return hifiasmBasePath + ".bp.r_utg.gfa";
	}
	std::cerr << "ERROR: Did not find hifiasm assembly graph" << std::endl;
	std::abort();
}

std::vector<std::vector<std::string>> getHifiasmReadNamesPerTangle(std::string hifiasmGraphFile, const std::vector<std::vector<std::string>>& nodesPerTangle)
{
	phmap::flat_hash_map<std::string, size_t> nodeToTangle;
	for (size_t i = 0; i < nodesPerTangle.size(); i++)
	{
		for (const std::string& node : nodesPerTangle[i])
		{
			assert(nodeToTangle.count(node) == 0);
			nodeToTangle[node] = i;
		}
	}
	std::vector<std::vector<std::string>> result;
	result.resize(nodesPerTangle.size());
	std::ifstream file { hifiasmGraphFile };
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (!file.good()) break;
		if (line[0] != 'A') continue;
		std::string dummy, nodename, readname;
		std::stringstream sstr { line };
		sstr >> dummy >> nodename >> dummy >> dummy >> readname;
		if (nodeToTangle.count(nodename) == 1)
		{
			size_t tangle = nodeToTangle.at(nodename);
			result[tangle].emplace_back(readname);
		}
	}
	return result;
}

std::vector<std::vector<std::string>> guessHifiasmRDNATangles(std::string hifiasmGraphFile, const std::vector<std::string>& referencePath)
{
	KmerMatcher matcher { 101 };
	for (auto file : referencePath)
	{
		FastQ::streamFastqFromFile(file, false, [&matcher](FastQ& fastq)
		{
			matcher.addReferenceKmers(fastq.sequence);
		});
	}
	auto result = guessTangles(matcher, hifiasmGraphFile);
	return result;
}
