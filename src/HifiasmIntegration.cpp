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
	std::cerr << "checking for " << hifiasmBasePath << ".hic.r_utg.gfa" << std::endl;
	if (fileExists(hifiasmBasePath + ".hic.r_utg.gfa"))
	{
		return hifiasmBasePath + ".hic.r_utg.gfa";
	}
	std::cerr << "checking for " << hifiasmBasePath << ".hic.bench.r_utg.gfa" << std::endl;
	if (fileExists(hifiasmBasePath + ".hic.bench.r_utg.gfa"))
	{
		return hifiasmBasePath + ".hic.bench.r_utg.gfa";
	}
	std::cerr << "checking for " << hifiasmBasePath << ".r_utg.gfa" << std::endl;
	if (fileExists(hifiasmBasePath + ".r_utg.gfa"))
	{
		return hifiasmBasePath + ".r_utg.gfa";
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
	phmap::flat_hash_map<std::string, size_t> nodelens;
	std::vector<std::tuple<std::string, std::string, size_t, size_t>> readsInTipNodes;
	while (file.good())
	{
		std::string line;
		getline(file, line);
		if (!file.good()) break;
		if (line[0] == 'S')
		{
			std::string dummy, nodename, sequence, lengthstr;
			std::stringstream sstr { line };
			sstr >> dummy >> nodename >> sequence >> lengthstr;
			size_t length = 0;
			if (sequence != "*")
			{
				length = sequence.size();
			}
			else
			{
				assert(lengthstr.substr(0, 5) == "LN:i:");
				length = std::stoull(lengthstr.substr(5));
			}
			assert(length > 0);
			nodelens[nodename] = length;
		}
		if (line[0] != 'A') continue;
		std::string dummy, nodename, readname;
		size_t position, readlen;
		std::stringstream sstr { line };
		sstr >> dummy >> nodename >> position >> dummy >> readname >> dummy >> readlen;
		if (nodeToTangle.count(nodename) == 1)
		{
			size_t tangle = nodeToTangle.at(nodename);
			result[tangle].emplace_back(readname);
		}
		else if (nodeToTangle.count(">" + nodename) == 1 || nodeToTangle.count("<" + nodename) == 1)
		{
			readsInTipNodes.emplace_back(readname, nodename, position, position+readlen);
		}
	}
	for (const auto& t : readsInTipNodes)
	{
		assert(nodeToTangle.count(">" + std::get<1>(t)) == 1 || nodeToTangle.count("<" + std::get<1>(t)) == 1);
		size_t validStart = 0;
		size_t validEnd = 50000;
		size_t tangle;
		if (nodeToTangle.count(">" + std::get<1>(t)) == 1)
		{
			tangle = nodeToTangle.at(">" + std::get<1>(t));
			validEnd = nodelens[std::get<1>(t)];
			assert(validEnd >= 50000);
			validStart = validEnd - 50000;
		}
		else
		{
			tangle = nodeToTangle.at("<" + std::get<1>(t));
		}
		if (std::get<3>(t) < validStart) continue;
		if (std::get<2>(t) > validEnd) continue;
		result[tangle].emplace_back(std::get<0>(t));
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
