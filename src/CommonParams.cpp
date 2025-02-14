#include <iostream>
#include <cassert>
#include "CommonParams.h"

CommonParams::CommonParams() :
	GraphAlignerPathpriv(""),
	MBGPath(""),
	winnowmapPath(""),
	samtoolsPath(""),
	liftoffPath(""),
	sampleName(""),
	orientReferencePath(""),
	maxResolveLength(0),
	morphClusterMaxDistance(0),
	morphReclusterMinDistance(0),
	numThreadspriv(0),
	k(0),
	annotationReferenceFastaAndGff3("", "")
{
}

void CommonParams::addParamsToOptions(cxxopts::Options& options)
{
	options.add_options()
		("t", "Number of threads", cxxopts::value<size_t>()->default_value("1"))
		("x", "Preset parameters", cxxopts::value<std::string>())
		("sample-name", "Name of the sample added to all morph names", cxxopts::value<std::string>())
		("orient-by-reference", "Rotate and possibly reverse complement the consensus to match the orientation of the given reference", cxxopts::value<std::string>())
		("approx-morphsize", "Approximate length of one morph", cxxopts::value<size_t>())
		("k", "k-mer size", cxxopts::value<size_t>()->default_value("101"))
		("annotation-reference-fasta", "Lift over the annotations from given reference fasta+gff3 (requires liftoff)", cxxopts::value<std::string>())
		("annotation-gff3", "Lift over the annotations from given reference fasta+gff3 (requires liftoff)", cxxopts::value<std::string>())
		("morph-cluster-maxedit", "Maximum edit distance between two morphs to assign them into the same cluster", cxxopts::value<size_t>()->default_value("200"))
		("morph-recluster-minedit", "Minimum edit distance to recluster morphs", cxxopts::value<size_t>()->default_value("5"))
		("mbg", "MBG path", cxxopts::value<std::string>())
		("graphaligner", "GraphAligner path", cxxopts::value<std::string>())
		("winnowmap", "winnowmap/minimap path", cxxopts::value<std::string>())
		("samtools", "samtools path", cxxopts::value<std::string>())
	;
}

bool CommonParams::parseParamsAndPrintErrors(const cxxopts::ParseResult& params)
{
	bool paramError = false;
	if (params.count("x") == 1)
	{
		if (params["x"].as<std::string>() != "human")
		{
			std::cerr << "No preset \"" << params["x"].as<std::string>() << "\"" << std::endl;
			paramError = true;
		}
	}
	if (params.count("mbg") == 0)
	{
		std::cerr << "checking for MBG" << std::endl;
		int foundMBG = system("which MBG");
		if (foundMBG != 0)
		{
			std::cerr << "MBG not found" << std::endl;
			std::cerr << "MBG path (--mbg) is required" << std::endl;
			paramError = true;
		}
		else
		{
			MBGPath = "MBG";
		}
	}
	else
	{
		MBGPath = params["mbg"].as<std::string>();
	}
	if (params.count("graphaligner") == 1)
	{
		GraphAlignerPathpriv = params["graphaligner"].as<std::string>();
	}
	if (params.count("graphaligner") == 0)
	{
		std::cerr << "checking for GraphAligner" << std::endl;
		int foundGraphAligner = system("which GraphAligner");
		if (foundGraphAligner != 0)
		{
			std::cerr << "GraphAligner not found" << std::endl;
		}
		else
		{
			GraphAlignerPathpriv = "GraphAligner";
		}
	}
	if (params.count("winnowmap") == 1)
	{
		winnowmapPath = params["winnowmap"].as<std::string>();
	}
	if (params.count("winnowmap") == 0)
	{
		std::cerr << "checking for winnowmap" << std::endl;
		int foundWinnowmap = system("which winnowmap");
		if (foundWinnowmap != 0)
		{
			std::cerr << "winnowmap not found" << std::endl;
		}
		else
		{
			winnowmapPath = "winnowmap";
		}
	}
	if (params.count("samtools") == 1)
	{
		samtoolsPath = params["samtools"].as<std::string>();
	}
	if (params.count("samtools") == 0)
	{
		std::cerr << "checking for samtools" << std::endl;
		int foundSamtools = system("which samtools");
		if (foundSamtools != 0)
		{
			std::cerr << "samtools not found" << std::endl;
		}
		else
		{
			samtoolsPath = "samtools";
		}
	}
	if (params.count("approx-morphsize") == 0 && params.count("x") == 0)
	{
		std::cerr << "Approximate size of one morph (--approx-morphsize) is required" << std::endl;
		paramError = true;
	}
	if (params.count("t") == 1 && params["t"].as<size_t>() == 0)
	{
		std::cerr << "number of threads can't be 0" << std::endl;
		paramError = true;
	}
	if (params.count("k") == 1 && params["k"].as<size_t>() < 31)
	{
		std::cerr << "k must be at least 31" << std::endl;
		paramError = true;
	}
	if (params.count("k") == 1 && params["k"].as<size_t>() % 2 == 0)
	{
		std::cerr << "k must be odd" << std::endl;
		paramError = true;
	}
	if (params.count("annotation-gff3") == 1 && params.count("annotation-reference-fasta") == 0)
	{
		std::cerr << "--annotation-reference-fasta is missing while --annotation-gff3 is used" << std::endl;
		paramError = true;
	}
	if (params.count("annotation-gff3") == 0 && params.count("annotation-reference-fasta") == 1)
	{
		std::cerr << "--annotation-gff3 is missing while --annotation-reference-fasta is used" << std::endl;
		paramError = true;
	}
	if (params.count("liftoff") == 0)
	{
		std::cerr << "checking for liftoff" << std::endl;
		int foundLiftoff = system("which liftoff");
		if (foundLiftoff != 0)
		{
			std::cerr << "liftoff not found" << std::endl;
		}
		else
		{
			liftoffPath = "liftoff";
		}
	}
	else
	{
		liftoffPath = params["liftoff"].as<std::string>();
	}
	if (liftoffPath == "")
	{
		if (params.count("annotation-gff3") == 1 || params.count("annotation-reference-fasta") == 1 || (params.count("x") == 1 && params["x"].as<std::string>() == "human"))
		{
			paramError = true;
		}
	}
	if (paramError) return false;
	k = params["k"].as<size_t>();
	morphClusterMaxDistance = params["morph-cluster-maxedit"].as<size_t>();
	morphReclusterMinDistance = params["morph-recluster-minedit"].as<size_t>();
	if (params.count("x") == 1)
	{
		if (params["x"].as<std::string>() == "human")
		{
			k = 101;
			morphClusterMaxDistance = 200;
			morphReclusterMinDistance = 5;
			maxResolveLength = 1000;
			orientReferencePath = std::string{RIBOTIN_TEMPLATE_PATH} + "/rDNA_one_unit.fasta";
			annotationReferenceFastaAndGff3 = std::make_pair(std::string{RIBOTIN_TEMPLATE_PATH} + "/rDNA_one_unit.fasta", std::string{RIBOTIN_TEMPLATE_PATH} + "/rDNA_annotation.gff3");
		}
		if (params.count("k") == 1) k = params["k"].as<size_t>();
		if (params.count("morph-cluster-maxedit") == 1) morphClusterMaxDistance = params["morph-cluster-maxedit"].as<size_t>();
		if (params.count("morph-recluster-minedit") == 1) morphReclusterMinDistance = params["morph-recluster-minedit"].as<size_t>();
	}
	numThreadspriv = params["t"].as<size_t>();
	if (params.count("approx-morphsize") == 1) maxResolveLength = 100 + params["approx-morphsize"].as<size_t>()/50;
	if (params.count("sample-name") == 1) sampleName = params["sample-name"].as<std::string>();
	if (params.count("orient-by-reference") == 1) orientReferencePath = params["orient-by-reference"].as<std::string>();
	if (params.count("annotation-reference-fasta") == 1)
	{
		assert(params.count("annotation-gff3") == 1);
		annotationReferenceFastaAndGff3 = std::make_pair(params["annotation-reference-fasta"].as<std::string>(), params["annotation-gff3"].as<std::string>());
	}
	else
	{
		assert(params.count("annotation-gff3") == 0);
	}
	return true;
}

void CommonParams::addToClusterOptions(ClusterParams& clusterParams)
{
	clusterParams.GraphAlignerPath = GraphAlignerPath();
	clusterParams.MBGPath = MBGPath;
	clusterParams.winnowmapPath = winnowmapPath;
	clusterParams.samtoolsPath = samtoolsPath;
	clusterParams.k = k;
	clusterParams.orientReferencePath = orientReferencePath;
	clusterParams.liftoffPath = liftoffPath;
	clusterParams.annotationFasta = annotationReferenceFastaAndGff3.first;
	clusterParams.annotationGff3 = annotationReferenceFastaAndGff3.second;
	clusterParams.namePrefix = sampleName;
	clusterParams.numThreads = numThreads();
	clusterParams.maxClusterDifference = morphClusterMaxDistance;
	clusterParams.minReclusterDistance = morphReclusterMinDistance;
	clusterParams.maxResolveLength = maxResolveLength;
}

size_t CommonParams::numThreads() const
{
	return numThreadspriv;
}

std::string CommonParams::GraphAlignerPath() const
{
	return GraphAlignerPathpriv;
}

bool CommonParams::hasGraphAligner() const
{
	return GraphAlignerPath() != "";
}

bool CommonParams::hasWinnowmap() const
{
	return winnowmapPath != "";
}

bool CommonParams::hasSamtools() const
{
	return samtoolsPath != "";
}
