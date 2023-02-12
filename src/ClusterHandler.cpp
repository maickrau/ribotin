#include <cstdlib>
#include <iostream>
#include "ClusterHandler.h"

void HandleCluster(std::string basePath, std::string readPath, std::string MBGPath)
{
	std::string mbgCommand;
	mbgCommand = "/usr/bin/time -v " + MBGPath + " -o " + basePath + "/graph.gfa -i " +readPath + " -k 51 -w 20 -a 2 -u 3 -r 15000 -R 4000 --error-masking=msat --output-sequence-paths " + basePath + "/paths.gaf --only-local-resolve 1> " + basePath + "/mbg_stdout.txt 2> " + basePath + "/mbg_stderr.txt";
	std::cerr << "running MBG with command:" << std::endl;
	std::cerr << mbgCommand << std::endl;
	system(mbgCommand.c_str());
}
