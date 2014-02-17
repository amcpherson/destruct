/*
 *  setcover.cpp
 *
 *  Created by Andrew McPherson on 10-08-24.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Parsers.h"
#include "Algorithms.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <tclap/CmdLine.h>

using namespace boost;
using namespace std;


int main(int argc, char* argv[])
{
	string clustersFilename;
	string outClustersFilename;
	bool distanceMin;
	
	try
	{
		TCLAP::CmdLine cmd("Set cover for maximum parsimony");
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> outClustersFilenameArg("o","outclust","Output Clusters Filename",true,"","string",cmd);
		TCLAP::SwitchArg distanceMinArg("d","distmin","Minimize distance",cmd);
		cmd.parse(argc,argv);
		
		clustersFilename = clustersFilenameArg.getValue();
		outClustersFilename = outClustersFilenameArg.getValue();
		distanceMin = distanceMinArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	double epsilon = 0.0;
	if (distanceMin)
	{
		epsilon = 0.00001;
	}
	
	cout << "Reading clusters" << endl;
	
	ClusterMembership clusterMembership(clustersFilename);
	
	IntegerVecMap clusters;
	IntegerMap distances;
	clusterMembership.Read(clusters, distances);
	
	double maxDistance = 1.0;
	for (IntegerMapConstIter distIter = distances.begin(); distIter != distances.end(); distIter++)
	{
		maxDistance = max(maxDistance, (double)distIter->second);
	}
	
	DoubleMap clusterWeights;
	for (IntegerVecMapConstIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		double clusterWeight = 1.0;
		
		IntegerMapConstIter distIter = distances.find(clusterIter->first);
		if (distIter != distances.end())
		{
			clusterWeight += epsilon * log((double)distIter->second);
		}
		else
		{
			clusterWeight += epsilon * log(maxDistance * 2.0);
		}
		
		clusterWeights[clusterIter->first] = clusterWeight;
	}
	
	cout << "Calculating set cover solution" << endl;

	IntegerVec solution;
	SetCover(clusters, clusterWeights, solution);
	
	cout << "Assigning fragments to solution sets" << endl;
	
	IntegerVecMap assignment;
	AssignInOrder(clusters, solution, assignment);
	
	cout << "Writing out clusters" << endl;
	
	clusterMembership.Write(outClustersFilename, assignment);
}

