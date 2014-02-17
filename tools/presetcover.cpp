/*
 *  presetcover.cpp
 *  tools
 *
 *  Created by Andrew McPherson.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Parsers.h"
#include "BinaryMinHeap.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/unordered_set.hpp>

using namespace boost;
using namespace boost::bimaps;
using namespace std;


void PreSetCover(const IntegerVecMap& sets, IntegerSet& solution)
{
	BinaryMinHeap minHeap;
	
	IntegerVecMap elementSets;
	IntegerMap setSizes;
	for (IntegerVecMapConstIter setIter = sets.begin(); setIter != sets.end(); setIter++)
	{
		setSizes[setIter->first] = (int)setIter->second.size();
		
		minHeap.Push(setIter->first, 1.0 / (double)setIter->second.size());
		
		for (IntegerVecConstIter elementIter = setIter->second.begin(); elementIter != setIter->second.end(); elementIter++)
		{
			elementSets[*elementIter].push_back(setIter->first);
		}
	}
	
	IntegerSet assigned;
	while (!minHeap.Empty())
	{
		IntegerSet nextSetIDs;
		
		double nextPriority = minHeap.MinPriority();
		
		nextSetIDs.insert(minHeap.Pop());
		
		while (!minHeap.Empty() && minHeap.MinPriority() == nextPriority)
		{
			nextSetIDs.insert(minHeap.Pop());
		}
		
		IntegerSet newAssigned;
		for (IntegerSetConstIter nextSetIter = nextSetIDs.begin(); nextSetIter != nextSetIDs.end(); nextSetIter++)
		{
			int nextSetID = *nextSetIter;
			
			solution.insert(nextSetID);
			
			const IntegerVec& set = sets.find(nextSetID)->second;
			
			newAssigned.insert(set.begin(), set.end());
		}
		
		IntegerSet alteredSets;
		for (IntegerSetConstIter elementIter = newAssigned.begin(); elementIter != newAssigned.end(); elementIter++)
		{
			if (assigned.insert(*elementIter).second)
			{
				for (IntegerVecConstIter setIter = elementSets[*elementIter].begin(); setIter != elementSets[*elementIter].end(); setIter++)
				{
					setSizes[*setIter]--;
					alteredSets.insert(*setIter);
				}
			}
		}
		
		for (IntegerSetConstIter setIter = alteredSets.begin(); setIter != alteredSets.end(); setIter++)
		{
			DebugCheck(setSizes[*setIter] >= 0);
			
			if (nextSetIDs.find(*setIter) != nextSetIDs.end())
			{
				continue;
			}
			
			if (setSizes[*setIter] == 0)
			{
				minHeap.Remove(*setIter);
			}
			else
			{
				minHeap.IncreaseKey(*setIter, 1.0 / (double)setSizes[*setIter]);
			}
		}
	}
}

int main(int argc, char* argv[])
{
	string clustersFilename;
	string outClustersFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Prefilter step for set cover for maximum parsimony");
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Clusters filename",true,"","string",cmd);
		TCLAP::ValueArg<string> outClustersFilenameArg("o","output","Output clusters filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		clustersFilename = clustersFilenameArg.getValue();
		outClustersFilename = outClustersFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	cout << "Reading clusters" << endl;
	
	ClusterMembership clusterMembership(clustersFilename);
	
	IntegerVecMap clusters;	
	clusterMembership.Read(clusters);
		
	cout << "Pre set cover" << endl;
	
	IntegerSet solution;
	PreSetCover(clusters, solution);
	
	cout << "Writing out clusters" << endl;
	
	clusterMembership.Write(outClustersFilename, solution);
}

