/*
 *  overlapclusters.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-08-11.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "CompactBreakRegion.h"
#include "DebugCheck.h"
#include "Matrix.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;


template<typename TIterator1, typename TIterator2>
bool intersects(TIterator1 first1, TIterator1 last1, TIterator2 first2, TIterator2 last2)
{
	while (first1 != last1 && first2 != last2)
	{
		if (*first1 < *first2)
		{
			++first1;
		}
		else if (*first2 < *first1)
		{
			++first2;
		}
		else
		{
			return true;
		}
	}
	
	return false;
}

bool Overlap(const CompactLocation& r1, const CompactLocation& r2)
{
	return (r1.refStrand.id == r2.refStrand.id && !(r1.region.end < r2.region.start || r2.region.end < r1.region.start));
}

bool Overlap(const CompactLocationVec& regions1, const CompactLocationVec& regions2)
{
	DebugCheck(regions1.size() == 2 && regions2.size() == 2);

	Matrix<int> overlap(2,2);
	overlap.Clear(0);

	for (int regionIndex1 = 0; regionIndex1 <= 1; regionIndex1++)
	{
		for (int regionIndex2 = 0; regionIndex2 <= 1; regionIndex2++)
		{
			if (Overlap(regions1[regionIndex1], regions2[regionIndex2]))
			{
				overlap(regionIndex1,regionIndex2) = 1;
			}
		}
	}
	
	return overlap(0,0) && overlap(1,1) || overlap(0,1) && overlap(1,0);
}

int CalculateNumDisjointComponents(unordered_map<int,IntegerVec>& overlapGraph)
{
	int componentCount = 0;
	
	IntegerVec searchStack;
	while (!overlapGraph.empty())
	{
		if (searchStack.empty())
		{
			searchStack.push_back(overlapGraph.begin()->first);
			componentCount++;
		}
		
		int currentElement = searchStack.back();
		searchStack.pop_back();
		
		if (overlapGraph.find(currentElement) != overlapGraph.end())
		{
			searchStack.insert(searchStack.end(), overlapGraph.find(currentElement)->second.begin(), overlapGraph.find(currentElement)->second.end());
			overlapGraph.erase(currentElement);
		}
	}
	
	return componentCount;
}

int main(int argc, char* argv[])
{
	string allClustersFilename;
	string solutionClustersFilename;
	string decisionsFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Overlap between break regions tool");
		TCLAP::ValueArg<string> allClustersFilenameArg("a","all","Input All Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> solutionClustersFilenameArg("s","solv","Input Solution Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> decisionsFilenameArg("d","decisions","Output Decision Metrics Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		allClustersFilename = allClustersFilenameArg.getValue();
		solutionClustersFilename = solutionClustersFilenameArg.getValue();
		decisionsFilename = decisionsFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	if (sizeof(long) != 8)
	{
		cerr << "Error: sizeof(long) != 8" << endl;		
		exit(1);
	}
	
	cout << "Reading all clusters" << endl;
	
	NameIndex referenceNames;
	
	IntegerVecMap allClusters;
	IntegerVecMap allFragmentClusters;
	unordered_map<int,CompactLocationVec> allClusterRegions;
	
	{
		ifstream allClustersFile(allClustersFilename.c_str());
		CheckFile(allClustersFile, allClustersFilename);
		
		ClusterReader clusterReader(allClustersFile);
		
		while (clusterReader.Next())
		{
			int clusterID = clusterReader.FetchClusterID();
			LocationVec clusterLocations = clusterReader.FetchLocations();
			IntegerVec fragmentIndices = clusterReader.FetchFragmentIndices();
			
			allClusterRegions[clusterID] = CompactLocationVec(2);
			
			for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
			{
				allClusterRegions[clusterID][clusterEnd].refStrand.referenceIndex = referenceNames.Index(clusterLocations[clusterEnd].refName);
				allClusterRegions[clusterID][clusterEnd].refStrand.strand = clusterLocations[clusterEnd].strand;
				allClusterRegions[clusterID][clusterEnd].region.start = clusterLocations[clusterEnd].start;
				allClusterRegions[clusterID][clusterEnd].region.end = clusterLocations[clusterEnd].end;
			}
			
			for (IntegerVecConstIter fragIter = fragmentIndices.begin(); fragIter != fragmentIndices.end(); fragIter++)
			{
				allClusters[clusterID].push_back(*fragIter);
				allFragmentClusters[*fragIter].push_back(clusterID);
			}
		}
	}
	
	cout << "Calculating number of competing clusters" << endl;
	
	ofstream decisionsFile(decisionsFilename.c_str());
	CheckFile(decisionsFile, decisionsFilename);
	
	{
		ifstream solutionClustersFile(solutionClustersFilename.c_str());
		CheckFile(solutionClustersFile, solutionClustersFilename);
		
		ClusterReader clusterReader(solutionClustersFile);
		
		while (clusterReader.Next())
		{
			int clusterID = clusterReader.FetchClusterID();
			LocationVec clusterLocations = clusterReader.FetchLocations();
			IntegerVec fragmentIndices = clusterReader.FetchFragmentIndices();
			
			CompactLocationVec location(2);
			for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
			{
				location[clusterEnd].refStrand.referenceIndex = referenceNames.Index(clusterLocations[clusterEnd].refName);
				location[clusterEnd].refStrand.strand = clusterLocations[clusterEnd].strand;
				location[clusterEnd].region.start = clusterLocations[clusterEnd].start;
				location[clusterEnd].region.end = clusterLocations[clusterEnd].end;
			}
			
			IntegerSet overlapClusters;
			for (IntegerVecConstIter fragIter = fragmentIndices.begin(); fragIter != fragmentIndices.end(); fragIter++)
			{
				overlapClusters.insert(allFragmentClusters[*fragIter].begin(), allFragmentClusters[*fragIter].end());
			}
			
			IntegerSet competingClusters;
			for (IntegerSetConstIter clusterIter = overlapClusters.begin(); clusterIter != overlapClusters.end(); clusterIter++)
			{
				IntegerVec& overlapFragments = allClusters[*clusterIter];
				CompactLocationVec& overlapLocation = allClusterRegions[*clusterIter];
				
				bool eq = fragmentIndices.size() == overlapFragments.size() && equal(fragmentIndices.begin(), fragmentIndices.end(), overlapFragments.begin());
				bool incl = includes(fragmentIndices.begin(), fragmentIndices.end(), overlapFragments.begin(), overlapFragments.end());
				bool olap = intersects(fragmentIndices.begin(), fragmentIndices.end(), overlapFragments.begin(), overlapFragments.end());
				bool equiv = Overlap(location, overlapLocation);
				
				if (!equiv && (eq || !incl && olap))
				{
					competingClusters.insert(*clusterIter);
				}
			}
			
			unordered_map<int,IntegerVec> overlapGraph;
			for (IntegerSetConstIter cluster1Iter = competingClusters.begin(); cluster1Iter != competingClusters.end(); cluster1Iter++)
			{
				CompactLocationVec& location1 = allClusterRegions[*cluster1Iter];
				
				overlapGraph.insert(make_pair(*cluster1Iter, IntegerVec()));
				
				IntegerSetConstIter cluster2Iter = cluster1Iter;
				cluster2Iter++;
				for (; cluster2Iter != competingClusters.end(); cluster2Iter++)
				{
					CompactLocationVec& location2 = allClusterRegions[*cluster2Iter];
					
					if (Overlap(location1,location2))
					{
						overlapGraph[*cluster1Iter].push_back(*cluster2Iter);
						overlapGraph[*cluster2Iter].push_back(*cluster1Iter);
					}
				}
			}
			
			int numComponents = CalculateNumDisjointComponents(overlapGraph);
			
			decisionsFile << clusterID << "\t" << numComponents + 1 << endl;
		}
	}
	
	decisionsFile.close();
}

