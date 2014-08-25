/*
 *  overlapclusters.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-08-11.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "Common.h"
#include "CompactBreakRegion.h"
#include "DebugCheck.h"

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


struct CompactSimpleRegion
{
	int ref;
	int start;
	int end;
};

typedef vector<CompactSimpleRegion> RegVec;
typedef vector<CompactSimpleRegion>::iterator RegVecIter;
typedef vector<CompactSimpleRegion>::const_iterator RegVecConstIter;

typedef unordered_map<int,IntegerVec> IntegerVecMap;
typedef unordered_map<int,IntegerVec>::iterator IntegerVecMapIter;
typedef unordered_map<int,IntegerVec>::const_iterator IntegerVecMapConstIter;

typedef unordered_set<IntegerPair> IntPairSet;
typedef unordered_set<IntegerPair>::const_iterator IntPairSetIter;

void ReadSimpleRegions(const string& regionsFilename, NameIndex& referenceNames, RegVec& regions)
{
	// Open regions file
	ifstream regionsFile(regionsFilename.c_str());
	if (!regionsFile)
	{
		cerr << "Error: unable to read from regions file " << regionsFilename << endl;		
		exit(1);
	}
	
	// Parse file contents
	string line;
	int lineNumber = 0;
	while (getline(regionsFile, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty break regions line " << lineNumber << endl;
			exit(1);
		}
		
		vector<string> regionFields;
		split(regionFields, line, is_any_of("\t"));
		
		if (regionFields.size() < 3)
		{
			cerr << "Error: Format error for regions line " << lineNumber << endl;
			exit(1);
		}
		
		string referenceName = regionFields[0];
		int start = SAFEPARSE(int, regionFields[1]);
		int end = SAFEPARSE(int, regionFields[2]);
		
		CompactSimpleRegion region;
		region.ref = referenceNames.Index(referenceName);
		region.start = start;
		region.end = end;
		
		regions.push_back(region);
	}
	
	regionsFile.close();
}

bool Overlap(const Region& r1, const Region& r2)
{
	return !(r1.end < r2.start || r2.end < r1.start);
}

class BinnedRegions
{
public:
	BinnedRegions(int binSpacing) : mBinSpacing(binSpacing) {}
	
	void Load(const RegVec& regions)
	{
		for (RegVecConstIter regIter = regions.begin(); regIter != regions.end(); regIter++)
		{
			IntegerVec bins;
			GetBins(regIter->start, regIter->end, bins);
			
			for (IntegerVecConstIter binIter = bins.begin(); binIter != bins.end(); binIter++)
			{
				mBinned[IntegerPair(regIter->ref,*binIter)].push_back(IntegerPair(regIter->start,regIter->end));
			}
		}
	}
	
	bool Contained(int ref, int start, int end) const
	{
		IntegerVec bins;
		GetBins(start, end, bins);
		
		for (IntegerVecConstIter binIter = bins.begin(); binIter != bins.end(); binIter++)
		{
			unordered_map<IntegerPair,IntegerPairVec>::const_iterator findIter = mBinned.find(IntegerPair(ref,*binIter));
			
			if (findIter != mBinned.end())
			{
				for (IntegerPairVecConstIter regIter = findIter->second.begin(); regIter != findIter->second.end(); regIter++)
				{
					if (Contains(*regIter,IntegerPair(start,end)))
					{
						return true;
					}
				}
			}			
		}
		
		return false;
	}
	
	void Clear()
	{
		mBinned.clear();
	}
	
private:
	void GetBins(int start, int end, IntegerVec& bins) const
	{
		int startBin = start / mBinSpacing;
		int endBin = end / mBinSpacing;
		
		if (startBin == endBin)
		{
			bins.push_back(startBin);
			bins.push_back(startBin + 1);
		}
		else
		{
			bins.push_back(startBin + 1);
		}	
	}
	
	bool Contains(IntegerPair region1, IntegerPair region2) const
	{
		if (region2.first >= region1.first && region2.second <= region1.second)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	
	int mBinSpacing;
	unordered_map<IntegerPair,IntegerPairVec> mBinned;
};

void SortClusters(IntegerTable& clusters)
{
	for (IntegerTableIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		sort(clusterIter->begin(), clusterIter->end());
	}
}

bool Overlaps(const IntegerVec& set1, const IntegerVec& set2)
{
	IntegerVecConstIter first1 = set1.begin();
	IntegerVecConstIter last1 = set1.end();
	
	IntegerVecConstIter first2 = set2.begin();
	IntegerVecConstIter last2 = set2.end();
	
	while (first1!=last1 && first2!=last2)
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

int main(int argc, char* argv[])
{
	string alignRegionsFilename;
	string clustersFilename;
	string repeatRegionsFilename;
	string outClustersFilename;
	int minClusterSize;

	try
	{
		TCLAP::CmdLine cmd("Filter repeat clusters tool");
		TCLAP::ValueArg<string> alignRegionsFilenameArg("a","aligns","Input Alignment Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Input Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> repeatRegionsFilenameArg("r","repeats","Input Repeat Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> outClustersFilenameArg("o","output","Output Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> minClusterSizeArg("m","minsize","Minimum Cluster Size",true,-1,"integer",cmd);
		cmd.parse(argc,argv);
		
		alignRegionsFilename = alignRegionsFilenameArg.getValue();
		clustersFilename = clustersFilenameArg.getValue();
		repeatRegionsFilename = repeatRegionsFilenameArg.getValue();
		outClustersFilename = outClustersFilenameArg.getValue();
		minClusterSize = minClusterSizeArg.getValue();
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
	
	cout << "Reading align regions" << endl;
	
	BrRegVec alignRegions;
	NameIndex referenceNames;
	ReadBreakRegions(alignRegionsFilename, referenceNames, alignRegions);
	
	if ((alignRegions.size() % 2) == 1)
	{
		cerr << "Error: align regions file has uneven number of align regions" << endl;		
		exit(1);
	}
	
	cout << "Reading clusters" << endl;
	
	IntegerTable clusters;
	ReadClusters(clustersFilename, clusters);
	SortClusters(clusters);
	
	cout << "Reading repeat regions" << endl;
	
	RegVec repeatRegions;
	ReadSimpleRegions(repeatRegionsFilename, referenceNames, repeatRegions);
	BinnedRegions binnedRepeatRegions(2000);
	binnedRepeatRegions.Load(repeatRegions);
	
	cout << "Finding clusters contained in repeats" << endl;
	
	unordered_set<int> repeatClusters;
	for (BrRegVecConstIter alignRegionIter = alignRegions.begin(); alignRegionIter != alignRegions.end(); alignRegionIter++)
	{
		if (binnedRepeatRegions.Contained(alignRegionIter->refStrand.referenceIndex, alignRegionIter->start, alignRegionIter->end))
		{
			repeatClusters.insert(alignRegionIter->clustEnd.clusterID);
		}
	}
	
	repeatRegions.clear();
	binnedRepeatRegions.Clear();
	
	cout << "Finding fragments contained in repeat clusters" << endl;
	
	unordered_set<int> repeatFragments;
	for (unordered_set<int>::const_iterator clusterIter = repeatClusters.begin(); clusterIter != repeatClusters.end(); clusterIter++)
	{
		for (IntegerVecConstIter fragmentIter = clusters[*clusterIter].begin(); fragmentIter != clusters[*clusterIter].end(); fragmentIter++)
		{
			repeatFragments.insert(*fragmentIter);
		}
	}
	
	cout << "Finding all clusters for repeat fragments" << endl;
	
	IntegerVecMap repeatFragmentClusters;
	for (int clusterIndex = 0; clusterIndex < clusters.size(); clusterIndex++)
	{
		for (IntegerVecConstIter fragmentIter = clusters[clusterIndex].begin(); fragmentIter != clusters[clusterIndex].end(); fragmentIter++)
		{
			if (repeatFragments.find(*fragmentIter) != repeatFragments.end())
			{
				repeatFragmentClusters[*fragmentIter].push_back(clusterIndex);
			}
		}
	}	
	
	cout << "Finding fragments contained only in repeat clusters" << endl;
	
	unordered_set<int> repeatOnlyFragments;
	for (IntegerVecMapIter fragmentIter = repeatFragmentClusters.begin(); fragmentIter != repeatFragmentClusters.end(); fragmentIter++)
	{
		bool foundNonRepeatCluster = false;
		for (IntegerVecConstIter clusterIter = fragmentIter->second.begin(); clusterIter != fragmentIter->second.end(); clusterIter++)
		{
			if (repeatClusters.find(*clusterIter) == repeatClusters.end())
			{
				foundNonRepeatCluster = true;
				break;
			}
		}
		
		if (!foundNonRepeatCluster)
		{
			repeatOnlyFragments.insert(fragmentIter->first);
		}
	}
	
	repeatFragmentClusters.clear();
	
	cout << "Finding clusters containing all repeat fragments" << endl;
	
	unordered_set<int> repeatOnlyClusters;
	for (unordered_set<int>::const_iterator clusterIter = repeatClusters.begin(); clusterIter != repeatClusters.end(); clusterIter++)
	{
		bool foundNonRepeatFragment = false;
		for (IntegerVecConstIter fragmentIter = clusters[*clusterIter].begin(); fragmentIter != clusters[*clusterIter].end(); fragmentIter++)
		{
			if (repeatOnlyFragments.find(*fragmentIter) == repeatOnlyFragments.end())
			{
				foundNonRepeatFragment = true;
				break;
			}
		}
		
		if (!foundNonRepeatFragment)
		{
			repeatOnlyClusters.insert(*clusterIter);
		}
	}
		
	cout << "Removing " << repeatOnlyClusters.size() << " clusters" << endl;
	
	for (int clusterIndex = 0; clusterIndex < clusters.size(); clusterIndex++)
	{
		if (repeatOnlyClusters.find(clusterIndex) != repeatOnlyClusters.end())
		{
			clusters[clusterIndex].clear();
		}
	}
	
	cout << "Writing out clusters" << endl;
	
	WriteClusters(clustersFilename, outClustersFilename, clusters, minClusterSize);
}



