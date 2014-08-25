/*
 *  clusterconnectivity.cpp
 *
 *  Created by Andrew McPherson on 06-06-12.
 *
 */


#include "BinaryMinHeap.h"
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


class RegionDB
{
public:
	explicit RegionDB(int binSpacing = 10000) : mBinSpacing(binSpacing) {}
	
	void Add(int id, int ref, int start, int end)
	{
		for (int bin = start / mBinSpacing; bin <= end / mBinSpacing; bin++)
		{
			mBinned[RefBin(ref,bin)].push_back(IDStartEnd(id,start,end));
		}
	}
	
	void Overlapping(int ref, int start, int end, IntegerSet& overlapping) const
	{
		for (int bin = start / mBinSpacing; bin <= end / mBinSpacing; bin++)
		{
			unordered_map<RefBin,IDStartEndVec>::const_iterator findIter = mBinned.find(RefBin(ref,bin));
			
			if (findIter != mBinned.end())
			{
				for (IDStartEndVec::const_iterator regionIter = findIter->second.begin(); regionIter != findIter->second.end(); regionIter++)
				{
					if (regionIter->Overlaps(start,end))
					{
						overlapping.insert(regionIter->ID());
					}
				}
			}
		}
	}
	
private:
	typedef pair<int,int> RefBin;
	
	class IDStartEnd
	{
	public:
		IDStartEnd(int id, int start, int end) : mID(id), mStart(start), mEnd(end) {}
		
		bool Overlaps(int start, int end) const
		{
			return !(mEnd < start || end < mStart);
		}
		
		int ID() const { return mID; }
		
	private:
		int mID;
		int mStart;
		int mEnd;
	};
	
	typedef vector<IDStartEnd> IDStartEndVec;
	
	int mBinSpacing;
	unordered_map<RefBin,IDStartEndVec> mBinned;
};

void UnweightedSetCover(const IntegerVecMap& sets, unordered_multimap<int,int>& solution)
{
	BinaryMinHeap minHeap;
	
	IntegerVecMap elementSets;
	IntegerMap setSizes;
	for (IntegerVecMapConstIter setIter = sets.begin(); setIter != sets.end(); setIter++)
	{
		setSizes[setIter->first] = (int)setIter->second.size();
		
		for (IntegerVecConstIter elementIter = setIter->second.begin(); elementIter != setIter->second.end(); elementIter++)
		{
			elementSets[*elementIter].push_back(setIter->first);
			
			minHeap.Push(setIter->first, 1.0 / (double)setIter->second.size());
		}
	}
	
	IntegerSet assigned;
	while (!minHeap.Empty())
	{
		int nextSetID = minHeap.MinID();
		
		const IntegerVec& set = sets.find(nextSetID)->second;
		
		IntegerSet alteredSets;
		for (IntegerVecConstIter elementIter = set.begin(); elementIter != set.end(); elementIter++)
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
			
			solution.insert(make_pair(*setIter, nextSetID));
			
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

int NumDisjointOverlapping(const RegionDB& clusterRegionDB, const unordered_multimap<int,int>& setCover, int ref, int start, int end)
{
	IntegerSet overlapping;
	clusterRegionDB.Overlapping(ref, start, end, overlapping);
	
	IntegerSet overlappingSetCover;
	for (IntegerSetConstIter overlapIter = overlapping.begin(); overlapIter != overlapping.end(); overlapIter++)
	{
		pair<unordered_multimap<int,int>::const_iterator, unordered_multimap<int,int>::const_iterator> range = setCover.equal_range(*overlapIter);
		for (unordered_multimap<int,int>::const_iterator iter = range.first; iter != range.second; iter++)
		{
			overlappingSetCover.insert(iter->second);
		}
	}
	
	return overlappingSetCover.size();
}

double CalculateApproximateConnectivity(const unordered_map<int,vector<double> >& regionConnectivity, int resolution, int refID, int start, int end)
{
	int length = end - start + 1;
	
	int startBin = start / resolution;
	int endBin = end / resolution;
	
	double clusterConn = 0.0;
	for (int bin = startBin; bin <= endBin; bin++)
	{
		int binStart = bin * resolution;
		int binEnd = binStart + resolution - 1;
		
		int overlapStart = max(start, binStart);
		int overlapEnd = min(end, binEnd);
		int overlapLength = overlapEnd - overlapStart + 1;
		
		assert(overlapLength > 0.0);
		
		double proportion = (double)overlapLength / (double)length;
		
		clusterConn += proportion * regionConnectivity.find(refID)->second[bin];
	}
	
	return clusterConn;
}

template<typename ttype>
void print(const ttype& seq)
{
	for (typename ttype::const_iterator iter = seq.begin(); iter != seq.end(); iter++)
	{
		cout << *iter << " ";
	}
	cout << endl;
}

void UpdateMax(int& currentMax, int newValue)
{
	currentMax = max(currentMax, newValue);
}

int main(int argc, char* argv[])
{
	string connectFilename;
	string nullFilename;
	int numSamples;
	string restrictChromosome;
	
	try
	{
		TCLAP::CmdLine cmd("Number of overlapping clusters with distinct fragment sets");
		TCLAP::ValueArg<string> connectFilenameArg("c","connect","Connectivity Counts Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> nullFilenameArg("n","null","Null Connectivity Counts Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> numSamplesArg("s","samples","Number of Samples from Null",true,-1,"integer",cmd);
		TCLAP::ValueArg<string> restrictChromosomeArg("r","restrict","Restrict to the Given Chromosome",false,"","string",cmd);
		cmd.parse(argc,argv);
		
		connectFilename = connectFilenameArg.getValue();
		nullFilename = nullFilenameArg.getValue();
		numSamples = numSamplesArg.getValue();
		restrictChromosome = restrictChromosomeArg.getValue();
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
	
	cerr << "Reading all clusters" << endl;
	
	NameIndex referenceNames;
	IntegerMap referenceMaxClusterEnd;
	IntegerMap clusterReference[2];
	IntegerMap clusterRegionStart[2];
	IntegerMap clusterRegionEnd[2];
	int clusterMaxRegionLength = 0;
	RegionDB clusterRegionDB;
	vector<ClusterEndID> clusterEndIDs;
	IntegerVecMap clusters;
	
	{
		ClusterReader clusterReader(cin);
		
		while (clusterReader.Next())
		{
			int clusterID = clusterReader.FetchClusterID();
			LocationVec clusterLocations = clusterReader.FetchLocations();
			IntegerVec fragmentIndices = clusterReader.FetchFragmentIndices();
			
			for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
			{
				if (restrictChromosome != "" && restrictChromosome != clusterLocations[clusterEnd].refName)
				{
					continue;
				}
				
				ClusterEndID clusterEndID;
				clusterEndID.clusterID = clusterID;
				clusterEndID.clusterEnd = clusterEnd;
				
				clusterEndIDs.push_back(clusterEndID);
				
				RefStrand refStrand;
				refStrand.referenceIndex = referenceNames.Index(clusterLocations[clusterEnd].refName);
				refStrand.strand = clusterLocations[clusterEnd].strand;
				
				clusterReference[clusterEnd][clusterID] = refStrand.referenceIndex;
				clusterRegionStart[clusterEnd][clusterID] = clusterLocations[clusterEnd].start;
				clusterRegionEnd[clusterEnd][clusterID] = clusterLocations[clusterEnd].end;			
				
				clusterRegionDB.Add(clusterID, refStrand.referenceIndex, clusterLocations[clusterEnd].start, clusterLocations[clusterEnd].end);
				
				UpdateMax(referenceMaxClusterEnd.insert(make_pair(refStrand.referenceIndex,0)).first->second, clusterLocations[clusterEnd].end);
				UpdateMax(clusterMaxRegionLength, clusterLocations[clusterEnd].end - clusterLocations[clusterEnd].start + 1);
			}
			
			for (IntegerVecConstIter fragIter = fragmentIndices.begin(); fragIter != fragmentIndices.end(); fragIter++)
			{
				clusters[clusterID].push_back(*fragIter);
			}
		}
	}
	
	cerr << "Unweighted set cover" << endl;
	
	unordered_multimap<int,int> setCover;
	UnweightedSetCover(clusters, setCover);
	
	cerr << "Calculating cluster connectivity" << endl;
	
	int resolution = clusterMaxRegionLength / 2;
	unordered_map<int,vector<double> > regionConnectivity;
	for (int refID = 0; refID < referenceNames.Size(); refID++)
	{
		for (int start = 0; start < referenceMaxClusterEnd[refID]; start += resolution)
		{
			double numDisjoint = (double)NumDisjointOverlapping(clusterRegionDB, setCover, refID, start, start + resolution - 1);
			
			regionConnectivity[refID].push_back(numDisjoint);
		}
	}
	
	ofstream connectFile(connectFilename.c_str());
	CheckFile(connectFile, connectFilename);
	
	for (int refID = 0; refID < referenceNames.Size(); refID++)
	{
		for (int start = 1; start < referenceMaxClusterEnd[refID]; start += 500)
		{
			int numDisjoint = NumDisjointOverlapping(clusterRegionDB, setCover, refID, start, start + 500 - 1);
		}
	}

	for (vector<ClusterEndID>::const_iterator clusterEndIDIter = clusterEndIDs.begin(); clusterEndIDIter != clusterEndIDs.end(); clusterEndIDIter++)
	{
		int clusterID = clusterEndIDIter->clusterID;
		int clusterEnd = clusterEndIDIter->clusterEnd;
		
		int refID = clusterReference[clusterEnd][clusterID];
		int start = clusterRegionStart[clusterEnd][clusterID];
		int end = clusterRegionEnd[clusterEnd][clusterID];
		
		double clusterConn = CalculateApproximateConnectivity(regionConnectivity, resolution, refID, start, end);
		
		connectFile << clusterID << "\t" << clusterEnd << "\t" << clusterConn << endl;
	}
	
	connectFile.close();
	
	cerr << "Calculating null connectivity" << endl;
	
	RandomNumberGenerator rng;
	
	ofstream nullFile(nullFilename.c_str());
	CheckFile(nullFile, nullFilename);
	
	int samples = 0;
	while (samples < numSamples)
	{
		if (referenceNames.Size() == 0)
		{
			break;
		}
		
		int refID = rng.Next(0, referenceNames.Size() - 1);
		int refLength = referenceMaxClusterEnd[refID];
		
		int start = rng.Next(0, refLength - clusterMaxRegionLength - 1);
		int end = start + clusterMaxRegionLength;
		
		double clusterConn = CalculateApproximateConnectivity(regionConnectivity, resolution, refID, start, end);
		
		if (clusterConn == 0.0)
		{
			continue;
		}
		
		nullFile << clusterConn << endl;
		
		samples++;
	}
	
	nullFile.close();
}

