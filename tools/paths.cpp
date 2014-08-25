/*
 *  paths.cpp
 */

#include "Common.h"
#include "DebugCheck.h"
#include "ShortestPath.h"
#include "BreakpointGraph.h"
#include "Indexer.h"

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_set.hpp>

using namespace boost;
using namespace std;


bool FindPath(BreakpointGraph* breakpointGraph, const LocationVec& rnaBreakpoint, NameIndex& references, int visitMax, double scoreMax, unordered_set<pair<unsigned int,int> > blocked, vector<pair<unsigned int,int> >& path, double& pathScore, int& numVisited)
{
	breakpointGraph->Reset();
	
	unsigned int refID1 = references.Index(rnaBreakpoint[0].refName);
	int strand1 = rnaBreakpoint[0].strand;
	unsigned int position1 = (strand1 == PlusStrand) ? rnaBreakpoint[0].end : rnaBreakpoint[0].start;
	
	unsigned int refID2 = references.Index(rnaBreakpoint[1].refName);
	int strand2 = rnaBreakpoint[1].strand;
	unsigned int position2 = (strand2 == PlusStrand) ? rnaBreakpoint[1].end : rnaBreakpoint[1].start;
	
	unsigned int startVertex;
	unsigned int endVertex;
	if (!breakpointGraph->SetPathSearch(refID1, strand1, position1, refID2, strand2, position2, startVertex, endVertex))
	{
		return false;
	}
	
	for (unordered_set<pair<unsigned int,int> >::const_iterator blockedIter = blocked.begin(); blockedIter != blocked.end(); blockedIter++)
	{
		breakpointGraph->BlockVertex(blockedIter->first, blockedIter->second);
	}
	
	vector<unsigned int> vertexPath;
	if (ShortestPath(breakpointGraph, startVertex, endVertex, visitMax, scoreMax, vertexPath, pathScore, numVisited))
	{
		for (vector<unsigned int>::const_iterator vertexIter = vertexPath.begin(); vertexIter != vertexPath.end(); vertexIter++)
		{
			if (startVertex == *vertexIter || endVertex == *vertexIter)
			{
				continue;
			}
			
			path.push_back(pair<unsigned int,int>(breakpointGraph->GetClusterID(*vertexIter),breakpointGraph->GetClusterEnd(*vertexIter)));
		}
		
		return true;
	}
	else
	{
		return false;
	}
}

void PrintPath(int rnaClusterID, const vector<pair<unsigned int,int> >& path, double pathScore, int numVisited)
{
	cerr << "Found path in " << numVisited << endl;
	
	cout << pathScore << "\t";
	
	cout << rnaClusterID << "\t";
	
	for (vector<pair<unsigned int,int> >::const_iterator clusterIDIter = path.begin(); clusterIDIter != path.end(); clusterIDIter++)
	{
		cout << clusterIDIter->first << "\t" << clusterIDIter->second << "\t";
	}
	
	cout << endl;	
}

inline bool operator==(const vector<pair<unsigned int,int> >& vec1, const vector<pair<unsigned int,int> >& vec2)
{
	return vec1.size() == vec2.size() && equal(vec1.begin(), vec1.end(), vec2.begin());
}

inline size_t hash_value(const vector<pair<unsigned int,int> >& vec)
{
	size_t seed = 0;
	for (vector<pair<unsigned int,int> >::const_iterator elementIter = vec.begin(); elementIter != vec.end(); elementIter++)
	{
		hash_combine(seed, *elementIter);
	}
    return seed;
}

inline bool operator==(const unordered_set<pair<unsigned int,int> >& set1, const unordered_set<pair<unsigned int,int> >& set2)
{
	return set1.size() == set2.size() && equal(set1.begin(), set1.end(), set2.begin());
}

struct unordered_set_hash
{
	inline size_t operator()(const unordered_set<pair<unsigned int,int> >& set) const
	{
		size_t seed = 0;
		for (unordered_set<pair<unsigned int,int> >::const_iterator elementIter = set.begin(); elementIter != set.end(); elementIter++)
		{
			hash_combine(seed, *elementIter);
		}
		return seed;
	}
};

template <typename TType>
void Print(const TType& a)
{
	for (typename TType::const_iterator elementIter = a.begin(); elementIter != a.end(); elementIter++)
	{
		cout << *elementIter << " ";
	}
	cout << endl;
}

int main(int argc, char* argv[])
{
	string dnaClustersFilename;
	string dnaProbsFilename;
	string rnaBreakpointsFilename;
	double distanceLambdaFwd;
	double distanceLambdaRev;
	int visitMax;
	double scoreMax;
	double altDiffMax;
	
	try
	{
		TCLAP::CmdLine cmd("DNA paths between RNA-Seq breakpoints");
		TCLAP::ValueArg<string> dnaClustersFilenameArg("d","dna","DNA Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> dnaProbsFilenameArg("p","probs","DNA Cluster Probabilities Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> rnaBreakpointsFilenameArg("r","rna","RNA-Seq Breakpoints Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> distanceLambdaFwdArg("y","lambdafwd","Distance Falloff Parameter Lambda Forward Direction",false,2000.0,"float",cmd);
		TCLAP::ValueArg<double> distanceLambdaRevArg("z","lambdarev","Distance Falloff Parameter Lambda Reverse Direction",false,20.0,"float",cmd);
		TCLAP::ValueArg<int> visitMaxArg("v","visitmax","Maximum Number of Vertices to Visit",false,100000,"integer",cmd);
		TCLAP::ValueArg<double> scoreMaxArg("s","scoremax","Maximum Score for Paths",false,30.0,"float",cmd);		
		TCLAP::ValueArg<double> altDiffMaxArg("a","altscorediffmax","Maximum Difference Between Best Path and Alternate Paths",false,-1.0,"float",cmd);		
		cmd.parse(argc,argv);
		
		dnaClustersFilename = dnaClustersFilenameArg.getValue();
		dnaProbsFilename = dnaProbsFilenameArg.getValue();
		rnaBreakpointsFilename = rnaBreakpointsFilenameArg.getValue();
		distanceLambdaFwd = distanceLambdaFwdArg.getValue();
		distanceLambdaRev = distanceLambdaRevArg.getValue();
		visitMax = visitMaxArg.getValue();
		scoreMax = scoreMaxArg.getValue();
		altDiffMax = altDiffMaxArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	cerr << "Reading DNA clusters" << endl;
	
	CompactLocationVecMap dnaClusters;
	NameIndex references;
	ReadClusters(dnaClustersFilename, dnaClusters, references);
	
	cerr << "Reading DNA cluster probabilities" << endl;
	
	DoubleMap dnaClusterProb;
	ReadDoubleMap(dnaProbsFilename, dnaClusterProb);
	
	cerr << "Creating DNA breakpoint graph" << endl;
	
	BreakpointGraph breakpointGraph(distanceLambdaFwd, distanceLambdaRev);
	for (CompactLocationVecMapConstIter dnaClusterIter = dnaClusters.begin(); dnaClusterIter != dnaClusters.end(); dnaClusterIter++)
	{
		DebugCheck(dnaClusterIter->second.size() == 2);
		
		unsigned int refID1 = dnaClusterIter->second[0].refStrand.referenceIndex;
		int strand1 = dnaClusterIter->second[0].refStrand.strand;
		unsigned int position1 = (strand1 == PlusStrand) ? dnaClusterIter->second[0].region.end : dnaClusterIter->second[0].region.start;
		
		unsigned int refID2 = dnaClusterIter->second[1].refStrand.referenceIndex;
		int strand2 = dnaClusterIter->second[1].refStrand.strand;
		unsigned int position2 = (strand2 == PlusStrand) ? dnaClusterIter->second[1].region.end : dnaClusterIter->second[1].region.start;
		
		breakpointGraph.AddBreakpoint(dnaClusterIter->first, refID1, strand1, position1, refID2, strand2, position2, -log(dnaClusterProb[dnaClusterIter->first]));
	}
	breakpointGraph.ConstructGraph();
	
	cerr << "Reading RNA breakpoints" << endl;
	
	LocationVecMap rnaBreakpoints;
	ReadAlignRegionPairs(rnaBreakpointsFilename, rnaBreakpoints);
	
	for (LocationVecMapConstIter rnaBreakpointIter = rnaBreakpoints.begin(); rnaBreakpointIter != rnaBreakpoints.end(); rnaBreakpointIter++)
	{
		cerr << "RNA cluster " << rnaBreakpointIter->first << endl;
		
		vector<unordered_set<pair<unsigned int,int> > > blockedStack;
		blockedStack.push_back(unordered_set<pair<unsigned int,int> >());
		
		double bestPathScore;
		
		unordered_set<vector<pair<unsigned int,int> > > allPaths;
		
		unordered_set<unordered_set<pair<unsigned int,int> >, unordered_set_hash> blockedSets;
		
		while (!blockedStack.empty())
		{
			unordered_set<pair<unsigned int,int> > blocked;
			swap(blockedStack.back(), blocked);
			blockedStack.pop_back();
			
			vector<pair<unsigned int,int> > path;
			double pathScore;
			int numVisited;
			
			if (!FindPath(&breakpointGraph, rnaBreakpointIter->second, references, visitMax, scoreMax, blocked, path, pathScore, numVisited) || pathScore > scoreMax)
			{
				continue;
			}
			else
			{
				if (blocked.empty())
				{
					bestPathScore = pathScore;
				}
				else if (pathScore > bestPathScore + altDiffMax)
				{
					continue;
				}
				
				if (allPaths.insert(path).second)
				{
					PrintPath(rnaBreakpointIter->first, path, pathScore, numVisited);
				}
				
				if (altDiffMax < 0.0)
				{
					break;
				}
				
				unordered_set<pair<unsigned int,int> > uniqueClusterIDs(path.begin(), path.end());
				for (unordered_set<pair<unsigned int,int> >::const_iterator clusterIDIter = uniqueClusterIDs.begin(); clusterIDIter != uniqueClusterIDs.end(); clusterIDIter++)
				{
					unordered_set<pair<unsigned int,int> > nextBlocked = blocked;
					nextBlocked.insert(*clusterIDIter);
					
					if (blockedSets.insert(nextBlocked).second)
					{
						blockedStack.push_back(unordered_set<pair<unsigned int,int> >());
						swap(blockedStack.back(), nextBlocked);
					}
				}
			}			
		}
	}
}


