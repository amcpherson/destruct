/*
 *  cycles.cpp
 */

#include "Common.h"
#include "DebugCheck.h"
#include "AlignmentRecord.h"
#include "ShortestPath.h"
#include "BreakpointGraph.h"
#include "Indexer.h"

#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_set.hpp>

using namespace boost;
using namespace std;


bool FindPath(BreakpointGraph* breakpointGraph, int clusterID, double clusterScore, int visitMax, double scoreMax, unordered_set<pair<unsigned int,int> > blocked, vector<pair<unsigned int,int> >& cycle, double& cycleScore, int& numVisited)
{
	breakpointGraph->Reset();
	
	unsigned int startVertex;
	unsigned int endVertex;
	if (!breakpointGraph->SetCycleSearch(clusterID, startVertex, endVertex))
	{
		return false;
	}
	
	for (unordered_set<pair<unsigned int,int> >::const_iterator blockedIter = blocked.begin(); blockedIter != blocked.end(); blockedIter++)
	{
		breakpointGraph->BlockVertex(blockedIter->first, blockedIter->second);
	}
	
	vector<unsigned int> path;
	double pathScore;
	if (ShortestPath(breakpointGraph, startVertex, endVertex, visitMax, scoreMax - clusterScore, path, pathScore, numVisited))
	{
		cycleScore = pathScore + clusterScore;
		
		if (path.size() <= 2)
		{
			return false;
		}
		
		for (vector<unsigned int>::const_iterator vertexIter = path.begin(); vertexIter != path.end(); vertexIter++)
		{
			cycle.push_back(pair<unsigned int,int>(breakpointGraph->GetClusterID(*vertexIter),breakpointGraph->GetClusterEnd(*vertexIter)));
		}
		
		return true;
	}
	else
	{
		return false;
	}
}

void PrintCycle(const vector<pair<unsigned int,int> >& cycle, double cycleScore, int numVisited)
{
	cerr << "Found cycle in " << numVisited << endl;
	
	cout << cycleScore;
	
	for (vector<pair<unsigned int,int> >::const_iterator clusterIDIter = cycle.begin(); clusterIDIter != cycle.end(); clusterIDIter++)
	{
		cout << "\t" << clusterIDIter->first << "\t" << clusterIDIter->second;
	}
	
	cout << endl;	
}

int main(int argc, char* argv[])
{
	string breakpointsFilename;
	string probsFilename;
	int startClusterID;
	string startClusterIDsFilename;
	double startProbThreshold;
	double distanceLambda;
	int visitMax;
	double scoreMax;
	
	try
	{
		TCLAP::CmdLine cmd("Search for cycles in the breakpoint graph");
		TCLAP::ValueArg<string> breakpointsFilenameArg("b","breakpoints","Breakpoints Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> probsFilenameArg("p","probs","Cluster Probabilities Filename",false,"","string",cmd);
		TCLAP::ValueArg<int> startClusterIDArg("i","id","Starting Cluster ID",false,-1,"integer");
		TCLAP::ValueArg<string> startClusterIDsFilenameArg("","idsfile","Starting Cluster IDs Filename",false,"","string");
		TCLAP::ValueArg<double> startProbThresholdArg("t","threshold","Starting Cluster Probability Threshold",false,0.0,"float");
		TCLAP::ValueArg<double> distanceLambdaArg("y","lambda","Distance Falloff Parameter Lambda",false,2000.0,"float",cmd);
		TCLAP::ValueArg<int> visitMaxArg("v","visitmax","Maximum Number of Vertices to Visit",false,100000,"integer",cmd);
		TCLAP::ValueArg<double> scoreMaxArg("s","scoremax","Maximum Score for Cycles",false,30.0,"float",cmd);
		
		vector<TCLAP::Arg*> clusterIDArgs;
		clusterIDArgs.push_back(&startClusterIDArg);
		clusterIDArgs.push_back(&startClusterIDsFilenameArg);
		clusterIDArgs.push_back(&startProbThresholdArg);
		cmd.xorAdd(clusterIDArgs);
		
		cmd.parse(argc,argv);
		
		breakpointsFilename = breakpointsFilenameArg.getValue();
		probsFilename = probsFilenameArg.getValue();
		startClusterID = startClusterIDArg.getValue();
		startClusterIDsFilename = startClusterIDsFilenameArg.getValue();
		startProbThreshold = startProbThresholdArg.getValue();
		distanceLambda = distanceLambdaArg.getValue();
		visitMax = visitMaxArg.getValue();
		scoreMax = scoreMaxArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}

	unordered_set<int> startClustersIDs;
	if (startClusterID != -1)
	{
		startClustersIDs.insert(startClusterID);
	}
	else if (startClusterIDsFilename != "")
	{
		ifstream startClusterIDsFile(startClusterIDsFilename.c_str());
		CheckFile(startClusterIDsFile, startClusterIDsFilename);
		
		string line;
		int lineNumber = 0;
		while (getline(startClusterIDsFile, line))
		{
			lineNumber++;
			
			vector<string> fields;
			split(fields, line, is_any_of("\t"));
			
			if (fields.size() < 1)
			{
				cerr << "Error: Empty clusters line " << lineNumber << " of " << startClusterIDsFilename << endl;
				exit(1);
			}
			
			int clusterID = SAFEPARSE(int, fields[0]);
			
			startClustersIDs.insert(clusterID);
		}
		
		startClusterIDsFile.close();
	}
	
	cerr << "Creating breakpoint graph" << endl;
	
	ifstream breakpointsFile(breakpointsFilename.c_str());
	CheckFile(breakpointsFile, breakpointsFilename);

	NameIndex chromosomeIndex;

	BreakpointGraph breakpointGraph(distanceLambda, distanceLambda);

	int numBreakpoints = 0;

	BreakpointRecord breakpointRecord;
	while (breakpointsFile >> breakpointRecord)
	{
		// REVISIT
		// This is currently a cludge.  The cycle detection should not care about
		// the breakpoint predictions, only the approximate breakpoint for each
		// cluster.
		if (breakpointRecord.predictionID != 0)
		{
			continue;
		}

		unsigned int refID1 = chromosomeIndex.Index(breakpointRecord.chromosome[0]);
		unsigned int refID2 = chromosomeIndex.Index(breakpointRecord.chromosome[1]);

		int strand1 = (breakpointRecord.strand[0] == "+") ? PlusStrand : MinusStrand;
		int strand2 = (breakpointRecord.strand[1] == "+") ? PlusStrand : MinusStrand;

		unsigned int position1 = breakpointRecord.position[0];
		unsigned int position2 = breakpointRecord.position[1];
		
		// REVISIT
		// We are now not using cluster probabilities.  Not sure if this matters.
		breakpointGraph.AddBreakpoint(breakpointRecord.clusterID, refID1, strand1, position1, refID2, strand2, position2, 0.0);

		numBreakpoints++;
	}
	breakpointGraph.ConstructGraph();
	
	cerr << "Starting search on " << numBreakpoints << " clusters" << endl;
	
	for (unordered_set<int>::const_iterator clusterIDIter = startClustersIDs.begin(); clusterIDIter != startClustersIDs.end(); clusterIDIter++)
	{
		int clusterID = *clusterIDIter;
		
		unordered_set<pair<unsigned int,int> > blockedEmpty;
		
		vector<pair<unsigned int,int> > cycle;
		double cycleScore;
		int numVisited;
		
		// REVISIT
		// Cluster score again assumed to be 0.0
		if (!FindPath(&breakpointGraph, clusterID, 0.0, visitMax, scoreMax, blockedEmpty, cycle, cycleScore, numVisited))
		{
			cerr << "No Cycle" << endl;
			continue;
		}
		else
		{
			PrintCycle(cycle, cycleScore, numVisited);			
		}
	}
}


