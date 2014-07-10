/*
 *  setcover.cpp
 *
 *  Created by Andrew McPherson on 10-08-24.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Algorithms.h"
#include "AlignmentRecord.h"

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
	string weightsFilename;
	string assignmentsFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Set cover for maximum parsimony");
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> weightsFilenameArg("w","weights","Weights Filename",false,"","string",cmd);
		TCLAP::ValueArg<string> assignmentsFilenameArg("a","assignments","Output Assignments Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		clustersFilename = clustersFilenameArg.getValue();
		weightsFilename = weightsFilenameArg.getValue();
		assignmentsFilename = assignmentsFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	cout << "Reading clusters" << endl;
	
	vector<int> ids;
	vector<vector<ReadRecord> > clusters;

	{
		ifstream clustersFile(clustersFilename.c_str());
		CheckFile(clustersFile, clustersFilename);

		GroupedRecordsStream<ClusterMemberRecord> memberStream(clustersFile);

		vector<ClusterMemberRecord> clusterRecords;
		while (memberStream.Next(clusterRecords, ClusterEqual<ClusterMemberRecord>))
		{
			DebugCheck(clusterRecords.size() > 0);

			ids.push_back(clusterRecords.front().clusterID);
			clusters.push_back(vector<ReadRecord>());

			for (vector<ClusterMemberRecord>::const_iterator recordIter = clusterRecords.begin(); recordIter != clusterRecords.end(); recordIter++)
			{
				if (recordIter->clusterEnd == 0)
				{
					clusters.back().push_back(recordIter->GetReadRecord());
				}
			}
		}
	}

	vector<double> weights;

	if (!weightsFilename.empty())
	{
		cout << "Reading weights" << endl;

		ifstream weightsFile(weightsFilename.c_str());
		CheckFile(weightsFile, weightsFilename);

		double weight;
		while (weightsFile >> weight)
		{
			weights.push_back(weight);
		}

		if (clusters.size() != weights.size())
		{
			cerr << "Error: read " << clusters.size() << " clusters and " << weights.size() << " weights" << endl;
			exit(1);
		}
	}
	else
	{
		cout << "Setting weights to 1" << endl;

		weights = vector<double>(clusters.size(), 1.0);
	}

	cout << "Calculating set cover solution" << endl;

	vector<int> solution;
	SetCover(clusters, weights, solution);
	
	cout << "Assigning fragments to solution sets" << endl;
	
	vector<vector<ReadRecord> > assignment;
	AssignInOrder(clusters, solution, assignment);

	DebugCheck(ids.size() == assignment.size());
	
	cout << "Writing out assignments" << endl;

	{
		ofstream assignmentsFile(assignmentsFilename.c_str());
		CheckFile(assignmentsFile, assignmentsFilename);

		ifstream clustersFile(clustersFilename.c_str());
		CheckFile(clustersFile, clustersFilename);

		GroupedRecordsStream<ClusterMemberRecord> memberStream(clustersFile);

		int fileIdx = 0;

		vector<ClusterMemberRecord> clusterRecords;
		while (memberStream.Next(clusterRecords, ClusterEqual<ClusterMemberRecord>))
		{
			DebugCheck(clusterRecords.size() > 0);
			DebugCheck(clusterRecords.front().clusterID == ids[fileIdx]);

			unordered_set<ReadRecord> assigned(clusters[fileIdx].begin(), clusters[fileIdx].end());

			for (vector<ClusterMemberRecord>::const_iterator recordIter = clusterRecords.begin(); recordIter != clusterRecords.end(); recordIter++)
			{
				if (assigned.find(recordIter->GetReadRecord()) != assigned.end())
				{
					assignmentsFile << (*recordIter);
				}
			}

			fileIdx++;
		}
	}
}

