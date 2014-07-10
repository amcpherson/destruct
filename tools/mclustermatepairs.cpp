/*
 *  mclustermatepairs.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "DebugCheck.h"
#include "Parsers.h"
#include "DiscordantAlignments.h"
#include "MatePairDelauny.h"
#include "MatePairGibbs.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


void ReadLibStats(const string& libStatsFilename, vector<double>& fragmentLengthMeans, vector<double>& fragmentLengthStdDevs)
{
	ifstream libStatsFile(libStatsFilename.c_str());
	CheckFile(libStatsFile, libStatsFilename);
	
	vector<double> libFragmentLengthMeans;
	vector<double> libFragmentLengthStdDevs;
	
	StringVec fields;
	while (ReadTSV(libStatsFile, fields))
	{
		int libId = SAFEPARSE(int, fields[0]);

		if (libId >= fragmentLengthMeans.size())
		{
			fragmentLengthMeans.resize(libId + 1);
			fragmentLengthStdDevs.resize(libId + 1);
		}

		fragmentLengthMeans[libId] = SAFEPARSE(double, fields[1]);
		fragmentLengthStdDevs[libId] = SAFEPARSE(double, fields[2]);
	}
}

void OutputClusterMember(ostream& out, int clusterID, int clusterEnd, const ReadInfo& readInfo)
{
	ClusterMemberRecord memberRecord;

	memberRecord.clusterID = clusterID;
	memberRecord.clusterEnd = clusterEnd;
	memberRecord.libID = readInfo.libID;
	memberRecord.readID = readInfo.readID;
	memberRecord.readEnd = readInfo.readEnd;
	memberRecord.alignID = readInfo.alignID;

	out << memberRecord;
}

int main(int argc, char* argv[])
{
	string alignmentsFilename;
	string libStatsFilename;
	string clustersFilename;
	int minClusterSize;
	int maxFragmentLength;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Pair Clustering Tool");
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","alignments","Paired End Alignments",true,"","string",cmd);
		TCLAP::ValueArg<string> libStatsFilenameArg("s","stats","Library Read Stats",true,"","string",cmd);
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> minClusterSizeArg("","clustmin","Minimum Cluster Size",true,-1,"integer",cmd);
		TCLAP::ValueArg<int> maxFragmentLengthArg("","fragmax","Maximum Fragment Length",true,-1,"integer",cmd);
		cmd.parse(argc,argv);
		
		alignmentsFilename = alignmentsFilenameArg.getValue();
		libStatsFilename = libStatsFilenameArg.getValue();
		clustersFilename = clustersFilenameArg.getValue();
		minClusterSize = minClusterSizeArg.getValue();
		maxFragmentLength = maxFragmentLengthArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	vector<double> fragmentLengthMeans;
	vector<double> fragmentLengthStdDevs;
	ReadLibStats(libStatsFilename, fragmentLengthMeans, fragmentLengthStdDevs);
	
	DiscordantAlignments discordantAlignments(fragmentLengthMeans, fragmentLengthStdDevs, maxFragmentLength);
	
	ifstream alignmentsFile(alignmentsFilename.c_str());
	CheckFile(alignmentsFile, alignmentsFilename);

	GroupedRecordsStream<SpanningAlignmentRecord> alignmentsStream(alignmentsFile);

	vector<SpanningAlignmentRecord> readAlignments;
	while (alignmentsStream.Next(readAlignments, ReadEqual<SpanningAlignmentRecord>))
	{
		discordantAlignments.AddFragmentAlignments(readAlignments);
	}

	cerr << "Read " << discordantAlignments.GetAlignmentCount() << " alignments for " << discordantAlignments.GetFragmentCount() << " fragments" << endl;
	
	// Open output clusters file
	ofstream clustersFile(clustersFilename.c_str());
	if (!clustersFile)
	{
		cerr << "Error: unable to write to clusters file" << endl;		
		exit(1);
	}
	
	cout << "Creating clusters" << endl;
	
	int clusterID = 0;

	vector<pair<uint32_t,uint32_t> > chrStrIdxPairs = discordantAlignments.GetChrStrIdxPairs();

	for (vector<pair<uint32_t,uint32_t> >::const_iterator chrStrIdxPairIter = chrStrIdxPairs.begin(); chrStrIdxPairIter != chrStrIdxPairs.end(); chrStrIdxPairIter++)
	{
		vector<MatePair> matePairs = discordantAlignments.CreateMatePairs(*chrStrIdxPairIter);
		vector<pair<ReadInfo,ReadInfo> > readInfos = discordantAlignments.CreateReadInfos(*chrStrIdxPairIter);

		if (matePairs.size() < minClusterSize)
		{
			continue;
		}
		
		if (matePairs.size() == 1)
		{
			OutputClusterMember(clustersFile, clusterID, 0, readInfos.front().first);
			OutputClusterMember(clustersFile, clusterID, 1, readInfos.front().second);
			clusterID++;
			continue;
		}
		
		MatePairDelauny delaunyClusterer;
		
		IntegerTable delaunyClusters;
		delaunyClusterer.DoClustering(matePairs, delaunyClusters);
		
		for (IntegerTableConstIter delaunyClusterIter = delaunyClusters.begin(); delaunyClusterIter != delaunyClusters.end(); delaunyClusterIter++)
		{
			MatePairVec delaunyMatePairs;
			IntegerVec alignPairIndices;
			for (IntegerVecConstIter elementIter = delaunyClusterIter->begin(); elementIter != delaunyClusterIter->end(); elementIter++)
			{
				delaunyMatePairs.push_back(matePairs[*elementIter]);
				alignPairIndices.push_back(*elementIter);
			}
			
			MatePairGibbs gibbsClusterer;
			
			IntegerTable gibbsClusters;
			gibbsClusterer.DoClustering(delaunyMatePairs, gibbsClusters);
			
			for (int clusterIndex = 0; clusterIndex < gibbsClusters.size(); clusterIndex++)
			{
				const IntegerVec& cluster = gibbsClusters[clusterIndex];
				
				unordered_set<int> clusterFragmentIndices;
				for (int elementIndex = 0; elementIndex < cluster.size(); elementIndex++)
				{
					int alignPairIndex = alignPairIndices[cluster[elementIndex]];

					// For fragments with multiple alignments supporting the same cluster, select first alignment
					if (!clusterFragmentIndices.insert(readInfos[alignPairIndex].first.readID).second)
					{
						continue;
					}
					
					OutputClusterMember(clustersFile, clusterID, 0, readInfos[alignPairIndex].first);
					OutputClusterMember(clustersFile, clusterID, 1, readInfos[alignPairIndex].second);
				}
				
				clusterID++;
			}
		}
	}
	
	clustersFile.close();
	
	cout << "Created " << clusterID << " clusters" << endl;
}

