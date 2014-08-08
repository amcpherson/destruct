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
	ClusterMemberRecord record;

	record.clusterID = clusterID;
	record.clusterEnd = clusterEnd;
	record.libID = readInfo.libID;
	record.readID = readInfo.readID;
	record.readEnd = readInfo.readEnd;
	record.alignID = readInfo.alignID;

	out << record;
}

void OutputBreakend(ostream& out, int clusterID, int clusterEnd, const string& chromosome, const string& strand, int position)
{
	BreakendRecord record;

	record.clusterID = clusterID;
	record.clusterEnd = clusterEnd;
	record.chromosome = chromosome;
	record.strand = strand;
	record.position = abs(position);

	out << record;
}

int main(int argc, char* argv[])
{
	string alignmentsFilename;
	string libStatsFilename;
	string clustersFilename;
	string breakendsFilename;
	int minClusterSize;
	int maxFragmentLength;
	string chromPair;
	string exclChromPairs;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Pair Clustering Tool");
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","alignments","Paired End Alignments",true,"","string",cmd);
		TCLAP::ValueArg<string> libStatsFilenameArg("s","stats","Library Read Stats",true,"","string",cmd);
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> breakendsFilenameArg("b","breakends","Breakends Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> minClusterSizeArg("","clustmin","Minimum Cluster Size",true,-1,"integer",cmd);
		TCLAP::ValueArg<int> maxFragmentLengthArg("","fragmax","Maximum Fragment Length",true,-1,"integer",cmd);
		TCLAP::ValueArg<string> chromPairArg("","inclchrompair","Include Chromosome Pair (comma separated)",false,"","string",cmd);
		TCLAP::ValueArg<string> exclChromPairsArg("","exclchrompairs","Exclude Chromosome Pairs from Set (comma separated)",false,"","string",cmd);
		cmd.parse(argc,argv);
		
		alignmentsFilename = alignmentsFilenameArg.getValue();
		libStatsFilename = libStatsFilenameArg.getValue();
		clustersFilename = clustersFilenameArg.getValue();
		breakendsFilename = breakendsFilenameArg.getValue();
		minClusterSize = minClusterSizeArg.getValue();
		maxFragmentLength = maxFragmentLengthArg.getValue();
		chromPair = chromPairArg.getValue();
		exclChromPairs = exclChromPairsArg.getValue();
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

	if (!chromPair.empty())
	{
		vector<string> chromPairFields;
		split(chromPairFields, chromPair, is_any_of(","));

		if (chromPairFields.size() != 2)
		{
			cerr << "Error: Require 2 chromosomes separated by comma for chrompair argument" << endl;
			exit(1);
		}

		cout << "Restricting analysis to alignments connecting chromosomes ";
		cout << chromPairFields[0] << " and " << chromPairFields[1] << endl;

		discordantAlignments.SetIncludedChromosomePair(chromPairFields[0], chromPairFields[1]);
	}

	if (!exclChromPairs.empty())
	{
		vector<string> exclChromPairsFields;
		split(exclChromPairsFields, exclChromPairs, is_any_of(","));

		cout << "Excluding analysis of alignments connecting pairs of chromosomes from ";
		for (vector<string>::const_iterator chromIter = exclChromPairsFields.begin(); chromIter != exclChromPairsFields.end(); chromIter++)
		{
			cout << *chromIter << " ";
		}
		cout << endl;

		discordantAlignments.SetExcludedChromosomePairs(exclChromPairsFields);
	}
	
	ifstream alignmentsFile(alignmentsFilename.c_str());
	CheckFile(alignmentsFile, alignmentsFilename);

	GroupedRecordsStream<SpanningAlignmentRecord> alignmentsStream(alignmentsFile);

	vector<SpanningAlignmentRecord> readAlignments;
	while (alignmentsStream.Next(readAlignments, ReadEqual<SpanningAlignmentRecord>))
	{
		discordantAlignments.AddFragmentAlignments(readAlignments);
	}

	cout << "Read " << discordantAlignments.GetAlignmentCount() << " alignments for " << discordantAlignments.GetFragmentCount() << " fragments" << endl;
	
	ofstream clustersFile(clustersFilename.c_str());
	CheckFile(clustersFile, clustersFilename);

	ofstream breakendsFile(breakendsFilename.c_str());
	CheckFile(breakendsFile, breakendsFilename);
	
	cout << "Creating clusters" << endl;
	
	int clusterID = 0;

	vector<pair<uint32_t,uint32_t> > chrStrIdxPairs = discordantAlignments.GetChrStrIdxPairs();

	for (vector<pair<uint32_t,uint32_t> >::const_iterator chrStrIdxPairIter = chrStrIdxPairs.begin(); chrStrIdxPairIter != chrStrIdxPairs.end(); chrStrIdxPairIter++)
	{
		pair<string,string> chromosomes = discordantAlignments.GetChromosomePair(*chrStrIdxPairIter);
		pair<string,string> strands = discordantAlignments.GetStrandPair(*chrStrIdxPairIter);
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

			OutputBreakend(breakendsFile, clusterID, 0, chromosomes.first, strands.first, (int)matePairs[0].x);
			OutputBreakend(breakendsFile, clusterID, 1, chromosomes.second, strands.second, (int)matePairs[0].y);

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

				int clusterX = -1;
				int clusterY = -1;
				
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

					if (clusterX == -1 && clusterY == -1)
					{
						clusterX = (int)matePairs[alignPairIndex].x;
						clusterY = (int)matePairs[alignPairIndex].y;
					}
					else
					{
						clusterX = max(clusterX, (int)matePairs[alignPairIndex].x);
						clusterY = max(clusterY, (int)matePairs[alignPairIndex].y);
					}
				}

				OutputBreakend(breakendsFile, clusterID, 0, chromosomes.first, strands.first, clusterX);
				OutputBreakend(breakendsFile, clusterID, 1, chromosomes.second, strands.second, clusterY);

				clusterID++;
			}
		}
	}
	
	clustersFile.close();
	
	cout << "Created " << clusterID << " clusters" << endl;
}

