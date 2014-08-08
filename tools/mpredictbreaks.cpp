#include "AlignmentStream.h"
#include "DebugCheck.h"
#include "Indexer.h"
#include "DiscordantAlignments.h"
#include "MatePairDelauny.h"
#include "MatePairGibbs.h"
#include "Parsers.h"
#include "Sequences.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


struct SplitDefinition
{
	int position[2];
	int numInserted;
};

inline bool operator==(const SplitDefinition& split1, const SplitDefinition& split2)
{
	return split1.position[0] == split2.position[0] &&
	       split1.position[1] == split2.position[1] &&
	       split1.numInserted == split2.numInserted;
}

inline size_t hash_value(const SplitDefinition& split)
{
	size_t seed = 0;
	hash_combine(seed, split.position[0]);
	hash_combine(seed, split.position[1]);
	hash_combine(seed, split.numInserted);
	return seed;
}

bool SplitScoreLessThan(const pair<SplitDefinition,int>& a, const pair<SplitDefinition,int>& b)
{
	return a.second < b.second;
}

void Complement(char& nucleotide)
{
	switch (nucleotide)
	{
		case 'A': nucleotide = 'T'; break;
		case 'C': nucleotide = 'G'; break;
		case 'T': nucleotide = 'A'; break;
		case 'G': nucleotide = 'C'; break;
		case 'a': nucleotide = 't'; break;
		case 'c': nucleotide = 'g'; break;
		case 't': nucleotide = 'a'; break;
		case 'g': nucleotide = 'c'; break;
	}
}

int CalculateOffset(const string& strand, int offset)
{
	int dir = (strand == "+") ? 1 : -1;
	return offset * dir;
}

int CalculateForwardHomology(int maxOffset, const Sequences& sequences, const string (&chromosome)[2], const string (&strand)[2], const int (&position)[2], bool flip=false)
{
	int idx1 = (flip) ? 1 : 0;
	int idx2 = 1 - idx1;

	const char* seqPtr1 = sequences.Get(chromosome[idx1], position[idx1]);
	const char* seqPtr2 = sequences.Get(chromosome[idx2], position[idx2]);

	int homology = 0;
	for (int offset = 1; offset <= maxOffset; offset++)
	{
		char nt1 = *(seqPtr1 + CalculateOffset(strand[idx1], offset));
		char nt2 = *(seqPtr2 + CalculateOffset(strand[idx2], 1 - offset));
		
		if (strand[idx1] != "+")
		{
			Complement(nt1);
		}
		
		if (strand[idx2] != "-")
		{
			Complement(nt2);
		}
		
		if (nt1 != nt2)
		{
			break;
		}
		
		homology = offset;
	}
	
	return homology;
}

void WriteBreakpoint(ostream& out, int clusterID, const string (&chromosome)[2], const string (&strand)[2], const SplitDefinition& breakpoint, int splitCount, int offset1, int offset2)
{
	BreakpointRecord record;
	
	record.clusterID = clusterID;
	record.chromosome[0] = chromosome[0];
	record.strand[0] = strand[0];
	record.position[0] = breakpoint.position[0] + CalculateOffset(strand[0], offset1);
	record.chromosome[1] = chromosome[1];
	record.strand[1] = strand[1];
	record.position[1] = breakpoint.position[1] + CalculateOffset(strand[1], offset2);
	record.numInserted = breakpoint.numInserted;
	record.splitCount = splitCount;

	out << record;
}

struct ClusterSplitInfo
{
	string chromosome[2];
	string strand[2];
	unordered_map<SplitDefinition,int> count;
	unordered_map<SplitDefinition,int> score;
};

int main(int argc, char* argv[])
{
	string alignmentsFilename;
	string clustersFilename;
	string breakpointsFilename;
	string referenceFasta;
	
	try
	{
		TCLAP::CmdLine cmd("Breakpoint Prediction Tool");
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","alignments","Split Alignments",true,"","string",cmd);
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> breakpointsFilenameArg("b","breakpoints","Output Breakpoints Filename",true,"","integer",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		alignmentsFilename = alignmentsFilenameArg.getValue();
		clustersFilename = clustersFilenameArg.getValue();
		breakpointsFilename = breakpointsFilenameArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	const int cMaxBreakpointHomology = 25;
	
	cerr << "Reading clusters" << endl;

	unordered_map<AlignmentPairKey,unordered_set<pair<int,bool> > > memberships;
	
	ifstream clustersFile(clustersFilename.c_str());
	CheckFile(clustersFile, clustersFilename);

	GroupedRecordsStream<ClusterMemberRecord> memberStream(clustersFile);

	int numClusters = 0;

	vector<ClusterMemberRecord> clusterRecords;
	while (memberStream.Next(clusterRecords, ClusterReadEqual<ClusterMemberRecord>))
	{
		DebugCheck(clusterRecords.size() == 2);
		DebugCheck(clusterRecords[0].clusterID == clusterRecords[1].clusterID);
		DebugCheck(clusterRecords[0].clusterEnd != clusterRecords[1].clusterEnd);
		DebugCheck(clusterRecords[0].libID == clusterRecords[1].libID);
		DebugCheck(clusterRecords[0].readID == clusterRecords[1].readID);
		DebugCheck(clusterRecords[0].readEnd != clusterRecords[1].readEnd);

		AlignmentPairKey alignPairKey;

		alignPairKey.libID = clusterRecords[0].libID;
		alignPairKey.readID = clusterRecords[0].readID;
		for (int idx = 0; idx < 2; idx++)
		{
			alignPairKey.alignID[clusterRecords[idx].readEnd] = clusterRecords[idx].alignID;
		}

		bool flip = (clusterRecords[0].readEnd != clusterRecords[0].clusterEnd);

		memberships[alignPairKey].insert(make_pair(clusterRecords[0].clusterID, flip));

		numClusters++;
	}

	cerr << "Read " << numClusters << " clusters" << endl;

	cerr << "Reading split alignments" << endl;

	ifstream alignmentsFile(alignmentsFilename.c_str());
	CheckFile(alignmentsFile, alignmentsFilename);

	unordered_map<int,ClusterSplitInfo> clusterSplitInfo;
	
	int numAlignments = 0;

	SplitAlignmentRecord splitRecord;
	while (alignmentsFile >> splitRecord)
	{
		AlignmentPairKey alignPairKey = splitRecord.GetAlignmentPairKey();
		
		unordered_map<AlignmentPairKey,unordered_set<pair<int,bool> > >::const_iterator clusterSetIter = memberships.find(alignPairKey);

		if (clusterSetIter != memberships.end())
		{
			for (unordered_set<pair<int,bool> >::const_iterator clusterIter = clusterSetIter->second.begin(); clusterIter != clusterSetIter->second.end(); clusterIter++)
			{
				int clusterID = clusterIter->first;
				bool flip = clusterIter->second;
				
				SplitDefinition split;
				
				split.position[0] = splitRecord.position[0];
				split.position[1] = splitRecord.position[1];
				split.numInserted = (int)splitRecord.inserted.size();

				for (int readEnd = 0; readEnd < 2; readEnd++)
				{
					int clusterEnd = readEnd;
					if (flip)
					{
						clusterEnd = 1 - readEnd;
					}

					clusterSplitInfo[clusterID].chromosome[clusterEnd] = splitRecord.chromosome[readEnd];
					clusterSplitInfo[clusterID].strand[clusterEnd] = splitRecord.strand[readEnd];

					split.position[clusterEnd] = splitRecord.position[readEnd];
				}
				
				clusterSplitInfo[clusterID].count.insert(make_pair(split, 0)).first->second++;
				clusterSplitInfo[clusterID].score.insert(make_pair(split, 0)).first->second += splitRecord.score;
			}
		}
		
		numAlignments++;
	}
	
	cerr << "Read " << numAlignments << " split alignments" << endl;

	cerr << "Reading reference fasta" << endl;
	
	Sequences referenceSequences(1000);
	referenceSequences.Read(referenceFasta);
	
	cerr << "Calculating breakpoint homology and outputting breakpoints" << endl;
	
	ofstream breakpointsFile(breakpointsFilename.c_str());
	CheckFile(breakpointsFile, breakpointsFilename);
	
	for (unordered_map<int,ClusterSplitInfo>::const_iterator splitInfoIter = clusterSplitInfo.begin(); splitInfoIter != clusterSplitInfo.end(); splitInfoIter++)
	{
		int clusterID = splitInfoIter->first;
		const ClusterSplitInfo& clusterSplitInfo = splitInfoIter->second;
		
		const SplitDefinition& breakpoint = max_element(clusterSplitInfo.score.begin(), clusterSplitInfo.score.end(), SplitScoreLessThan)->first;
		int splitCount = clusterSplitInfo.count.find(breakpoint)->second;
		
		int maxOffset1 = CalculateForwardHomology(cMaxBreakpointHomology, referenceSequences, clusterSplitInfo.chromosome, clusterSplitInfo.strand, breakpoint.position, false);
		int maxOffset2 = CalculateForwardHomology(cMaxBreakpointHomology, referenceSequences, clusterSplitInfo.chromosome, clusterSplitInfo.strand, breakpoint.position, true);
		
		if (maxOffset1 + maxOffset2 + 1 > cMaxBreakpointHomology)
		{
			continue;
		}
		
		for (int offset = -maxOffset2; offset <= maxOffset1; offset++)
		{
			WriteBreakpoint(breakpointsFile, clusterID, clusterSplitInfo.chromosome, clusterSplitInfo.strand, breakpoint, splitCount, offset, -offset);
		}
	}
}

