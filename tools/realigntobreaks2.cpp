

#include "Common.h"
#include "DebugCheck.h"
#include "AlignRead.h"
#include "AlignmentRecord.h"
#include "ReadStream.h"
#include "Sequences.h"

#include <boost/unordered_map.hpp>
#include <tclap/CmdLine.h>

using namespace boost;
using namespace std;


int CalculateRealignedScore(SimpleAligner& aligner, PreppedReads& reads, int readID, int readEnd, const string& breakpointSequence, const string& strand, int position)
{
	reads.SetCurrentRead(readID);

	const char* refPtr = &breakpointSequence[position];

	if (strand == "+")
	{
		if (position < 16 || position > breakpointSequence.size() / 2)
		{
			return 0;
		}
		
		return aligner.AlignBandedSSE2BW7ScoreFwd(refPtr, reads.StartPtr(readEnd, PlusStrand), reads.EndPtr(readEnd, PlusStrand));
	}
	else
	{
		if (position < breakpointSequence.size() / 2 || position > breakpointSequence.size() - 16)
		{
			return 0;
		}
		
		return aligner.AlignBandedSSE2BW7ScoreRev(refPtr, reads.StartPtr(readEnd, MinusStrand), reads.EndPtr(readEnd, MinusStrand));
	}
}

int main(int argc, char* argv[])
{
	int matchScore;
	int misMatchScore;
	int gapScore;
	int maxFragmentLength;
	string referenceFasta;
	string readSeqsFilename;
	string spanningFilename;
	string clustersFilename;
	string breakpointsFilename;
	string realignmentsFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Realignment to Breakpoints Tool");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> maxFragmentLengthArg("","flmax","Maximum Fragment Length",true,0,"int",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqsFilenameArg("s","seqs","Read Sequences Fastq",true,"","string",cmd);
		TCLAP::ValueArg<string> spanningFilenameArg("","span","Spanning Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> breakpointsFilenameArg("b","breakpoints","Breakpoints Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> realignmentsFilenameArg("","realignments","Realignment Scores Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		maxFragmentLength = maxFragmentLengthArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		readSeqsFilename = readSeqsFilenameArg.getValue();
		spanningFilename = spanningFilenameArg.getValue();
		clustersFilename = clustersFilenameArg.getValue();
		breakpointsFilename = breakpointsFilenameArg.getValue();
		realignmentsFilename = realignmentsFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}

	SimpleAligner aligner(matchScore, misMatchScore, gapScore);

	cerr << "Reading reference fasta" << endl;
	
	Sequences referenceSequences(1000);
	referenceSequences.Read(referenceFasta);
	
	cerr << "Reading fastq sequences" << endl;
	
	ifstream readSeqsFile(readSeqsFilename.c_str());
	CheckFile(readSeqsFile, readSeqsFilename);
	
	FastqReadStream readSeqsStream(readSeqsFile);
	
	PreppedReads preppedReads;
	preppedReads.Prep(readSeqsStream);
	
	cerr << "Reading spanning alignments" << endl;

	unordered_map<AlignmentKey,SpanningAlignmentRecord> spanningAlignments;
	{
		ifstream spanningFile(spanningFilename.c_str());
		CheckFile(spanningFile, spanningFilename);

		SpanningAlignmentRecord spanningRecord;
		while (spanningFile >> spanningRecord)
		{
			spanningAlignments[spanningRecord.GetAlignmentKey()] = spanningRecord;
		}
	}

	cerr << "Reading cluster memberships" << endl;

	unordered_map<int,vector<ClusterMemberRecord> > clusters;
	{
		ifstream clustersFile(clustersFilename.c_str());
		CheckFile(clustersFile, clustersFilename);

		ClusterMemberRecord memberRecord;
		while (clustersFile >> memberRecord)
		{
			// Ignore memberships of unknown alignments
			if (spanningAlignments.find(memberRecord.GetAlignmentKey()) == spanningAlignments.end())
			{
				continue;
			}

			clusters[memberRecord.clusterID].push_back(memberRecord);
		}
	}

	cerr << "Iterating breakpoints" << endl;

	ifstream breakpointsFile(breakpointsFilename.c_str());
	CheckFile(breakpointsFile, breakpointsFilename);

	ofstream realignmentsFile(realignmentsFilename.c_str());
	CheckFile(realignmentsFile, realignmentsFilename);

	BreakpointRecord breakpointRecord;
	while (breakpointsFile >> breakpointRecord)
	{
		// Create breakend sequences of specific lengths, and record the start and end of the region from which the
		// sequence originated in the reference genome.
		string breakendSequence[2];
		int breakendStart[2];
		int breakendEnd[2];
		for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
		{
			if (breakpointRecord.strand[clusterEnd] == "+")
			{
				breakendStart[clusterEnd] = breakpointRecord.position[clusterEnd] - maxFragmentLength + 1;
				breakendEnd[clusterEnd] = breakpointRecord.position[clusterEnd];
			}
			else
			{
				breakendStart[clusterEnd] = breakpointRecord.position[clusterEnd];
				breakendEnd[clusterEnd] = breakpointRecord.position[clusterEnd] + maxFragmentLength - 1;
			}

			referenceSequences.Get(breakpointRecord.chromosome[clusterEnd],
			                       breakendStart[clusterEnd],
			                       breakendEnd[clusterEnd],
			                       breakendSequence[clusterEnd]);
		}

		// Create a breakpoint sequence for each cluster end with the sequence of the reference maintained (not reverse
		// complemented) for that breakend.  Also calculate an offset which will be used for calculating the 0-based
		// positions of alignments in the breakpoint sequence.
		string breakpointSequence[2];
		int breakpointOffset[2];
		for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
		{
			string otherBreakendSequence = breakendSequence[1-clusterEnd];

			if (breakpointRecord.strand[1-clusterEnd] != breakpointRecord.strand[clusterEnd])
			{
				ReverseComplement(otherBreakendSequence);
			}

			breakpointOffset[clusterEnd] = -breakendStart[clusterEnd];

			if (breakpointRecord.strand[clusterEnd] == "+")
			{
				breakpointSequence[clusterEnd] = breakendSequence[clusterEnd] + otherBreakendSequence;
			}
			else
			{
				breakpointSequence[clusterEnd] = otherBreakendSequence + breakendSequence[clusterEnd];
				breakpointOffset[clusterEnd] += maxFragmentLength;
			}
		}

		// Iterate cluster reads and their alignments
		unordered_map<int,vector<ClusterMemberRecord> >::const_iterator clusterIter = clusters.find(breakpointRecord.clusterID);

		if (clusterIter == clusters.end())
		{
			cerr << "Error: No matching cluster for breakpoint " << breakpointRecord.clusterID << endl;
			exit(1);
		}

		const vector<ClusterMemberRecord>& memberships = clusterIter->second;

		for (vector<ClusterMemberRecord>::const_iterator memberIter = memberships.begin(); memberIter != memberships.end(); memberIter++)
		{
			const ClusterMemberRecord& memberRecord = *memberIter;

			unordered_map<AlignmentKey,SpanningAlignmentRecord>::const_iterator alignIter = spanningAlignments.find(memberRecord.GetAlignmentKey());
			DebugCheck(alignIter != spanningAlignments.end());

			const SpanningAlignmentRecord& spanningRecord = alignIter->second;

			// Calculate position of alignment in breakpoint sequence
			int adjustedPosition = spanningRecord.GetOuterPosition() + breakpointOffset[memberRecord.clusterEnd];

			int score = CalculateRealignedScore(aligner, preppedReads, memberRecord.readID, memberRecord.readEnd, breakpointSequence[memberRecord.clusterEnd], spanningRecord.strand, adjustedPosition);

			ClusterMemberScoreRecord scoreRecord;

			scoreRecord.clusterID = memberRecord.clusterID;
			scoreRecord.clusterEnd = memberRecord.clusterEnd;
			scoreRecord.libID = memberRecord.libID;
			scoreRecord.readID = memberRecord.readID;
			scoreRecord.readEnd = memberRecord.readEnd;
			scoreRecord.alignID = memberRecord.alignID;
			scoreRecord.score = score;

			realignmentsFile << scoreRecord;
		}
	}
}


