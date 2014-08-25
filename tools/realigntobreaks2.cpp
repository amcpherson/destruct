

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
		for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
		{
			string selfbreakendSequence = breakendSequence[clusterEnd];
			string mateBreakendSequence = breakendSequence[1-clusterEnd];

			if (breakpointRecord.strand[clusterEnd] == "-")
			{
				ReverseComplement(selfbreakendSequence);
			}

			if (breakpointRecord.strand[1-clusterEnd] == "+")
			{
				ReverseComplement(mateBreakendSequence);
			}

			string insertedSequence = breakpointRecord.inserted;

			if (clusterEnd == 1)
			{
				ReverseComplement(insertedSequence);
			}

			breakpointSequence[clusterEnd] = selfbreakendSequence + insertedSequence + mateBreakendSequence;
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

			// Calculate the length between the breakend and the start of the read
			int templateLength = abs(spanningRecord.GetOuterPosition() - breakpointRecord.position[memberRecord.clusterEnd]) + 1;

			// Calculate position of alignment in breakpoint sequence
			int adjustedPosition = maxFragmentLength - templateLength;

			// Set current read in prepped read set
			preppedReads.SetCurrentRead(memberRecord.readID);

			int score = 0;
			if (adjustedPosition >= 16 && adjustedPosition <= breakpointSequence[memberRecord.clusterEnd].size() / 2)
			{
				// Pointer to first position to which the read should align
				const char* refPtr = &breakpointSequence[memberRecord.clusterEnd][adjustedPosition];

				// Align to forward strand of breakpoint sequence
				score = aligner.AlignBandedSSE2BW7ScoreFwd(refPtr, preppedReads.StartPtr(memberRecord.readEnd, PlusStrand), preppedReads.EndPtr(memberRecord.readEnd, PlusStrand));
			}

			BreakAlignScoreRecord scoreRecord;

			scoreRecord.clusterID = breakpointRecord.clusterID;
			scoreRecord.predictionID = breakpointRecord.predictionID;
			scoreRecord.clusterEnd = memberRecord.clusterEnd;
			scoreRecord.libID = memberRecord.libID;
			scoreRecord.readID = memberRecord.readID;
			scoreRecord.readEnd = memberRecord.readEnd;
			scoreRecord.alignID = memberRecord.alignID;
			scoreRecord.alignedLength = preppedReads.ReadLength(memberRecord.readEnd);
			scoreRecord.templateLength = templateLength;
			scoreRecord.score = score;

			realignmentsFile << scoreRecord;
		}
	}
}


