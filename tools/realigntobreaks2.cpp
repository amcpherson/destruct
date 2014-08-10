

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
	string referenceFasta;
	string readSeqsFilename;
	string spanningFilename;
	string clustersFilename;
	string breakpointsFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Realignment to Breakpoints Tool");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqsFilenameArg("s","seqs","Read Sequences Fastq",true,"","string",cmd);
		TCLAP::ValueArg<string> spanningFilenameArg("","span","Spanning Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> breakpointsFilenameArg("b","breakpoints","Output Breakpoints Filename",true,"","integer",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		readSeqsFilename = readSeqsFilenameArg.getValue();
		spanningFilename = spanningFilenameArg.getValue();
		clustersFilename = clustersFilenameArg.getValue();
		breakpointsFilename = breakpointsFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}

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

	ifstream spanningFile(spanningFilename.c_str());
	CheckFile(spanningFile, spanningFilename);

	unordered_map<AlignmentKey,SpanningAlignmentRecord> spanningAlignments;

	SpanningAlignmentRecord spanningRecord;
	while (spanningFile >> spanningRecord)
	{
		spanningAlignments[spanningRecord.GetAlignmentKey()] = spanningRecord;
	}

	cerr << "Reading cluster memberships" << endl;

	ifstream clustersFile(clustersFilename.c_str());
	CheckFile(clustersFile, clustersFilename);

	unordered_map<int,vector<ClusterMemberRecord> > clusters;

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

	cerr << "Iterating breakpoints" << endl;

	ifstream breakpointsFile(breakpointsFilename.c_str());
	CheckFile(breakpointsFile, breakpointsFilename);

	BreakpointRecord breakpointRecord;
	while (breakpointsFile >> breakpointRecord)
	{
		// Reconstruct breakpoint sequence
		

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
			unordered_map<AlignmentKey,SpanningAlignmentRecord>::const_iterator alignIter = spanningAlignments.find(memberRecord.GetAlignmentKey());
			DebugCheck(alignIter != spanningAlignments.end());

			const SpanningAlignmentRecord& spanningRecord = alignIter->second;

			// Calculate position of alignment in breakpoint sequence

			preppedReads.SetCurrentRead(spanningRecord.readID);

		}
	}
}


