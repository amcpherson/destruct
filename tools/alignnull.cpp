/*
 *  alignnull.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "AlignmentStream.h"
#include "ReadStream.h"
#include "Sequences.h"
#include "SimpleAligner.h"
#include "AlignRead.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


bool Coincident(const Location& location, const RawAlignment& alignment)
{
	if (location.refName != alignment.reference)
	{
		return false;
	}
	
	if (location.strand != alignment.strand)
	{
		return false;
	}
	
	if (location.start > alignment.region.end || location.end < alignment.region.start)
	{
		return false;
	}
	
	return true;
}

int main(int argc, char* argv[])
{
	int matchScore;
	int misMatchScore;
	int gapScore;
	string referenceFasta;
	string readSeqsFilename;
	string readSeqs2Filename;
	string alignmentsFilename;
	string concordantFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Realignment Null Distribution");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqsFilenameArg("s","seq","Read Sequences",true,"","string",cmd);
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Sam Alignments",true,"","string",cmd);
		TCLAP::ValueArg<string> concordantFilenameArg("c","concordant","Concordant Sam Alignments",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		readSeqsFilename = readSeqsFilenameArg.getValue();
		alignmentsFilename = alignmentsFilenameArg.getValue();
		concordantFilename = concordantFilenameArg.getValue();
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
	
	cerr << "Reading concordant alignments" << endl;
	
	unordered_map<int,Location> concordantAlignments;
	{
		SamAlignmentStream alignmentStream(concordantFilename);
		
		RawAlignment alignment;
		while (alignmentStream.GetNextAlignment(alignment))
		{
			ReadID readID;
			readID.fragmentIndex = SAFEPARSE(int, alignment.fragment);
			readID.readEnd = alignment.readEnd;
			
			concordantAlignments[readID.id].refName = alignment.reference;
			concordantAlignments[readID.id].strand = alignment.strand;
			concordantAlignments[readID.id].start = alignment.region.start;
			concordantAlignments[readID.id].end = alignment.region.end;
		}
	}
	
	cerr << "Realigning" << endl;
	
	SimpleAligner aligner(matchScore, misMatchScore, gapScore);
	
	SamAlignmentStream alignmentStream(alignmentsFilename);
	FragmentAlignmentStream fragmentAlignmentStream(&alignmentStream);
	
	RawAlignmentVec alignments;
	while (fragmentAlignmentStream.GetNextAlignments(alignments))
	{
		int fragmentIndex = SAFEPARSE(int, alignments.front().fragment);
		
		preppedReads.SetCurrentRead(fragmentIndex);
		
		int bestScore[2] = {numeric_limits<int>::min(), numeric_limits<int>::min()};
		int bestAlignmentIndex[2] = {-1, -1};
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			ReadID readID;
			readID.fragmentIndex = fragmentIndex;
			readID.readEnd = alignment.readEnd;
			
			unordered_map<int,Location>::const_iterator concordantIter = concordantAlignments.find(readID.id);
			if (concordantIter != concordantAlignments.end() && Coincident(concordantIter->second, alignment))
			{
				continue;
			}
			
			int score = AlignSelfScoreSSE(aligner, alignment, referenceSequences, preppedReads);
			
			if (score > bestScore[alignment.readEnd])
			{
				bestScore[alignment.readEnd] = score;
				bestAlignmentIndex[alignment.readEnd] = alignmentIndex;
			}
		}
		
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			if (bestAlignmentIndex[readEnd] == -1)
			{
				continue;
			}
			
			const RawAlignment& alignment = alignments[bestAlignmentIndex[readEnd]];
			
			AlignInfo alignInfo = AlignSelfFullSSE(aligner, alignment, referenceSequences, preppedReads);
			
			for (int alignedLength = 10; alignedLength <= preppedReads.ReadLength(alignment.readEnd); alignedLength++)
			{
				cout << alignedLength << "\t" << alignInfo.SeqScores()[alignedLength] << endl;
			}
		}
	}
}


