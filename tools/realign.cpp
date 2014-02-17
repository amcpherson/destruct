/*
 *  realign.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Indexer.h"
#include "ReadStream.h"
#include "AlignmentStream.h"
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


int main(int argc, char* argv[])
{
	int matchScore;
	int misMatchScore;
	int gapScore;
	int maxPaired;
	int maxUnpaired;
	string referenceFasta;
	string readSeqsFilename;
	string alignmentsFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Realignment Tool");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> maxPairedArg("p","pairedmax","Max Paired Alignments",true,0,"int",cmd);
		TCLAP::ValueArg<int> maxUnpairedArg("u","unpairedmax","Max Unpaired Alignments",true,0,"int",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqsFilenameArg("s","seq","Read Sequences",true,"","string",cmd);
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Sam Alignments",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		maxPaired = maxPairedArg.getValue();
		maxUnpaired = maxUnpairedArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		readSeqsFilename = readSeqsFilenameArg.getValue();
		alignmentsFilename = alignmentsFilenameArg.getValue();
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
	
	cerr << "Realigning" << endl;
	
	SimpleAligner aligner(matchScore, misMatchScore, gapScore);
	
	SamAlignmentStream alignmentStream(alignmentsFilename);
	FragmentAlignmentStream fragmentAlignmentStream(&alignmentStream);
	
	RawAlignmentVec alignments;
	while (fragmentAlignmentStream.GetNextAlignments(alignments))
	{
		preppedReads.SetCurrentRead(SAFEPARSE(int, alignments.front().fragment));
		
		map<int,vector<int> > alignmentScores[2];
		
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			int score = AlignSelfScoreSSE(aligner, alignment, referenceSequences, preppedReads);
			
			alignmentScores[alignment.readEnd][score].push_back(alignmentIndex);
		}
		
		IntegerVec alignedEnds;
		if (alignmentScores[0].size() == 0)
		{
			if (alignmentScores[1].rbegin()->second.size() > maxUnpaired)
			{
				continue;
			}
			
			alignedEnds.push_back(1);
		}
		else if (alignmentScores[1].size() == 0)
		{
			if (alignmentScores[0].rbegin()->second.size() > maxUnpaired)
			{
				continue;
			}
			
			alignedEnds.push_back(0);
		}
		else
		{
			if (alignmentScores[0].rbegin()->second.size() * alignmentScores[1].rbegin()->second.size() > maxPaired)
			{
				continue;
			}
			
			alignedEnds.push_back(0);
			alignedEnds.push_back(1);
		}
		
		for (IntegerVecConstIter readEndIter = alignedEnds.begin(); readEndIter != alignedEnds.end(); readEndIter++)
		{
			int readEnd = *readEndIter;
			const vector<int>& topScoring = alignmentScores[readEnd].rbegin()->second;
			for (vector<int>::const_iterator alignIter = topScoring.begin(); alignIter != topScoring.end(); alignIter++)
			{
				const RawAlignment& alignment = alignments[*alignIter];
				
				AlignInfo alignInfo = AlignSelfFullSSE(aligner, alignment, referenceSequences, preppedReads);
				
				int seqLength = alignInfo.BestPartialSeqLength();
				int score = alignInfo.SeqScores()[seqLength];
				
				DebugCheck(alignment.readEnd == readEnd);
				
				cout << alignment.fragment << "\t";
				cout << alignment.readEnd << "\t";
				cout << alignment.reference << "\t";
				cout << ((alignment.strand == PlusStrand) ? "+" : "-") << "\t";
				cout << alignInfo.AlignmentStart(seqLength) << "\t";
				cout << alignInfo.AlignmentEnd(seqLength) << "\t";
				cout << preppedReads.ReadLength(alignment.readEnd) << "\t";
				cout << seqLength << "\t";
				cout << score << endl;
			}
		}
	}
}

