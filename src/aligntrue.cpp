/*
 *  aligntrue.cpp
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Indexer.h"
#include "AlignmentStream.h"
#include "AlignRead.h"
#include "Sequences.h"
#include "SimpleAligner.h"

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
	int minFragmentLength;
	int maxFragmentLength;
	string referenceFasta;
	string reads1Filename;
	string reads2Filename;
	string alignmentsFilename;
	string scoresFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Read Alignment True Distribution");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> minFragmentLengthArg("","flmin","Minimum Fragment Length",true,0,"int",cmd);
		TCLAP::ValueArg<int> maxFragmentLengthArg("","flmax","Maximum Fragment Length",true,0,"int",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> reads1FilenameArg("1","reads1","Read End 1 Fastq",true,"","string",cmd);
		TCLAP::ValueArg<string> reads2FilenameArg("2","reads2","Read End 2 Fastq",true,"","string",cmd);
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Sam Alignments",true,"","string",cmd);
		TCLAP::ValueArg<string> scoresFilenameArg("s","scores","Output Scores Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		minFragmentLength = minFragmentLengthArg.getValue();
		maxFragmentLength = maxFragmentLengthArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		reads1Filename = reads1FilenameArg.getValue();
		reads2Filename = reads2FilenameArg.getValue();
		alignmentsFilename = alignmentsFilenameArg.getValue();
		scoresFilename = scoresFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	const int cSeedScoreThreshold = 8;

	cerr << "Reading reference fasta" << endl;
	
	Sequences referenceSequences(1000);
	referenceSequences.Read(referenceFasta);
	
	cerr << "Reading fastq sequences" << endl;
	
	FastqReadStream reads1Stream(reads1Filename);
	FastqReadStream reads2Stream(reads2Filename);
	
	PreppedReads preppedReads;
	preppedReads.Prep(reads1Stream);
	preppedReads.Prep(reads2Stream);
	
	cerr << "Realigning" << endl;
	
	ofstream scoresFile(scoresFilename.c_str());
	CheckFile(scoresFile, scoresFilename);
	
	SimpleAligner aligner(matchScore, misMatchScore, gapScore);
	
	SamAlignmentStream alignmentStream(alignmentsFilename);
	FragmentAlignmentStream fragmentAlignmentStream(&alignmentStream);
	
	RawAlignmentVec alignments;
	while (fragmentAlignmentStream.GetNextAlignments(alignments))
	{
		int readID = SAFEPARSE(int, alignments.front().fragment);

		preppedReads.SetCurrentRead(readID);
		
		// Check that both ends are mapped
		bool readEndMapped[2] = {false,false};
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			readEndMapped[alignment.readEnd] = true;
		}
		if (!readEndMapped[0] || !readEndMapped[1])
		{
			continue;
		}
		
		//
		// Realignments
		//
		// Calculate the 'self' alignment, the alignment of the read
		// to its seed match location
		//
		// Calculate the 'mate' seed alignment, the alignment of a 16 nt seed
		// from the mate read to the region near the self seed location
		//
		// Given a reasonable mate seed alignment, calculate the 'forward'
		// mate alignment, the alignment of the mate read forward from the
		// seed match location
		//
		// Also calculate the 'reverse' mate alignment, the alignment of the
		// mate starting from the ending point of the forward mate alignment
		//
		// Keep the alignment info for the best self alignment and reverse
		// mate alignment where best is that with the highest total score.
		//
		int bestScore = -1;
		AlignInfo bestAlignments[2];
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			AlignInfo selfAlignInfo = AlignSelfFullSSE(aligner, alignment, referenceSequences, preppedReads);
			
			int selfScore = selfAlignInfo.SeqScores()[preppedReads.ReadLength(alignment.readEnd)];

			int mateEnd = OtherReadEnd(alignment.readEnd);
			
			int seedScore;
			int seedPosition;
			AlignMate3PrimeSeed16SSE(aligner, alignment, referenceSequences, preppedReads, minFragmentLength, maxFragmentLength, seedScore, seedPosition);
			
			if (seedScore >= cSeedScoreThreshold)
			{
				AlignInfo mateFwdAlignInfo = AlignFwdMateFullSSE(aligner, alignment, referenceSequences, preppedReads, seedPosition);
				
				AlignInfo mateRevAlignInfo = AlignRevMateFullSSE(aligner, alignment, referenceSequences, preppedReads, mateFwdAlignInfo.AlignmentPosition(preppedReads.ReadLength(mateEnd)));

				int mateScore = mateRevAlignInfo.SeqScores()[preppedReads.ReadLength(mateEnd)];

				if (selfScore + mateScore > bestScore)
				{
					bestScore = selfScore + mateScore;
					bestAlignments[alignment.readEnd] = selfAlignInfo;
					bestAlignments[mateEnd] = mateRevAlignInfo;
				}
			}
		}

		// Skip reads that dont have any reasonable concordant alignment
		if (bestScore < 0)
		{
			continue;
		}

		// Output alignment scores for all read lengths for both ends
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			int readLength = preppedReads.ReadLength(readEnd);

			for (int alignedLength = 10; alignedLength <= readLength; alignedLength++)
			{
				int score = bestAlignments[readEnd].SeqScores()[alignedLength];

				scoresFile << alignedLength << "\t" << score << endl;
			}
		}
	}
}



