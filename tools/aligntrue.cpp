/*
 *  aligntrue.cpp
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Indexer.h"
#include "AlignmentStream.h"
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
	string alignmentsFilename;
	string referenceFasta;
	
	try
	{
		TCLAP::CmdLine cmd("Read Alignment True Distribution");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Sam Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		alignmentsFilename = alignmentsFilenameArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
		
	cerr << "Read reference fasta" << endl;
	
	Sequences referenceSequences(100);
	referenceSequences.Read(referenceFasta);
	
	SimpleAligner aligner(matchScore, misMatchScore, gapScore);
	
	AlignmentStream* alignmentStream = new SamAlignmentStream(alignmentsFilename);
	FragmentAlignmentStream fragmentAlignments(alignmentStream);
	
	cerr << "Realigning to calculate partial scores" << endl;
	
	RawAlignmentVec alignments;
	while (fragmentAlignments.GetNextAlignments(alignments))
	{
		if (alignments.size() != 2)
		{
			cerr << "Error: Must be paired concordant alignments" << endl;
		}
		
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			int start = alignments[readEnd].region.start - 20;
			int end = alignments[readEnd].region.end + 20;
			
			string referenceSequence;
			referenceSequences.Get(alignments[readEnd].reference, start, end, referenceSequence);
			
			IntegerVec scores;
			aligner.AlignPartial(referenceSequence, alignments[readEnd].sequence, scores);
			
			for (int alignedLength = 10; alignedLength <= (int)alignments[readEnd].sequence.size(); alignedLength++)
			{
				cout << alignedLength << "\t" << scores[alignedLength] << endl;
			}
		}
	}
}

