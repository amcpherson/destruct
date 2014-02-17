/*
 *  calcpriorgenomic.cpp
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Indexer.h"
#include "ReadStream.h"
#include "AlignmentStream.h"
#include "Sequences.h"
#include "SimpleAligner.h"
#include "AlignmentProbability.h"
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
	string referenceFasta;
	string readSeqsFilename;
	string alignmentsFilename;
	string statsFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Realignment Tool");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqsFilenameArg("s","seq","Read Sequences",true,"","string",cmd);
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Sam Alignments",true,"","string",cmd);
		TCLAP::ValueArg<string> statsFilenameArg("z","stats","Stats Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		readSeqsFilename = readSeqsFilenameArg.getValue();
		alignmentsFilename = alignmentsFilenameArg.getValue();
		statsFilename = statsFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	cerr << "Reading alignment stats" << endl;
	
	AlignmentProbability alignProbability(matchScore);
	alignProbability.ReadDistributions(statsFilename);
	
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
	
	IntegerVec seqLengths;
	IntegerVec scores;
	
	RawAlignmentVec alignments;
	while (fragmentAlignmentStream.GetNextAlignments(alignments))
	{
		preppedReads.SetCurrentRead(SAFEPARSE(int, alignments.front().fragment));
		
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
		
		// Quick scoring of all alignments
		multimap<int,int> alignmentScores[2];
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			int score = AlignSelfScoreSSE(aligner, alignment, referenceSequences, preppedReads);
			
			alignmentScores[alignment.readEnd].insert(make_pair(score, alignmentIndex));
		}
		
		// Select best scoring end
		int bestReadEnd = (alignmentScores[0].rbegin()->first > alignmentScores[1].rbegin()->first) ? 0 : 1;
		
		const RawAlignment& alignment = alignments[alignmentScores[bestReadEnd].rbegin()->second];
		
		AlignInfo alignInfo = AlignSelfFullSSE(aligner, alignment, referenceSequences, preppedReads);
		
		int seqLength = preppedReads.ReadLength(alignment.readEnd);
		int score = alignInfo.SeqScores()[seqLength];
		
		seqLengths.push_back(seqLength);
		scores.push_back(score);
	}
	
	double priorGenomic = 0.5;
	for (int i = 0; i < 20; i++)
	{
		double probGenomicSum = 0.0;
		for (int alignIndex = 0; alignIndex < seqLengths.size(); alignIndex++)
		{
			probGenomicSum += alignProbability.Classify(seqLengths[alignIndex], scores[alignIndex], priorGenomic);
		}
		
		priorGenomic = probGenomicSum / (double)seqLengths.size();
	}
	
	cout << priorGenomic << endl;
}




