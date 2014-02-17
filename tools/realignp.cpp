/*
 *  realignp.cpp
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
	double chimericPrior;
	double alignmentThreshold;
	double chimericThreshold;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Realignment Tool");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqsFilenameArg("s","seqs","Read Sequences Fastq",true,"","string",cmd);
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Sam Alignments",true,"","string",cmd);
		TCLAP::ValueArg<string> statsFilenameArg("z","stats","Stats Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> chimericPriorArg("","pchimer","Prior Probility of Chimeric Read",true,0.05,"float",cmd);
		TCLAP::ValueArg<double> alignmentThresholdArg("","talign","Alignment Posterior Threshold",true,0.1,"float",cmd);
		TCLAP::ValueArg<double> chimericThresholdArg("","tchimer","Chimeric Posterior Threshold",true,0.1,"float",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		readSeqsFilename = readSeqsFilenameArg.getValue();
		alignmentsFilename = alignmentsFilenameArg.getValue();
		statsFilename = statsFilenameArg.getValue();
		chimericPrior = chimericPriorArg.getValue();
		alignmentThreshold = alignmentThresholdArg.getValue();
		chimericThreshold = chimericThresholdArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	const int cSeedScoreThreshold = 8;
	const int cMateSearchLength = 1000;
	const int cBreakEndAdjust = 5;
	
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
		
		// Realignments
		pair<int,int> bestAlignment[2] = {pair<int,int>(0,0),pair<int,int>(0,0)};
		vector<AlignInfo> selfAlignments;
		unordered_map<int,AlignInfo> mateRevAlignments;
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			AlignInfo selfAlignInfo = AlignSelfFullSSE(aligner, alignment, referenceSequences, preppedReads);
			
			int selfSeqLength = selfAlignInfo.BestPartialSeqLength();
			int selfScore = selfAlignInfo.SeqScores()[selfSeqLength];
			
			bestAlignment[alignment.readEnd] = max(bestAlignment[alignment.readEnd], pair<int,int>(selfScore, selfSeqLength));
			
			selfAlignments.push_back(selfAlignInfo);
			
			int mateEnd = OtherReadEnd(alignment.readEnd);
			
			int seedScore;
			int seedPosition;
			AlignMate3PrimeSeed16SSE(aligner, alignment, referenceSequences, preppedReads, cMateSearchLength, seedScore, seedPosition);
			
			if (seedScore >= cSeedScoreThreshold)
			{
				AlignInfo mateFwdAlignInfo = AlignFwdMateFullSSE(aligner, alignment, referenceSequences, preppedReads, seedPosition);
				
				AlignInfo mateRevAlignInfo = AlignRevMateFullSSE(aligner, alignment, referenceSequences, preppedReads, mateFwdAlignInfo.AlignmentPosition(preppedReads.ReadLength(mateEnd)));
				
				mateRevAlignments[alignmentIndex] = mateRevAlignInfo;
			}
		}
		
		// Initialize posterior calculation for all alignments based on length of best alignment
		int alignedLength[2] = {0, 0};
		AlignmentPosterior alignPosteriors[2];
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			alignedLength[readEnd] = bestAlignment[readEnd].second - cBreakEndAdjust;
			alignPosteriors[readEnd].Initialize(&alignProbability, alignedLength[readEnd]);
		}
		
		// Identify best alignment score
		// Add alignment scores to posterior calculation
		int bestSelfScore[2] = {0, 0};
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			int selfScore = selfAlignments[alignmentIndex].SeqScores()[alignedLength[alignment.readEnd]];
			
			alignPosteriors[alignment.readEnd].AddAlignment(selfScore);
			bestSelfScore[alignment.readEnd] = max(bestSelfScore[alignment.readEnd], selfScore);
			
			unordered_map<int,AlignInfo>::const_iterator mateAlignIter = mateRevAlignments.find(alignmentIndex);
			if (mateAlignIter != mateRevAlignments.end())
			{
				int mateEnd = OtherReadEnd(alignment.readEnd);
				
				alignPosteriors[mateEnd].AddAlignment(mateAlignIter->second.SeqScores()[alignedLength[mateEnd]]);
			}
		}
		
		// Calculate probability for best discordant score
		double probBestDiscordant = alignPosteriors[0].Posterior(bestSelfScore[0]) * alignPosteriors[1].Posterior(bestSelfScore[1]);
		
		// Calculate probability for best concordant score
		double probBestConcordant = 0.0;
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			unordered_map<int,AlignInfo>::const_iterator mateAlignIter = mateRevAlignments.find(alignmentIndex);
			if (mateAlignIter != mateRevAlignments.end())
			{
				int mateEnd = OtherReadEnd(alignment.readEnd);
				
				int selfScore = selfAlignments[alignmentIndex].SeqScores()[alignedLength[alignment.readEnd]];
				int mateScore = mateAlignIter->second.SeqScores()[alignedLength[mateEnd]];
				
				probBestConcordant = max(probBestConcordant, alignPosteriors[alignment.readEnd].Posterior(selfScore) * alignPosteriors[mateEnd].Posterior(mateScore));
			}
		}
		
		// Calculate probability of concordance
		double chimericPosterior = probBestDiscordant * chimericPrior / (probBestDiscordant * chimericPrior + probBestConcordant * (1.0 - chimericPrior));
		
		if (chimericPosterior < chimericThreshold)
		{
			continue;
		}
		
		// Distribute chimeric probability between each end by multiplying square root
		double chimericPosteriorPerEnd = sqrt(chimericPosterior);
		
		// Output alignments exceeding posterior threshold
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			AlignInfo alignInfo = selfAlignments[alignmentIndex];
			
			int seqLength = alignedLength[alignment.readEnd];
			int score = alignInfo.SeqScores()[alignedLength[alignment.readEnd]];
			
			double alignmentPosterior = alignPosteriors[alignment.readEnd].Posterior(score);
			
			if (alignmentPosterior < alignmentThreshold)
			{
				continue;
			}
			
			double probability = alignmentPosterior * chimericPosteriorPerEnd;
			
			cout << alignment.fragment << "\t";
			cout << alignment.readEnd << "\t";
			cout << alignment.reference << "\t";
			cout << ((alignment.strand == PlusStrand) ? "+" : "-") << "\t";
			cout << alignInfo.AlignmentStart(seqLength) << "\t";
			cout << alignInfo.AlignmentEnd(seqLength) << "\t";
			cout << preppedReads.ReadLength(alignment.readEnd) << "\t";
			cout << seqLength << "\t";
			cout << score << "\t";
			cout << probability << endl;
		}
	}
}




