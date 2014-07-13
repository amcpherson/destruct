/*
 *  realign2.cpp
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
#include "AlignmentRecord.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


int BestSplitAlignment(const short int* scores1Fwd, int scores1FwdLength, const short int* scores2Rev, int scores2RevLength, int breakInsertScore, IntegerVec& seq1Lengths, IntegerVec& seq2Lengths)
{
	IntegerVec scores1FwdBreakInsert(scores1FwdLength);
	IntegerVec scores1FwdPrevMax(scores1FwdLength);
	
	scores1FwdBreakInsert[0] = scores1Fwd[0];
	scores1FwdPrevMax[0] = 0;
	for (int idx = 1; idx < scores1FwdLength; idx++)
	{
		if (scores1FwdBreakInsert[idx-1] + breakInsertScore >= scores1Fwd[idx])
		{
			scores1FwdBreakInsert[idx] = scores1FwdBreakInsert[idx-1] + breakInsertScore;
			scores1FwdPrevMax[idx] = scores1FwdPrevMax[idx-1];
		}
		else
		{
			scores1FwdBreakInsert[idx] = scores1Fwd[idx];
			scores1FwdPrevMax[idx] = idx;
		}
	}
	
	IntegerVec seq1LengthMax;
	
	int score = numeric_limits<int>::min();
	for (int seq1Length = 0; seq1Length < scores1FwdLength; seq1Length++)
	{
		int seq2Length = scores2RevLength - seq1Length - 1;
		
		if (scores1FwdBreakInsert[seq1Length] + scores2Rev[seq2Length] > score)
		{
			score = scores1FwdBreakInsert[seq1Length] + scores2Rev[seq2Length];
			seq1LengthMax.clear();
		}
		
		if (scores1FwdBreakInsert[seq1Length] + scores2Rev[seq2Length] >= score)
		{
			seq1LengthMax.push_back(seq1Length);
		}
	}
	
	for (int idx = 0; idx < seq1LengthMax.size(); idx++)
	{
		int seq1Length = seq1LengthMax[idx];
		int seq2Length = scores2RevLength - seq1Length - 1;
		
		if (scores1FwdBreakInsert[seq1Length] == scores1Fwd[seq1Length])
		{
			seq1Lengths.push_back(seq1Length);
			seq2Lengths.push_back(seq2Length);
		}
		
		if (scores1FwdPrevMax[seq1Length] != seq1Length)
		{
			seq1Lengths.push_back(scores1FwdPrevMax[seq1Length]);
			seq2Lengths.push_back(seq2Length);
		}
	}
	
	return score;
}

int main(int argc, char* argv[])
{
	int matchScore;
	int misMatchScore;
	int gapScore;
	int minFragmentLength;
	int maxFragmentLength;
	string referenceFasta;
	string readSeqsFilename;
	string alignmentsFilename;
	string statsFilename;
	double validReadThreshold;
	double chimericPrior;
	double chimericThreshold;
	double alignmentThreshold;
	int libID;
	string spanningFilename;
	string splitFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Realignment Tool");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> minFragmentLengthArg("","flmin","Minimum Fragment Length",true,0,"int",cmd);
		TCLAP::ValueArg<int> maxFragmentLengthArg("","flmax","Maximum Fragment Length",true,0,"int",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqsFilenameArg("s","seqs","Read Sequences Fastq",true,"","string",cmd);
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Sam Alignments",true,"","string",cmd);
		TCLAP::ValueArg<string> statsFilenameArg("z","stats","Stats Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> validReadThresholdArg("","tvalid","Valid Read Threshold",true,0.01,"float",cmd);
		TCLAP::ValueArg<double> chimericPriorArg("","pchimer","Prior Probility of Chimeric Read",true,0.05,"float",cmd);
		TCLAP::ValueArg<double> chimericThresholdArg("","tchimer","Chimeric Posterior Threshold",true,0.1,"float",cmd);
		TCLAP::ValueArg<double> alignmentThresholdArg("","talign","Alignment Posterior Threshold",true,0.1,"float",cmd);
		TCLAP::ValueArg<int> libIDArg("l","lib","Library ID",true,0,"int",cmd);
		TCLAP::ValueArg<string> spanningFilenameArg("","span","Spanning Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> splitFilenameArg("","split","Splits Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		minFragmentLength = minFragmentLengthArg.getValue();
		maxFragmentLength = maxFragmentLengthArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		readSeqsFilename = readSeqsFilenameArg.getValue();
		alignmentsFilename = alignmentsFilenameArg.getValue();
		statsFilename = statsFilenameArg.getValue();
		validReadThreshold = validReadThresholdArg.getValue();
		chimericPrior = chimericPriorArg.getValue();
		chimericThreshold = chimericThresholdArg.getValue();
		alignmentThreshold = alignmentThresholdArg.getValue();
		libID = libIDArg.getValue();
		spanningFilename = spanningFilenameArg.getValue();
		splitFilename = splitFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	ofstream spanningFile(spanningFilename.c_str());
	CheckFile(spanningFile, spanningFilename);
	
	ofstream splitFile(splitFilename.c_str());
	CheckFile(splitFile, splitFilename);
	
	const int cSeedScoreThreshold = 8;
	const int cMaxBreakpointHomology = 25;
	const int cMinAnchor = 8;
	const int cBreakEndAdjust = 5;
	
	cerr << "Reading alignment stats" << endl;
	
	AlignmentProbability alignProbability(matchScore);
	alignProbability.ReadDistributions(statsFilename, validReadThreshold);
	
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
		pair<int,int> bestAlignment[2] = {pair<int,int>(0,0),pair<int,int>(0,0)};
		vector<AlignInfo> selfAlignments;
		unordered_map<int,AlignInfo> mateFwdAlignments;
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
			AlignMate3PrimeSeed16SSE(aligner, alignment, referenceSequences, preppedReads, minFragmentLength, maxFragmentLength, seedScore, seedPosition);
			
			if (seedScore >= cSeedScoreThreshold)
			{
				AlignInfo mateFwdAlignInfo = AlignFwdMateFullSSE(aligner, alignment, referenceSequences, preppedReads, seedPosition);
				
				mateFwdAlignments[alignmentIndex] = mateFwdAlignInfo;
				
				AlignInfo mateRevAlignInfo = AlignRevMateFullSSE(aligner, alignment, referenceSequences, preppedReads, mateFwdAlignInfo.AlignmentPosition(preppedReads.ReadLength(mateEnd)));
				
				mateRevAlignments[alignmentIndex] = mateRevAlignInfo;
			}
		}

		//
		// Calculate adjusted aligned length for each end
		// Also create an array of read lengths
		//
		vector<int> alignedLength(2);
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			alignedLength[readEnd] = bestAlignment[readEnd].second - cBreakEndAdjust;
		}
		
		//
		// Identify best full alignment scores
		//
		// Add partial alignment scores to posterior calculation including:
		//  - self alignments
		//  - reverse mate alignments, if they exist
		// Calculate best full alignment scores
		//
		int bestSelfScoreFull[2] = {0, 0};
		AlignmentPosterior alignPosteriorPartial(alignProbability, chimericPrior, alignedLength);
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			int selfScorePart = selfAlignments[alignmentIndex].SeqScores()[alignedLength[alignment.readEnd]];
			int selfScoreFull = selfAlignments[alignmentIndex].SeqScores()[preppedReads.ReadLength(alignment.readEnd)];

			bestSelfScoreFull[alignment.readEnd] = max(bestSelfScoreFull[alignment.readEnd], selfScoreFull);
			
			unordered_map<int,AlignInfo>::const_iterator mateAlignIter = mateRevAlignments.find(alignmentIndex);
			if (mateAlignIter == mateRevAlignments.end())
			{
				alignPosteriorPartial.AppendAlignment(alignment.readEnd, selfScorePart);
			}
			else
			{
				int mateEnd = OtherReadEnd(alignment.readEnd);

				int mateScore = mateAlignIter->second.SeqScores()[alignedLength[mateEnd]];
				
				alignPosteriorPartial.AppendAlignmentWithMate(alignment.readEnd, selfScorePart, mateScore);
			}
		}

		// 
		// Filter concordant alignments based on threshold on 
		// the posterior probability that any of the concordant
		// alignments are above threshold
		// 
		if (alignPosteriorPartial.PosteriorConcordant() >= chimericThreshold)
		{
			continue;
		}
		
		//
		// Create list of alignment indices for each end, filter
		// based on partial alignment posterior
		//
		vector<int> alignmentIndices[2];
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];

			double alignmentPosterior = alignPosteriorPartial.Posterior(alignmentIndex);

			if (alignmentPosterior < alignmentThreshold)
			{
				continue;
			}
			
			alignmentIndices[alignment.readEnd].push_back(alignmentIndex);
		}

		//
		// Filter reads unless there is at least 1 alignment of each end exceeding
		// partial alignment posterior threshold
		//
		if (alignmentIndices[0].empty() || alignmentIndices[1].empty())
		{
			continue;
		}
		
		//
		// Calculate split alignments for selected alignment pairs
		//
		// Update best full self score
		// 
		// Store split scores and length score arrays
		//
		unordered_map<int,int> splitScores[2];
		unordered_map<int,IntegerVec> splitSeq1Length[2];
		unordered_map<int,IntegerVec> splitSeq2Length[2];
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			int mateEnd = OtherReadEnd(readEnd);
			
			for (vector<int>::const_iterator selfAlignmentIter = alignmentIndices[readEnd].begin(); selfAlignmentIter != alignmentIndices[readEnd].end(); selfAlignmentIter++)
			{
				int selfAlignmentIndex = *selfAlignmentIter;
				const AlignInfo& selfAlignInfo = selfAlignments[selfAlignmentIndex];
				const RawAlignment& selfAlignment = alignments[selfAlignmentIndex];
				
				for (vector<int>::const_iterator mateAlignmentIter = alignmentIndices[mateEnd].begin(); mateAlignmentIter != alignmentIndices[mateEnd].end(); mateAlignmentIter++)
				{
					int mateAlignmentIndex = *mateAlignmentIter;
					unordered_map<int,AlignInfo>::const_iterator mateAlignInfoIter = mateFwdAlignments.find(mateAlignmentIndex);
					
					if (mateAlignInfoIter == mateFwdAlignments.end())
					{
						continue;
					}
					
					const AlignInfo& mateAlignInfo = mateAlignInfoIter->second;
					const RawAlignment& mateAlignment = alignments[mateAlignmentIndex];
					
					IntegerVec seq1Length;
					IntegerVec seq2Length;
					int score = BestSplitAlignment(selfAlignInfo.SeqScores(), selfAlignInfo.SeqScoresLength(), mateAlignInfo.SeqScores(), mateAlignInfo.SeqScoresLength(), -1, seq1Length, seq2Length);

					splitScores[readEnd][selfAlignmentIndex] = score;
					swap(seq1Length, splitSeq1Length[readEnd][selfAlignmentIndex]);
					swap(seq2Length, splitSeq2Length[readEnd][selfAlignmentIndex]);

					bestSelfScoreFull[selfAlignment.readEnd] = max(bestSelfScoreFull[selfAlignment.readEnd], score);
				}
			}
		}

		//
		// Filter reads based on the cdf of the score likelihood
		//
		if (!alignProbability.AboveThreshold(preppedReads.ReadLength(0), bestSelfScoreFull[0]) ||
			!alignProbability.AboveThreshold(preppedReads.ReadLength(1), bestSelfScoreFull[1]))
		{
			continue;
		}
		
		// Calculate split alignments exceeding posterior threshold
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			int mateEnd = OtherReadEnd(readEnd);
			
			for (vector<int>::const_iterator selfAlignmentIter = alignmentIndices[readEnd].begin(); selfAlignmentIter != alignmentIndices[readEnd].end(); selfAlignmentIter++)
			{
				int selfAlignmentIndex = *selfAlignmentIter;
				const AlignInfo& selfAlignInfo = selfAlignments[selfAlignmentIndex];
				const RawAlignment& selfAlignment = alignments[selfAlignmentIndex];
				
				DebugCheck(selfAlignment.readEnd == readEnd);

				for (vector<int>::const_iterator mateAlignmentIter = alignmentIndices[mateEnd].begin(); mateAlignmentIter != alignmentIndices[mateEnd].end(); mateAlignmentIter++)
				{
					int mateAlignmentIndex = *mateAlignmentIter;
					unordered_map<int,AlignInfo>::const_iterator mateAlignInfoIter = mateFwdAlignments.find(mateAlignmentIndex);
					
					if (mateAlignInfoIter == mateFwdAlignments.end())
					{
						continue;
					}
					
					const AlignInfo& mateAlignInfo = mateAlignInfoIter->second;
					const RawAlignment& mateAlignment = alignments[mateAlignmentIndex];
					
					DebugCheck(mateAlignment.readEnd == mateEnd);

					int score = splitScores[readEnd][selfAlignmentIndex];
					const IntegerVec& seq1Length = splitSeq1Length[readEnd][selfAlignmentIndex];
					const IntegerVec& seq2Length = splitSeq2Length[readEnd][selfAlignmentIndex];

					string readSeq = preppedReads.Sequence(readEnd);
					
					for (int i = 0; i < seq1Length.size(); i++)
					{
						if (seq1Length[i] < cMinAnchor || seq2Length[i] < cMinAnchor)
						{
							continue;
						}

						SplitAlignmentRecord record;
						record.libID = libID;
						record.readID = readID;
						record.readEnd = selfAlignment.readEnd;
						record.alignID[readEnd] = selfAlignmentIndex;
						record.chromosome[readEnd] = selfAlignment.reference;
						record.strand[readEnd] = ((selfAlignment.strand == PlusStrand) ? "+" : "-");
						record.position[readEnd] = selfAlignInfo.BreakPosition(seq1Length[i]);
						record.alignID[mateEnd] = mateAlignmentIndex;
						record.chromosome[mateEnd] = mateAlignment.reference;
						record.strand[mateEnd] = ((mateAlignment.strand == PlusStrand) ? "+" : "-");
						record.position[mateEnd] = mateAlignInfo.BreakPosition(seq2Length[i]);
						record.inserted = readSeq.substr(seq1Length[i], readSeq.size() - seq2Length[i] - seq1Length[i]);
						record.score = score;

						splitFile << record;
					}
				}
			}
		}

		// Output spanning alignments exceeding posterior threshold
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			for (vector<int>::const_iterator alignmentIter = alignmentIndices[readEnd].begin(); alignmentIter != alignmentIndices[readEnd].end(); alignmentIter++)
			{
				int alignmentIndex = *alignmentIter;

				const RawAlignment& alignment = alignments[alignmentIndex];
				
				AlignInfo alignInfo = selfAlignments[alignmentIndex];
				
				int seqLength = alignedLength[alignment.readEnd];
				int score = alignInfo.SeqScores()[alignedLength[alignment.readEnd]];
				
				SpanningAlignmentRecord record;
				record.libID = libID;
				record.readID = readID;
				record.readEnd = alignment.readEnd;
				record.alignID = alignmentIndex;
				record.chromosome = alignment.reference;
				record.strand = ((alignment.strand == PlusStrand) ? "+" : "-");
				record.start = alignInfo.AlignmentStart(seqLength);
				record.end = alignInfo.AlignmentEnd(seqLength);
				record.score = score;

				spanningFile << record;
			}
		}
	}
}


