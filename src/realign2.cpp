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


bool BestSplitAlignment(const short int* scores1Fwd, int scores1FwdLength,
                        const short int* scores2Rev, int scores2RevLength,
                        int breakInsertScore, int minAnchor, int& score,
                        int& seq1Length, int& seq2Length)
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
	
	score = numeric_limits<int>::min();
	for (seq1Length = 0; seq1Length < scores1FwdLength; seq1Length++)
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
		seq1Length = seq1LengthMax[idx];
		seq2Length = scores2RevLength - seq1Length - 1;
		
		if (scores1FwdBreakInsert[seq1Length] == scores1Fwd[seq1Length])
		{
			if ((seq1Length >= minAnchor) && (seq2Length >= minAnchor))
			{
				return true;
			}
		}
		
		if (scores1FwdPrevMax[seq1Length] != seq1Length)
		{
			seq1Length = scores1FwdPrevMax[seq1Length];

			if ((seq1Length >= minAnchor) && (seq2Length >= minAnchor))
			{
				return true;
			}
		}
	}
	
	return false;
}


void Complement(char& nucleotide)
{
	switch (nucleotide)
	{
		case 'A': nucleotide = 'T'; break;
		case 'C': nucleotide = 'G'; break;
		case 'T': nucleotide = 'A'; break;
		case 'G': nucleotide = 'C'; break;
		case 'a': nucleotide = 't'; break;
		case 'c': nucleotide = 'g'; break;
		case 't': nucleotide = 'a'; break;
		case 'g': nucleotide = 'c'; break;
	}
}


int CalculateOffset(const string& strand, int offset)
{
	int dir = (strand == "+") ? 1 : -1;
	return offset * dir;
}


int CalculateForwardHomology(const Sequences& sequences,
                             const string (&chromosome)[2],
                             const string (&strand)[2],
                             const int (&position)[2],
                             int maxOffset,
                             bool flip=false)
{
	int idx1 = (flip) ? 1 : 0;
	int idx2 = 1 - idx1;

	const char* seqPtr1 = sequences.Get(chromosome[idx1], position[idx1]);
	const char* seqPtr2 = sequences.Get(chromosome[idx2], position[idx2]);

	int homology = 0;
	for (int offset = 1; offset <= maxOffset; offset++)
	{
		char nt1 = *(seqPtr1 + CalculateOffset(strand[idx1], offset));
		char nt2 = *(seqPtr2 + CalculateOffset(strand[idx2], 1 - offset));
		
		if (strand[idx1] != "+")
		{
			Complement(nt1);
		}
		
		if (strand[idx2] != "-")
		{
			Complement(nt2);
		}
		
		if (nt1 != nt2)
		{
			break;
		}
		
		homology = offset;
	}
	
	return homology;
}


void HomologyConsistentBreakpoint(const Sequences& sequences,
                                  const string (&chromosome)[2],
                                  const string (&strand)[2],
                                  int (&position)[2],
                                  int& homology,
                                  int maxOffset)
{
	int maxOffsetA = CalculateForwardHomology(sequences, chromosome, strand, position, maxOffset, false);
	int maxOffsetB = CalculateForwardHomology(sequences, chromosome, strand, position, maxOffset, true);

	// Ensure that the same breakpoint is selected among the multiple
	// breakpoints possible when there is breakpoint homology.  Always
	// select the breakpoint for which the minimum of the two breakend
	// positions is minimal

	int positionA1 = position[0] + CalculateOffset(strand[0], maxOffsetA);
	int positionA2 = position[1] + CalculateOffset(strand[1], -maxOffsetA);

	int positionB1 = position[0] + CalculateOffset(strand[0], -maxOffsetB);
	int positionB2 = position[1] + CalculateOffset(strand[1], maxOffsetB);

	if (min(positionA1, positionA2) < min(positionB1, positionB2))
	{
		position[0] = positionA1;
		position[1] = positionA2;
	}
	else
	{
		position[0] = positionB1;
		position[1] = positionB2;
	}

	homology = maxOffsetA + maxOffsetB;
}


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
		TCLAP::ValueArg<string> reads1FilenameArg("1","reads1","Read End 1 Fastq",true,"","string",cmd);
		TCLAP::ValueArg<string> reads2FilenameArg("2","reads2","Read End 2 Fastq",true,"","string",cmd);
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
		reads1Filename = reads1FilenameArg.getValue();
		reads2Filename = reads2FilenameArg.getValue();
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
	const int cMinAnchor = 8;
	const int cBreakEndAdjust = 5;
	const int cInsertedPenalty = -1;
	
	cerr << "Reading alignment stats" << endl;
	
	AlignmentProbability alignProbability(matchScore);
	alignProbability.ReadDistributions(statsFilename, validReadThreshold);

	if (!alignProbability.Good())
	{
	        cerr << "Warning: empty stats" << endl;

		// Make sure we parse alignments from std
		if (alignmentsFilename == "-")
		{
			for (string line; getline(cin, line);) {}
		}

		exit(0);
	}

	int minAlignedLength = alignProbability.GetMinAlignedLength();
	
	cerr << "Reading reference fasta" << endl;
	
	Sequences referenceSequences(4000);
	referenceSequences.Read(referenceFasta);
	
	cerr << "Reading fastq sequences" << endl;
	
	FastqReadStream reads1Stream(reads1Filename);
	FastqReadStream reads2Stream(reads2Filename);
	
	PreppedReads preppedReads;
	preppedReads.Prep(reads1Stream);
	preppedReads.Prep(reads2Stream);
	
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

		// Check for very poor alignments
		if (min(alignedLength[0], alignedLength[1]) < minAlignedLength)
		{
			continue;
		}
		
		//
		// Add partial alignment scores to posterior calculation including:
		//  - self alignments
		//  - reverse mate alignments, if they exist
		// Calculate best full alignment scores
		//
		AlignmentPosterior alignPosteriorPartial(alignProbability, chimericPrior, alignedLength);
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			int selfScorePart = selfAlignments[alignmentIndex].SeqScores()[alignedLength[alignment.readEnd]];
			
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
		// Calculate the full alignment score of each alignment and threshold
		// on the CDF of the likelihood of attaining this alignment.  Mark each
		// alignment index as passing or not passing the CDF threshold.  Also 
		// mark read ends as having an alignment that passes the CDF threshold.
		//
		// Create list of alignment indices for each end, filter
		// based on partial alignment posterior
		//
		bool validSpanningReadEnd[2] = {false, false};
		vector<bool> validSpanningAlignment(alignments.size(), false);
		vector<int> alignmentIndices[2];
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];

			int selfScoreFull = selfAlignments[alignmentIndex].SeqScores()[preppedReads.ReadLength(alignment.readEnd)];

			if (alignProbability.AboveThreshold(preppedReads.ReadLength(alignment.readEnd), selfScoreFull))
			{
				validSpanningReadEnd[alignment.readEnd] = true;
				validSpanningAlignment[alignmentIndex] = true;
			}

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
		// Calculate the best split length and score where best split length equals
		// read length unless there are insertions at the breakpoint, in which case
		// the insertion length is subtracked from the read length and the insertion
		// penalty is subtracted from the score.  Check whether the length and score
		// pass the CDF threshold.  Output the split read only it passes the 
		// threshold.  Also store a boolean as True for split reads that pass the
		// threshold to be used to determine if a spanning read should be output 
		// regardless of the fully aligned score.  Store a boolean for each read end
		// if a valid split alignment was found for that read end.
		//
		bool validSplitReadEnd[2] = {false, false};
		vector<bool> validSplitAlignment(alignments.size(), false);
		vector<SplitAlignmentRecord> splitRecords;
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
					
					DebugCheck(mateAlignment.readEnd == mateEnd);

					int score;
					int seq1Length;
					int seq2Length;
					bool hasSplit = BestSplitAlignment(selfAlignInfo.SeqScores(), selfAlignInfo.SeqScoresLength(),
					                                   mateAlignInfo.SeqScores(), mateAlignInfo.SeqScoresLength(),
					                                   cInsertedPenalty, cMinAnchor, score, seq1Length, seq2Length);

					if (!hasSplit)
					{
						continue;
					}

					string readSeq = preppedReads.Sequence(readEnd);

					string inserted = readSeq.substr(seq1Length, readSeq.size() - seq2Length - seq1Length);

					int alignedLength = readSeq.size() - inserted.size();
					int alignedScore = score - cInsertedPenalty * (int)inserted.size();

					if (!alignProbability.AboveThreshold(alignedLength, alignedScore))
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
					record.position[readEnd] = selfAlignInfo.BreakPosition(seq1Length);
					record.alignID[mateEnd] = mateAlignmentIndex;
					record.chromosome[mateEnd] = mateAlignment.reference;
					record.strand[mateEnd] = ((mateAlignment.strand == PlusStrand) ? "+" : "-");
					record.position[mateEnd] = mateAlignInfo.BreakPosition(seq2Length);
					record.homology = 0;
					record.inserted = readSeq.substr(seq1Length, readSeq.size() - seq2Length - seq1Length);
					record.score = score;

					if (record.inserted.empty())
					{
						HomologyConsistentBreakpoint(referenceSequences, record.chromosome, record.strand,
						                             record.position, record.homology, readSeq.size());
					}

					splitRecords.push_back(record);

					validSplitReadEnd[readEnd] = true;
					validSplitAlignment[selfAlignmentIndex] = true;
				}
			}
		}

		//
		// Filter reads based on the cdf of the score likelihood
		//
		if ((!validSpanningReadEnd[0] && !validSplitReadEnd[0]) || (!validSpanningReadEnd[1] && !validSplitReadEnd[1]))
		{
			continue;
		}

		//
		// Output split alignment records after cdf test
		//
		for (vector<SplitAlignmentRecord>::const_iterator recordIter = splitRecords.begin(); recordIter != splitRecords.end(); recordIter++)
		{
			splitFile << *recordIter;
		}

		// Output spanning alignments exceeding posterior threshold
		// Reads must also have a valid spanning or split read
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			for (vector<int>::const_iterator alignmentIter = alignmentIndices[readEnd].begin(); alignmentIter != alignmentIndices[readEnd].end(); alignmentIter++)
			{
				int alignmentIndex = *alignmentIter;

				const RawAlignment& alignment = alignments[alignmentIndex];
				
				DebugCheck(alignment.readEnd == readEnd);
				
				int mateEnd = OtherReadEnd(readEnd);

				AlignInfo alignInfo = selfAlignments[alignmentIndex];
				
				int selfSeqLength = alignInfo.BestPartialSeqLength();

				int mateScore = 0;

				unordered_map<int,AlignInfo>::const_iterator mateAlignIter = mateRevAlignments.find(alignmentIndex);
				if (mateAlignIter != mateRevAlignments.end())
				{
					mateScore = max((short)0, mateAlignIter->second.SeqScores()[alignedLength[mateEnd]]);
				}

				if (!validSpanningAlignment[alignmentIndex] && !validSplitAlignment[alignmentIndex])
				{
					continue;
				} 
				
				SpanningAlignmentRecord record;
				record.libID = libID;
				record.readID = readID;
				record.readEnd = alignment.readEnd;
				record.alignID = alignmentIndex;
				record.chromosome = alignment.reference;
				record.strand = ((alignment.strand == PlusStrand) ? "+" : "-");
				record.position = alignInfo.OuterPosition();
				record.alignedLength = selfSeqLength;
				record.mateLength = preppedReads.ReadLength(mateEnd);
				record.mateScore = mateScore;

				spanningFile << record;
			}
		}
	}
}


