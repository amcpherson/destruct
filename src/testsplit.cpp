/*
 *  testsplit.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "AlignmentStream.h"
#include "ReadStream.h"
#include "Sequences.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


IntegerVec mCurrent;
IntegerVec mPrevious;
int mMatchScore = 2;
int mMisMatchScore = -3;
int mGapScore = -4;

void AlignPartial(const string& reference, const string& sequence, IntegerVec& seqScores, IntegerVec& seqCounts, IntegerVec& refPositions)
{
	if (sequence.size() == 0)
	{
		return;
	}
	
	int matrixLength = reference.size() + 1;
	int matrixHeight = sequence.size() + 1;
	
	if (matrixLength != mCurrent.size())
	{
		mCurrent.resize(matrixLength);
		mPrevious.resize(matrixLength);
	}
	
	seqScores.resize(sequence.size() + 1);
	seqCounts.resize(sequence.size() + 1);
	refPositions.resize(sequence.size() + 1);
	
	for (int i = 0; i < matrixLength; i++) 
	{
		mPrevious[i] = 0;
	}
	
	seqScores[0] = 0;
	seqCounts[0] = matrixLength;
	refPositions[0] = 0;
	
	for (int j = 1; j < matrixHeight - 1; j++)
	{
		int seqPos = j - 1;
		
		mCurrent[0] = mPrevious[0] + mGapScore;
		
		seqScores[j] = mPrevious[0] + mGapScore;
		seqCounts[j] = 1;
		
		for (int i = 1; i < matrixLength; i++) 
		{
			int refPos = i - 1;
			
			int matchScore = mPrevious[i-1] + ((reference[refPos] == sequence[seqPos]) ? mMatchScore : mMisMatchScore);
			int gapRefScore = mCurrent[i-1] + mGapScore;
			int gapReadScore = mPrevious[i] + mGapScore;
			int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
			
			if (bestScore > seqScores[j])
			{
				seqScores[j] = bestScore;
				seqCounts[j] = 1;
				refPositions[j] = refPos;
			}
			else if (bestScore == seqScores[j])
			{
				seqCounts[j]++;
			}
			
			mCurrent[i] = bestScore;
		}
		
		swap(mCurrent, mPrevious);
	}
	
	int seqPos = matrixHeight - 2;
	
	mCurrent[0] = mPrevious[0] + mGapScore;
	
	seqScores[matrixHeight - 1] = mPrevious[0] + mGapScore;
	seqCounts[matrixHeight - 1] = 1;
	
	for (int i = 1; i < matrixLength; i++) 
	{
		int refPos = i - 1;
		
		int matchScore = mPrevious[i-1] + ((reference[refPos] == sequence[seqPos]) ? mMatchScore : mMisMatchScore);
		int gapRefScore = mCurrent[i-1] + mGapScore;
		int gapReadScore = mPrevious[i] + mGapScore;
		int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
		
		if (bestScore > seqScores[matrixHeight - 1])
		{
			seqScores[matrixHeight - 1] = bestScore;
			seqCounts[matrixHeight - 1] = 1;
			refPositions[matrixHeight - 1] = refPos;
		}
		else if (bestScore == seqScores[matrixHeight - 1])
		{
			seqCounts[matrixHeight - 1]++;
		}
		
		mCurrent[i] = bestScore;
	}
}

void AlignPartialFwd(const char* refStart, const char* refEnd, const char* seqStart, const char* seqEnd, IntegerVec& seqScores, IntegerVec& seqCounts, IntegerVec& refPositions)
{
	if (seqStart >= seqEnd)
	{
		return;
	}
	
	int matrixLength = refEnd - refStart + 1;
	int matrixHeight = seqEnd - seqStart + 1;
	
	if (matrixLength != mCurrent.size())
	{
		mCurrent.resize(matrixLength);
		mPrevious.resize(matrixLength);
	}
	
	seqScores.resize(matrixHeight);
	seqCounts.resize(matrixHeight);
	refPositions.resize(matrixHeight);
	
	for (int i = 0; i < matrixLength; i++) 
	{
		mPrevious[i] = 0;
	}
	
	seqScores[0] = 0;
	seqCounts[0] = matrixLength;
	refPositions[0] = 0;
	
	for (int j = 1; j < matrixHeight - 1; j++)
	{
		int seqPos = j - 1;
		
		mCurrent[0] = mPrevious[0] + mGapScore;
		
		seqScores[j] = mPrevious[0] + mGapScore;
		seqCounts[j] = 1;
		
		for (int i = 1; i < matrixLength; i++) 
		{
			int refPos = i - 1;
			
			int matchScore = mPrevious[i-1] + ((refStart[refPos] == seqEnd[-seqPos - 1]) ? mMatchScore : mMisMatchScore);
			int gapRefScore = mCurrent[i-1] + mGapScore;
			int gapReadScore = mPrevious[i] + mGapScore;
			int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
			
			if (bestScore > seqScores[j])
			{
				seqScores[j] = bestScore;
				seqCounts[j] = 1;
				refPositions[j] = refPos;
			}
			else if (bestScore == seqScores[j])
			{
				seqCounts[j]++;
			}
			
			mCurrent[i] = bestScore;
		}
		
		swap(mCurrent, mPrevious);
	}
	
	int seqPos = matrixHeight - 2;
	
	mCurrent[0] = mPrevious[0] + mGapScore;
	
	seqScores[matrixHeight - 1] = mPrevious[0] + mGapScore;
	seqCounts[matrixHeight - 1] = 1;
	
	for (int i = 1; i < matrixLength; i++) 
	{
		int refPos = i - 1;
		
		int matchScore = mPrevious[i-1] + ((refStart[refPos] == seqEnd[-seqPos - 1]) ? mMatchScore : mMisMatchScore);
		int gapRefScore = mCurrent[i-1] + mGapScore;
		int gapReadScore = mPrevious[i] + mGapScore;
		int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
		
		if (bestScore > seqScores[matrixHeight - 1])
		{
			seqScores[matrixHeight - 1] = bestScore;
			seqCounts[matrixHeight - 1] = 1;
			refPositions[matrixHeight - 1] = refPos;
		}
		else if (bestScore == seqScores[matrixHeight - 1])
		{
			seqCounts[matrixHeight - 1]++;
		}
		
		mCurrent[i] = bestScore;
	}
}

void AlignPartialRev(const char* refStart, const char* refEnd, const char* seqStart, const char* seqEnd, IntegerVec& seqScores, IntegerVec& seqCounts, IntegerVec& refPositions)
{
	if (seqStart >= seqEnd)
	{
		return;
	}
	
	int matrixLength = refEnd - refStart + 1;
	int matrixHeight = seqEnd - seqStart + 1;
	
	if (matrixLength != mCurrent.size())
	{
		mCurrent.resize(matrixLength);
		mPrevious.resize(matrixLength);
	}
	
	seqScores.resize(matrixHeight);
	seqCounts.resize(matrixHeight);
	refPositions.resize(matrixHeight);
	
	for (int i = 0; i < matrixLength; i++) 
	{
		mPrevious[i] = 0;
	}
	
	seqScores[0] = 0;
	seqCounts[0] = matrixLength;
	refPositions[0] = 0;
	
	for (int j = 1; j < matrixHeight - 1; j++)
	{
		int seqPos = j - 1;
		
		mCurrent[0] = mPrevious[0] + mGapScore;
		
		seqScores[j] = mPrevious[0] + mGapScore;
		seqCounts[j] = 1;
		
		for (int i = 1; i < matrixLength; i++) 
		{
			int refPos = i - 1;
			
			int matchScore = mPrevious[i-1] + ((refEnd[-refPos - 1] == seqStart[seqPos]) ? mMatchScore : mMisMatchScore);
			int gapRefScore = mCurrent[i-1] + mGapScore;
			int gapReadScore = mPrevious[i] + mGapScore;
			int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
			
			if (bestScore > seqScores[j])
			{
				seqScores[j] = bestScore;
				seqCounts[j] = 1;
				refPositions[j] = refPos;
			}
			else if (bestScore == seqScores[j])
			{
				seqCounts[j]++;
			}
			
			mCurrent[i] = bestScore;
		}
		
		swap(mCurrent, mPrevious);
	}
	
	int seqPos = matrixHeight - 2;
	
	mCurrent[0] = mPrevious[0] + mGapScore;
	
	seqScores[matrixHeight - 1] = mPrevious[0] + mGapScore;
	seqCounts[matrixHeight - 1] = 1;
	
	for (int i = 1; i < matrixLength; i++) 
	{
		int refPos = i - 1;
		
		int matchScore = mPrevious[i-1] + ((refStart[-refPos] == seqStart[seqPos]) ? mMatchScore : mMisMatchScore);
		int gapRefScore = mCurrent[i-1] + mGapScore;
		int gapReadScore = mPrevious[i] + mGapScore;
		int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
		
		if (bestScore > seqScores[matrixHeight - 1])
		{
			seqScores[matrixHeight - 1] = bestScore;
			seqCounts[matrixHeight - 1] = 1;
			refPositions[matrixHeight - 1] = refPos;
		}
		else if (bestScore == seqScores[matrixHeight - 1])
		{
			seqCounts[matrixHeight - 1]++;
		}
		
		mCurrent[i] = bestScore;
	}
}

int BestSplitAlignment(const IntegerVec& scores1Fwd, const IntegerVec& scores2Rev, int breakInsertScore)
{
	IntegerVec scores1FwdBreakInsert(scores1Fwd.size());
	scores1FwdBreakInsert[0] = scores1Fwd[0];
	for (int idx = 1; idx < scores1Fwd.size(); idx++)
	{
		scores1FwdBreakInsert[idx] = max(scores1FwdBreakInsert[idx-1] + breakInsertScore, scores1Fwd[idx]);
	}
	
	int score = numeric_limits<int>::min();
	for (int idx1 = 0; idx1 < scores1FwdBreakInsert.size(); idx1++)
	{
		int idx2 = (int)scores2Rev.size() - idx1 - 1;
		score = max(score, scores1FwdBreakInsert[idx1] + scores2Rev[idx2]);
	}
	
	return score;
}

int BestSplitAlignment(const IntegerVec& scores1Fwd, const IntegerVec& scores2Rev, int breakInsertScore, IntegerVec& seq1Lengths, IntegerVec& seq2Lengths)
{
	IntegerVec scores1FwdBreakInsert(scores1Fwd.size());
	IntegerVec scores1FwdPrevMax(scores1Fwd.size());
	
	scores1FwdBreakInsert[0] = scores1Fwd[0];
	scores1FwdPrevMax[0] = 0;
	for (int idx = 1; idx < scores1Fwd.size(); idx++)
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
	for (int seq1Length = 0; seq1Length < scores1Fwd.size(); seq1Length++)
	{
		int seq2Length = scores2Rev.size() - seq1Length - 1;
		
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
		int seq2Length = scores2Rev.size() - seq1Length - 1;
		
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

bool Equal(const IntegerVec& a, const IntegerVec& b)
{
	return a.size() == b.size() && equal(a.begin(), a.end(), b.begin());
}

int testing()
{
	string refSeq1 = "XX12345XX";
	string refSeq2 = "YY5678Y5678Y";
	string readSeq = "12345678";
	
	string revReadSeq = readSeq;
	reverse(revReadSeq.begin(), revReadSeq.end());
	
	string revRefSeq1 = refSeq1;
	reverse(revRefSeq1.begin(), revRefSeq1.end());
	
	string revRefSeq2 = refSeq2;
	reverse(revRefSeq2.begin(), revRefSeq2.end());
	
	IntegerVec scores1Fwd;
	IntegerVec counts1Fwd;
	IntegerVec positions1Fwd;
	AlignPartial(refSeq1, readSeq, scores1Fwd, counts1Fwd, positions1Fwd);
	
	IntegerVec scores1Fwd2;
	IntegerVec counts1Fwd2;
	IntegerVec positions1Fwd2;
	AlignPartialFwd(&*refSeq1.begin(), &*refSeq1.end(), &*revReadSeq.begin(), &*revReadSeq.end(), scores1Fwd2, counts1Fwd2, positions1Fwd2);
	
	DebugCheck(Equal(scores1Fwd, scores1Fwd2));
	DebugCheck(Equal(counts1Fwd, counts1Fwd2));
	DebugCheck(Equal(positions1Fwd, positions1Fwd2));
	
	IntegerVec scores1Fwd3;
	IntegerVec counts1Fwd3;
	IntegerVec positions1Fwd3;
	AlignPartialRev(&*revRefSeq1.begin(), &*revRefSeq1.end(), &*readSeq.begin(), &*readSeq.end(), scores1Fwd3, counts1Fwd3, positions1Fwd3);
	
	DebugCheck(Equal(scores1Fwd, scores1Fwd3));
	DebugCheck(Equal(counts1Fwd, counts1Fwd3));
	DebugCheck(Equal(positions1Fwd, positions1Fwd3));
	
	IntegerVec scores2Rev;
	IntegerVec counts2Rev;
	IntegerVec positions2Rev;
	AlignPartial(revRefSeq2, revReadSeq, scores2Rev, counts2Rev, positions2Rev);
	
	IntegerVec seq1Length;
	IntegerVec seq2Length;
	int score1 = BestSplitAlignment(scores1Fwd, scores2Rev, -1, seq1Length, seq2Length);
	
	cout << score1 << endl;
	
	for (int i = 0; i < seq1Length.size(); i++)
	{
		cout << seq1Length[i] << "\t" << seq2Length[i] << endl;
		cout << scores1Fwd[seq1Length[i]] << "\t" << scores2Rev[seq2Length[i]] << endl;
		cout << counts1Fwd[seq1Length[i]] << "\t" << counts2Rev[seq2Length[i]] << endl;
		cout << positions1Fwd[seq1Length[i]] << "\t" << positions2Rev[seq2Length[i]] << endl;
		cout << readSeq.substr(0, seq1Length[i]) << "\t" << readSeq.substr(seq1Length[i], readSeq.size() - seq2Length[i] - seq1Length[i]) << "\t" << revReadSeq.substr(0, seq2Length[i]) << endl;
	}
}

int testing2()
{
	string refSeq = "XXXX12345678XXXXXX";
	string readSeq = "12345678";
	
	string revReadSeq = readSeq;
	reverse(revReadSeq.begin(), revReadSeq.end());
	
	IntegerVec scores1;
	IntegerVec counts1;
	IntegerVec positions1;
	AlignPartialRev(&*refSeq.begin(), &*refSeq.end(), &*revReadSeq.begin() + 4, &*revReadSeq.end(), scores1, counts1, positions1);
	
	int offset = (int)refSeq.size() - positions1.back() - 1;
	
	cout << offset << endl;
	cout << string(refSeq.begin() + offset, refSeq.end()) << endl;
	
	IntegerVec scores2;
	IntegerVec counts2;
	IntegerVec positions2;
	AlignPartialFwd(&*refSeq.begin() + offset, &*refSeq.end(), &*revReadSeq.begin(), &*revReadSeq.end(), scores2, counts2, positions2);
	
	Print(scores2);
	
	return 0;
}

const int cPadding = 16;

int main(int argc, char* argv[])
{
	string referenceFasta;
	string readSeqsFilename;
	string alignmentsFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Realignment Tool");
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqsFilenameArg("s","seq","Read Sequences",true,"","string",cmd);
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Sam Alignments",true,"","string",cmd);
		cmd.parse(argc,argv);
		
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
	
	Sequences referenceSequences;
	referenceSequences.Read(referenceFasta);
	
	cerr << "Reading fastq sequences" << endl;
	
	FastqReadStream readSeqsStream(readSeqsFilename);
	
	unordered_map<ReadID,string> readSequences[2];
	unordered_map<ReadID,string> rawReadSeqs;
	
	RawRead rawRead;
	while (readSeqsStream.GetNextRead(rawRead))
	{
		ReadID readID;
		readID.fragmentIndex = SAFEPARSE(int, rawRead.fragment);
		readID.readEnd = rawRead.readEnd;
		
		rawReadSeqs[readID] = rawRead.sequence;
		
		string paddedSequence = string(cPadding,'X') + rawRead.sequence + string(cPadding,'X');
		
		readSequences[0][readID] = paddedSequence;
		reverse(readSequences[0][readID].begin(), readSequences[0][readID].end());
		
		readSequences[1][readID] = paddedSequence;
		ReverseComplement(readSequences[1][readID]);
		reverse(readSequences[1][readID].begin(), readSequences[1][readID].end());
	}
	
	cerr << "Realigning" << endl;
	
	SamAlignmentStream alignmentStream(alignmentsFilename);
	FragmentAlignmentStream fragmentAlignmentStream(&alignmentStream);
	
	RawAlignmentVec alignments;
	while (fragmentAlignmentStream.GetNextAlignments(alignments))
	{
		string rawReadSeq[2];
		ReadID readID[2];
		const char* seqStartPtr[2][2];
		const char* seqEndPtr[2][2];
		
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			for (int strand = 0; strand <= 1; strand++)
			{
				readID[readEnd].fragmentIndex = SAFEPARSE(int, alignments.front().fragment);
				readID[readEnd].readEnd = readEnd;
				
				if (readSequences[strand].find(readID[readEnd]) == readSequences[strand].end())
				{
					cerr << "Error: Could not find sequence for read " << readID[readEnd].fragmentIndex << endl;
					exit(1);
				}
				
				seqStartPtr[strand][readEnd] = readSequences[strand][readID[readEnd]].c_str() + cPadding;
				seqEndPtr[strand][readEnd] = readSequences[strand][readID[readEnd]].c_str() + readSequences[strand][readID[readEnd]].length() - cPadding;
				rawReadSeq[readEnd] = rawReadSeqs[readID[readEnd]];
			}
		}
		
		RawAlignmentVec endAlignments[2];
		
		IntegerTable seqScoresSelf[2];
		IntegerTable seqCountsSelf[2];
		IntegerTable refPositionsSelf[2];
		
		IntegerTable seqScoresMate[2];
		IntegerTable seqCountsMate[2];
		IntegerTable refPositionsMate[2];
		
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			int readEnd = alignment.readEnd;
			
			seqScoresSelf[readEnd].push_back(IntegerVec());
			seqCountsSelf[readEnd].push_back(IntegerVec());
			refPositionsSelf[readEnd].push_back(IntegerVec());
			if (alignment.strand == PlusStrand)
			{
				const char* refPtr = referenceSequences.Get(alignment.reference, alignment.region.start);
				AlignPartialFwd(refPtr, refPtr + 200, seqStartPtr[0][readEnd], seqEndPtr[0][readEnd], seqScoresSelf[readEnd].back(), seqCountsSelf[readEnd].back(), refPositionsSelf[readEnd].back());
			}
			else
			{
				const char* refPtr = referenceSequences.Get(alignment.reference, alignment.region.end);
				AlignPartialRev(refPtr - 200, refPtr, seqStartPtr[1][readEnd], seqEndPtr[1][readEnd], seqScoresSelf[readEnd].back(), seqCountsSelf[readEnd].back(), refPositionsSelf[readEnd].back());
			}
			
			int mateEnd = OtherReadEnd(alignment.readEnd);
			
			seqScoresMate[readEnd].push_back(IntegerVec());
			seqCountsMate[readEnd].push_back(IntegerVec());
			refPositionsMate[readEnd].push_back(IntegerVec());
			if (alignment.strand == PlusStrand)
			{
				const char* refPtr = referenceSequences.Get(alignment.reference, alignment.region.start);
				AlignPartialFwd(refPtr, refPtr + 1000, seqStartPtr[1][mateEnd], seqEndPtr[1][mateEnd], seqScoresMate[readEnd].back(), seqCountsMate[readEnd].back(), refPositionsMate[readEnd].back());
			}
			else
			{
				const char* refPtr = referenceSequences.Get(alignment.reference, alignment.region.end);
				AlignPartialRev(refPtr - 1000, refPtr, seqStartPtr[0][mateEnd], seqEndPtr[0][mateEnd], seqScoresMate[readEnd].back(), seqCountsMate[readEnd].back(), refPositionsMate[readEnd].back());
			}
			
			endAlignments[readEnd].push_back(alignment);
		}
		
		if (endAlignments[0].size() == 0 || endAlignments[1].size() == 0)
		{
			continue;
		}
		
		int alignIndex[2];
		for (alignIndex[0] = 0; alignIndex[0] < endAlignments[0].size(); alignIndex[0]++)
		{
			for (alignIndex[1] = 0; alignIndex[1] < endAlignments[1].size(); alignIndex[1]++)
			{
				for (int readEnd = 0; readEnd <= 1; readEnd++)
				{
					int mateEnd = OtherReadEnd(readEnd);
					
					IntegerVec seq1Length;
					IntegerVec seq2Length;
					int score = BestSplitAlignment(seqScoresSelf[readEnd][alignIndex[readEnd]], seqScoresMate[mateEnd][alignIndex[mateEnd]], -1, seq1Length, seq2Length);
					
					cout << score << endl;
					
					string readSeq = rawReadSeq[readEnd];
					
					for (int i = 0; i < seq1Length.size(); i++)
					{
						cout << seq1Length[i] << "\t" << seq2Length[i] << endl;
						cout << seqScoresSelf[readEnd][alignIndex[readEnd]][seq1Length[i]] << "\t" << seqScoresMate[mateEnd][alignIndex[mateEnd]][seq2Length[i]] << endl;
						cout << seqCountsSelf[readEnd][alignIndex[readEnd]][seq1Length[i]] << "\t" << seqCountsMate[mateEnd][alignIndex[mateEnd]][seq2Length[i]] << endl;
						cout << refPositionsSelf[readEnd][alignIndex[readEnd]][seq1Length[i]] << "\t" << refPositionsMate[mateEnd][alignIndex[mateEnd]][seq2Length[i]] << endl;
						cout << readSeq.substr(0, seq1Length[i]) << "\t" << readSeq.substr(seq1Length[i], readSeq.size() - seq2Length[i] - seq1Length[i]) << "\t" << readSeq.substr(readSeq.length() - seq2Length[i]) << endl;
					}
				}
			}
		}
	}
}


