/*
 *  SimpleAligner.cpp
 *
 *  Created by Andrew McPherson on 05/09/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "SimpleAligner.h"
#include "DebugCheck.h"

#include <iostream>
#include <list>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>

using namespace std;


struct SSEConvert16
{
	SSEConvert16() {}
	
	SSEConvert16(short int a)
	{
		for (int i = 0; i < 8; i++)
		{
			array16[i] = a;
		}
	}
	
	union
	{
		__m128i sse;
		short int array16[8];
	};
	
	void Print()
	{
		for (int i = 0; i < 8; i++)
		{
			cout << array16[i] << "\t";
		}
		cout << endl;
	}
};

struct SSEConvert8
{
	SSEConvert8() {}
	
	SSEConvert8(unsigned char a)
	{
		for (int i = 0; i < 16; i++)
		{
			array8[i] = a;
		}
	}
	
	union
	{
		__m128i sse;
		unsigned char array8[16];
	};
	
	void Print()
	{
		for (int i = 0; i < 16; i++)
		{
			cout << (int)array8[i] << "\t";
		}
		cout << endl;
	}
};

void PrintSSE16(__m128i sse)
{
	SSEConvert16 csse;
	csse.sse = sse;
	csse.Print();
}

void PrintSSE8(__m128i sse)
{
	SSEConvert8 csse;
	csse.sse = sse;
	csse.Print();
}

SimpleAligner::SimpleAligner(int matchScore, int misMatchScore, int gapScore)
: mMatchScore(matchScore), mMisMatchScore(misMatchScore), mGapScore(gapScore), mScoreBuffer8(2000), mScoreBuffer16(2000)
{
	mIncrement16 = SSEConvert16(1).sse;
	
	const short int badScore16 = numeric_limits<short int>::min() / 2;
	const short int bestScore16 = numeric_limits<short int>::max();
	
	mMatchBonus16 = SSEConvert16(matchScore).sse;
	mMismatchPenalty16 = SSEConvert16(abs(misMatchScore)).sse;
	mGapPenalty16 = SSEConvert16(abs(gapScore)).sse;
	mMaxScoreInit16 = SSEConvert16(badScore16).sse;
	
	SSEConvert16 shiftRightSSEConvert16((short int)0);
	shiftRightSSEConvert16.array16[7] = badScore16;
	mShiftRight16 = shiftRightSSEConvert16.sse;
	
	SSEConvert16 shiftLeftSSEConvert16((short int)0);
	shiftLeftSSEConvert16.array16[0] = badScore16;
	mShiftLeft16 = shiftLeftSSEConvert16.sse;
	
	SSEConvert16 scoreOddInit16((short int)badScore16);
	scoreOddInit16.array16[3] = 0;
	mScoreOddInit16 = scoreOddInit16.sse;
	
	SSEConvert16 scoreEvenInit16((short int)badScore16);
	scoreEvenInit16.array16[3] = (short int)gapScore;
	scoreEvenInit16.array16[4] = (short int)gapScore;
	mScoreEvenInit16 = scoreEvenInit16.sse;
	
	SSEConvert16 clampCornerScoreOdd16(bestScore16);
	clampCornerScoreOdd16.array16[7] = badScore16;
	mClampCornerScoreOdd16 = clampCornerScoreOdd16.sse;
	
	// Initialization for AlignBandedSSE2BW7ScoreFwd
	{
		SSEConvert16 lengthInit(0);
		for (int i = 0; i < 8; i++)
		{
			lengthInit.array16[i] = i - 3;
		}
		
		__m128i curLen = lengthInit.sse;
		
		__m128i refLen = SSEConvert16(0).sse;
		
		__m128i scoreOdd = mScoreOddInit16;
		__m128i scoreEven = mScoreEvenInit16;
		__m128i maxScore = mMaxScoreInit16;
		
		__m128i scoreCmpOdd = _mm_cmpgt_epi16(scoreOdd,maxScore);
		refLen = _mm_or_si128(_mm_and_si128(scoreCmpOdd,curLen),_mm_andnot_si128(scoreCmpOdd,refLen));
		
		maxScore = _mm_max_epi16(maxScore,scoreOdd);
		
		maxScore = _mm_slli_si128(maxScore,2);
		maxScore = _mm_or_si128(maxScore,mShiftLeft16);
		
		refLen = _mm_slli_si128(refLen,2);
		
		__m128i scoreCmpEven = _mm_cmpgt_epi16(scoreEven,maxScore);
		refLen = _mm_or_si128(_mm_and_si128(scoreCmpEven,curLen),_mm_andnot_si128(scoreCmpEven,refLen));
		
		maxScore = _mm_max_epi16(maxScore,scoreEven);
		
		curLen = _mm_add_epi16(curLen,mIncrement16);
		
		mCurLenFwdInit16 = curLen;
		mRefLenFwdInit16 = refLen;
		mMaxScoreFwdInit16 = maxScore;
	}
	
	// Initialization for AlignBandedSSE2BW7ScoreRev
	{
		SSEConvert16 lengthInit(0);
		for (int i = 0; i < 8; i++)
		{
			lengthInit.array16[i] = 3 - i;
		}
		
		__m128i curLen = lengthInit.sse;
		
		__m128i refLen = SSEConvert16(0).sse;
		
		__m128i scoreOdd = mScoreOddInit16;
		__m128i scoreEven = mScoreEvenInit16;
		__m128i maxScore = mMaxScoreInit16;
		
		maxScore = _mm_srli_si128(maxScore,2);
		maxScore = _mm_or_si128(maxScore,mShiftRight16);
		
		__m128i scoreCmpOdd = _mm_cmpgt_epi16(scoreOdd,maxScore);
		refLen = _mm_or_si128(_mm_and_si128(scoreCmpOdd,curLen),_mm_andnot_si128(scoreCmpOdd,refLen));
		
		maxScore = _mm_max_epi16(maxScore,scoreOdd);
		
		curLen = _mm_add_epi16(curLen,mIncrement16);
		
		__m128i scoreCmpEven = _mm_cmpgt_epi16(scoreEven,maxScore);
		refLen = _mm_or_si128(_mm_and_si128(scoreCmpEven,curLen),_mm_andnot_si128(scoreCmpEven,refLen));
		
		maxScore = _mm_max_epi16(maxScore,scoreEven);
		
		mCurLenRevInit16 = curLen;
		mRefLenRevInit16 = refLen;
		mMaxScoreRevInit16 = maxScore;
	}
	
	mOffset8 = 40;
	mZero8 = 80;
	
	//scores
	mMatchBonus8 = SSEConvert8(matchScore).sse;
	mMismatchPenalty8 = SSEConvert8(abs(misMatchScore)).sse;
	mGapPenalty8 = SSEConvert8(abs(gapScore)).sse;
	
	// Initialization for Align16baseSSE2Fwd
	{
		SSEConvert8 shiftRefGapSSEConvert8(0);
		shiftRefGapSSEConvert8.array8[15] = mOffset8 + mZero8;
		mShiftRefGapFwd8 = shiftRefGapSSEConvert8.sse;
		
		SSEConvert8 scoreDiagonalInit8(mOffset8);
		scoreDiagonalInit8.array8[15] = mZero8 + mOffset8;
		mScoreDiagonalInitFwd8 = scoreDiagonalInit8.sse;
		
		SSEConvert8 scoreOrthogonalInit8(mOffset8);
		scoreOrthogonalInit8.array8[15] = mZero8 - abs(gapScore) + mOffset8;
		mScoreOrthogonalInitFwd8 = scoreOrthogonalInit8.sse;
	}
	
	// Initialization for Align16baseSSE2Rev
	{
		SSEConvert8 shiftRefGapSSEConvert8(0);
		shiftRefGapSSEConvert8.array8[0] = mOffset8 + mZero8;
		mShiftRefGapRev8 = shiftRefGapSSEConvert8.sse;
		
		SSEConvert8 scoreDiagonalInit8(mOffset8);
		scoreDiagonalInit8.array8[0] = mZero8 + mOffset8;
		mScoreDiagonalInitRev8 = scoreDiagonalInit8.sse;
		
		SSEConvert8 scoreOrthogonalInit8(mOffset8);
		scoreOrthogonalInit8.array8[0] = mZero8 - abs(gapScore) + mOffset8;
		mScoreOrthogonalInitRev8 = scoreOrthogonalInit8.sse;
	}
}

void SimpleAligner::Align16baseSSE2Fwd(const char* refStart, const char* refEnd, const char* seqEnd, short int* refScores)
{
	__m128i scoreDiagonal = mScoreDiagonalInitFwd8;
	__m128i scoreOrthogonal = mScoreOrthogonalInitFwd8;
	
	const char* fwdPtr = refStart - 15;
	const char* revPtr = seqEnd;
	
	int refScoreIndex = -14;
	
	while ( fwdPtr <= refEnd - 1 )
	{
		__m128i fwd = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i rev = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmp8 = _mm_cmpeq_epi8(fwd,rev);
		
		// diagonalScore + match/mismatch
		__m128i matchDiagonal = _mm_sub_epi8(_mm_add_epi8(scoreDiagonal,_mm_and_si128(cmp8,mMatchBonus8)),_mm_andnot_si128(cmp8,mMismatchPenalty8));
		
		// left score + gap
		__m128i gapLeft = _mm_sub_epi8(scoreOrthogonal, mGapPenalty8);
		
		// up score + gap
		__m128i scoreOrthogonalShift = _mm_or_si128( _mm_srli_si128(scoreOrthogonal, 1), mShiftRefGapFwd8);
		__m128i gapUp = _mm_sub_epi8( scoreOrthogonalShift, mGapPenalty8);
		__m128i nextScores = _mm_max_epu8(_mm_max_epu8(gapLeft, gapUp), matchDiagonal);
		
		// max
		scoreDiagonal = scoreOrthogonalShift;
		scoreOrthogonal = nextScores;
		
		fwdPtr++;
		
		refScores[max(0,refScoreIndex)] = _mm_extract_epi16(_mm_srli_si128(_mm_slli_si128(nextScores, 15), 15), 0) - mZero8 - mOffset8;
		refScoreIndex++;
	}
}

void SimpleAligner::Align16baseSSE2Rev(const char* refStart, const char* refEnd, const char* seqStart, short int* refScores)
{
	__m128i scoreDiagonal = mScoreDiagonalInitRev8;
	__m128i scoreOrthogonal = mScoreOrthogonalInitRev8;
	
	const char* fwdPtr = seqStart;
	const char* revPtr = refEnd - 1;
	
	int refScoreIndex = -14;
	
	while ( revPtr >= refStart - 15 )
	{
		__m128i fwd = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i rev = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmp8 = _mm_cmpeq_epi8(fwd,rev);
		
		// diagonalScore + match/mismatch
		__m128i matchDiagonal = _mm_sub_epi8(_mm_add_epi8(scoreDiagonal,_mm_and_si128(cmp8,mMatchBonus8)),_mm_andnot_si128(cmp8,mMismatchPenalty8));
		
		// left score + gap
		__m128i gapLeft = _mm_sub_epi8(scoreOrthogonal, mGapPenalty8);
		
		// up score + gap
		__m128i scoreOrthogonalShift = _mm_or_si128( _mm_slli_si128(scoreOrthogonal, 1), mShiftRefGapRev8);
		__m128i gapUp = _mm_sub_epi8( scoreOrthogonalShift, mGapPenalty8);
		
		__m128i nextScores = _mm_max_epu8(_mm_max_epu8(gapLeft, gapUp), matchDiagonal);
		
		// max
		scoreDiagonal = scoreOrthogonalShift;
		scoreOrthogonal = nextScores;
		
		revPtr--;
		
		refScores[max(0,refScoreIndex)] = _mm_extract_epi16(_mm_srli_si128(nextScores, 15), 0) - mZero8 - mOffset8;
		refScoreIndex++;
	}
}

int SimpleAligner::AlignEndToEnd(const string& reference, const string& sequence)
{
	if (sequence.size() == 0)
	{
		return 0;
	}
	
	int matrixLength = reference.size() + 1;
	int matrixHeight = sequence.size() + 1;
	
	if (matrixLength != mCurrent.size())
	{
		mCurrent.resize(matrixLength);
		mPrevious.resize(matrixLength);
	}
	
	for (int i = 0; i < matrixLength; i++) 
	{
		mPrevious[i] = 0;
	}
	
	for (int j = 1; j < matrixHeight - 1; j++)
	{
		int seqPos = j - 1;
		
		mCurrent[0] = mPrevious[0] + mGapScore;
		
		for (int i = 1; i < matrixLength; i++) 
		{
			int refPos = i - 1;
			
			int matchScore = mPrevious[i-1] + ((reference[refPos] == sequence[seqPos]) ? mMatchScore : mMisMatchScore);
			int gapRefScore = mCurrent[i-1] + mGapScore;
			int gapReadScore = mPrevious[i] + mGapScore;
			int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
			
			mCurrent[i] = bestScore;
		}
		
		swap(mCurrent, mPrevious);
	}
	
	int seqPos = matrixHeight - 2;
	
	int score = numeric_limits<int>::min();
	
	mCurrent[0] = mPrevious[0] + mGapScore;
	
	for (int i = 1; i < matrixLength; i++) 
	{
		int refPos = i - 1;
		
		int matchScore = mPrevious[i-1] + ((reference[refPos] == sequence[seqPos]) ? mMatchScore : mMisMatchScore);
		int gapRefScore = mCurrent[i-1] + mGapScore;
		int gapReadScore = mPrevious[i] + mGapScore;
		int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
		
		score = max(score, bestScore);
		
		mCurrent[i] = bestScore;
	}
	
	return score;
}

void SimpleAligner::AlignEndToEndFwd(const char* refStart, const char* refEnd, const char* seqStart, const char* seqEnd, IntegerVec& refLengthScores)
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
	
	refLengthScores.resize(matrixLength);
	
	for (int i = 0; i < matrixLength; i++) 
	{
		mPrevious[i] = 0;
	}
	
	for (int j = 1; j < matrixHeight - 1; j++)
	{
		int seqPos = j - 1;
		
		mCurrent[0] = mPrevious[0] + mGapScore;
		
		for (int i = 1; i < matrixLength; i++) 
		{
			int refPos = i - 1;
			
			int matchScore = mPrevious[i-1] + ((refStart[refPos] == seqEnd[-seqPos - 1]) ? mMatchScore : mMisMatchScore);
			int gapRefScore = mCurrent[i-1] + mGapScore;
			int gapReadScore = mPrevious[i] + mGapScore;
			int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
			
			mCurrent[i] = bestScore;
		}
		
		swap(mCurrent, mPrevious);
	}
	
	int seqPos = matrixHeight - 2;
	
	mCurrent[0] = mPrevious[0] + mGapScore;
	
	refLengthScores[0] = mPrevious[0] + mGapScore;
	
	for (int i = 1; i < matrixLength; i++) 
	{
		int refPos = i - 1;
		
		int matchScore = mPrevious[i-1] + ((refStart[refPos] == seqEnd[-seqPos - 1]) ? mMatchScore : mMisMatchScore);
		int gapRefScore = mCurrent[i-1] + mGapScore;
		int gapReadScore = mPrevious[i] + mGapScore;
		int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
		
		refLengthScores[i] = bestScore;
		
		mCurrent[i] = bestScore;
	}
}

void SimpleAligner::AlignEndToEndRev(const char* refStart, const char* refEnd, const char* seqStart, const char* seqEnd, IntegerVec& refLengthScores)
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
	
	refLengthScores.resize(matrixLength);
	
	for (int i = 0; i < matrixLength; i++) 
	{
		mPrevious[i] = 0;
	}
	
	for (int j = 1; j < matrixHeight - 1; j++)
	{
		int seqPos = j - 1;
		
		mCurrent[0] = mPrevious[0] + mGapScore;
		
		for (int i = 1; i < matrixLength; i++) 
		{
			int refPos = i - 1;
			
			int matchScore = mPrevious[i-1] + ((refEnd[-refPos - 1] == seqStart[seqPos]) ? mMatchScore : mMisMatchScore);
			int gapRefScore = mCurrent[i-1] + mGapScore;
			int gapReadScore = mPrevious[i] + mGapScore;
			int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
			
			mCurrent[i] = bestScore;
		}
		
		swap(mCurrent, mPrevious);
	}
	
	int seqPos = matrixHeight - 2;
	
	mCurrent[0] = mPrevious[0] + mGapScore;
	
	refLengthScores[0] = mPrevious[0] + mGapScore;
	
	for (int i = 1; i < matrixLength; i++) 
	{
		int refPos = i - 1;
		
		int matchScore = mPrevious[i-1] + ((refEnd[-refPos - 1] == seqStart[seqPos]) ? mMatchScore : mMisMatchScore);
		int gapRefScore = mCurrent[i-1] + mGapScore;
		int gapReadScore = mPrevious[i] + mGapScore;
		int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
		
		refLengthScores[i] = bestScore;
		
		mCurrent[i] = bestScore;
	}
}

void SimpleAligner::AlignPartial(const string& reference, const string& sequence, IntegerVec& seqLengthScores)
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
	
	seqLengthScores.resize(sequence.size() + 1);
	
	for (int i = 0; i < matrixLength; i++) 
	{
		mPrevious[i] = 0;
	}

	seqLengthScores[0] = 0;
	
	for (int j = 1; j < matrixHeight - 1; j++)
	{
		int seqPos = j - 1;
		
		mCurrent[0] = mPrevious[0] + mGapScore;
		
		int seqLengthScore = mPrevious[0] + mGapScore;
		
		for (int i = 1; i < matrixLength; i++) 
		{
			int refPos = i - 1;
			
			int matchScore = mPrevious[i-1] + ((reference[refPos] == sequence[seqPos]) ? mMatchScore : mMisMatchScore);
			int gapRefScore = mCurrent[i-1] + mGapScore;
			int gapReadScore = mPrevious[i] + mGapScore;
			int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
			
			seqLengthScore = max(seqLengthScore, bestScore);
			
			mCurrent[i] = bestScore;
		}
		
		seqLengthScores[j] = seqLengthScore;
		
		swap(mCurrent, mPrevious);
	}
	
	int seqPos = matrixHeight - 2;
	
	mCurrent[0] = mPrevious[0] + mGapScore;
	
	int seqLengthScore = mPrevious[0] + mGapScore;
	
	for (int i = 1; i < matrixLength; i++) 
	{
		int refPos = i - 1;
		
		int matchScore = mPrevious[i-1] + ((reference[refPos] == sequence[seqPos]) ? mMatchScore : mMisMatchScore);
		int gapRefScore = mCurrent[i-1] + mGapScore;
		int gapReadScore = mPrevious[i] + mGapScore;
		int bestScore = max(matchScore,max(gapRefScore,gapReadScore));
		
		seqLengthScore = max(seqLengthScore, bestScore);
		
		mCurrent[i] = bestScore;
	}
	
	seqLengthScores[sequence.size()] = seqLengthScore;
}

struct matprint
{
	matprint()
	{
		diag = 1;
	}
	
	void add(int a, int b, int score)
	{
		matrix[pair<int,int>(a,b)] = score;
	}
	
	void next16(__m128i nextScoreOdd, __m128i nextScoreEven)
	{
		SSEConvert16 nodd;
		nodd.sse = nextScoreOdd;
		SSEConvert16 neven;
		neven.sse = nextScoreEven;
		
		for (int i = 0; i <= 7; i++)
		{
			int a = diag - i + 3;
			int b = diag + i - 3;
			if (a < 0 || b < 0)
			{
				continue;
			}
			//			cout << a << "\t" << b << "\t" << nodd.array[i] - 0x30 << endl;
			matrix[pair<int,int>(a,b)] = nodd.array16[i];
		}
		
		for (int i = 0; i <= 7; i++)
		{
			int a = diag - i + 4;
			int b = diag + i - 3;
			if (a < 0 || b < 0)
			{
				continue;
			}
			//			cout << a << "\t" << b << "\t" << nodd.array[i] - 0x30 << endl;
			matrix[pair<int,int>(a,b)] = neven.array16[i];
		}
		
		diag++;
	}
	
	void printfwd(int bandWidth, const char* refStart, const char* seqStart, const char* seqEnd)
	{
		int matrixLength = seqEnd - seqStart + 1 + bandWidth;
		int matrixHeight = seqEnd - seqStart + 1;
		
		cout << "\t\t";
		for (int j = 1; j < matrixLength; j++)
		{
			cout << refStart[j-1] << "\t";
		}
		cout << endl;
		
		for (int i = 0; i < matrixHeight; i++)
		{
			if (i == 0)
			{
				cout << "\t";
			}
			else
			{
				cout << seqEnd[-i] << "\t";
			}
			
			for (int j = 0; j < matrixLength; j++)
			{
				if (matrix.find(pair<int,int>(j,i)) == matrix.end())
				{
					cout << " \t";
				}
				else
				{
					cout << matrix.find(pair<int,int>(j,i))->second << "\t";
				}
			}
			
			cout << endl;
		}
	}
	
	void print(int maxsize)
	{
		for (int i = 0; i < maxsize; i++)
		{
			for (int j = 0; j < maxsize; j++)
			{
				if (matrix.find(pair<int,int>(i,j)) == matrix.end())
				{
					cout << " \t";
				}
				else
				{
					cout << matrix.find(pair<int,int>(i,j))->second << "\t";
				}
			}
			
			cout << endl;
		}
	}
	
	unordered_map<pair<int,int>,int> matrix;
	int diag;
};

int SimpleAligner::AlignBandedSSE2BW7ScoreFwd(const char* refStart, const char* seqStart, const char* seqEnd)
{
	__m128i scoreOdd = mScoreOddInit16;
	__m128i scoreEven = mScoreEvenInit16;
	__m128i maxScore = {0,0};
	
	const char* fwdPtr = refStart - 3;
	const char* revPtr = seqEnd - 4;
	
	while (revPtr >= seqStart - 7)
	{
		__m128i fwdOdd = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i revOdd = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmpOdd8 = _mm_cmpeq_epi8(fwdOdd,revOdd);
		__m128i cmpOdd = _mm_unpacklo_epi8(cmpOdd8,cmpOdd8);
		
		// matchOdd = scoreOdd + cmpOdd ? mMatchBonus16 : -mMismatchPenalty16
		__m128i matchOdd = _mm_sub_epi16(_mm_add_epi16(scoreOdd,_mm_and_si128(cmpOdd,mMatchBonus16)),_mm_andnot_si128(cmpOdd,mMismatchPenalty16));
		
		// gapOdd = scoreEven - gapPenalty
		__m128i gapOdd = _mm_sub_epi16(scoreEven,mGapPenalty16);
		
		// max
		__m128i nextScoreOdd = _mm_min_epi16(_mm_max_epi16(matchOdd,_mm_max_epi16(_mm_or_si128(_mm_srli_si128(gapOdd,2),mShiftRight16),gapOdd)),mClampCornerScoreOdd16);
		
		revPtr--;
		
		__m128i fwdEven = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i revEven = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmpEven8 = _mm_cmpeq_epi8(fwdEven,revEven);
		__m128i cmpEven = _mm_unpacklo_epi8(cmpEven8,cmpEven8);
		
		// matchEven = scoreEven + cmpEven ? mMatchBonus16 : -mMismatchPenalty16
		__m128i matchEven = _mm_sub_epi16(_mm_add_epi16(scoreEven,_mm_and_si128(cmpEven,mMatchBonus16)),_mm_andnot_si128(cmpEven,mMismatchPenalty16));
		
		// gapEven = scoreOdd - mGapPenalty16
		__m128i gapEven = _mm_sub_epi16(nextScoreOdd,mGapPenalty16);
		
		// max
		__m128i nextScoreEven = _mm_max_epi16(matchEven,_mm_max_epi16(_mm_or_si128(_mm_slli_si128(gapEven,2),mShiftLeft16),gapEven));
		
		fwdPtr++;
		
		scoreOdd = nextScoreOdd;
		scoreEven = nextScoreEven;
		
		maxScore = _mm_max_epi16(maxScore,nextScoreOdd);
		maxScore = _mm_max_epi16(maxScore,nextScoreEven);
	}
	
	maxScore = _mm_max_epi16(_mm_srli_si128(maxScore,8),maxScore);
	maxScore = _mm_max_epi16(_mm_srli_si128(maxScore,4),maxScore);
	maxScore = _mm_max_epi16(_mm_srli_si128(maxScore,2),maxScore);
	
	return _mm_extract_epi16(maxScore,0);
}

void SimpleAligner::AlignBandedSSE2BW7ScoreFwd(const char* refStart, const char* seqStart, const char* seqEnd, short int* seqScores)
{
	__m128i scoreOdd = mScoreOddInit16;
	__m128i scoreEven = mScoreEvenInit16;
	__m128i maxScore = mMaxScoreFwdInit16;
	
	const char* fwdPtr = refStart - 3;
	const char* revPtr = seqEnd - 4;
	
	int seqScoreIndex = -2;
	while (revPtr >= seqStart - 7)
	{
		__m128i fwdOdd = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i revOdd = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmpOdd8 = _mm_cmpeq_epi8(fwdOdd,revOdd);
		__m128i cmpOdd = _mm_unpacklo_epi8(cmpOdd8,cmpOdd8);
		
		// matchOdd = scoreOdd + cmpOdd ? mMatchBonus16 : -mMismatchPenalty16
		__m128i matchOdd = _mm_sub_epi16(_mm_add_epi16(scoreOdd,_mm_and_si128(cmpOdd,mMatchBonus16)),_mm_andnot_si128(cmpOdd,mMismatchPenalty16));
		
		// gapOdd = scoreEven - mGapPenalty16
		__m128i gapOdd = _mm_sub_epi16(scoreEven,mGapPenalty16);
		
		// max
		__m128i nextScoreOdd = _mm_min_epi16(_mm_max_epi16(matchOdd,_mm_max_epi16(_mm_or_si128(_mm_srli_si128(gapOdd,2),mShiftRight16),gapOdd)),mClampCornerScoreOdd16);
		
		revPtr--;
		
		__m128i fwdEven = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i revEven = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmpEven8 = _mm_cmpeq_epi8(fwdEven,revEven);
		__m128i cmpEven = _mm_unpacklo_epi8(cmpEven8,cmpEven8);
		
		// matchEven = scoreEven + cmpEven ? mMatchBonus16 : -mMismatchPenalty16
		__m128i matchEven = _mm_sub_epi16(_mm_add_epi16(scoreEven,_mm_and_si128(cmpEven,mMatchBonus16)),_mm_andnot_si128(cmpEven,mMismatchPenalty16));
		
		// gapEven = scoreOdd - mGapPenalty16
		__m128i gapEven = _mm_sub_epi16(nextScoreOdd,mGapPenalty16);
		
		// max
		__m128i nextScoreEven = _mm_max_epi16(matchEven,_mm_max_epi16(_mm_or_si128(_mm_slli_si128(gapEven,2),mShiftLeft16),gapEven));
		
		fwdPtr++;
		
		scoreOdd = nextScoreOdd;
		scoreEven = nextScoreEven;
		
		maxScore = _mm_max_epi16(maxScore,scoreOdd);
		maxScore = _mm_slli_si128(maxScore,2);
		maxScore = _mm_or_si128(maxScore,mShiftLeft16);
		maxScore = _mm_max_epi16(maxScore,scoreEven);
		
		seqScores[max(0,seqScoreIndex)] = _mm_extract_epi16(maxScore,7);
		seqScoreIndex++;
	}
}

void SimpleAligner::AlignBandedSSE2BW7ScoreFwd(const char* refStart, const char* seqStart, const char* seqEnd, short int* seqScores, short int* refLengths)
{
	__m128i scoreOdd = mScoreOddInit16;
	__m128i scoreEven = mScoreEvenInit16;
	__m128i maxScore = mMaxScoreFwdInit16;
	__m128i curLen = mCurLenFwdInit16;
	__m128i refLen = mRefLenFwdInit16;
	
	const char* fwdPtr = refStart - 3;
	const char* revPtr = seqEnd - 4;
	
	int seqScoreIndex = -2;
	while (revPtr >= seqStart - 7)
	{
		__m128i fwdOdd = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i revOdd = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmpOdd8 = _mm_cmpeq_epi8(fwdOdd,revOdd);
		__m128i cmpOdd = _mm_unpacklo_epi8(cmpOdd8,cmpOdd8);
		
		// matchOdd = scoreOdd + cmpOdd ? mMatchBonus16 : -mMismatchPenalty16
		__m128i matchOdd = _mm_sub_epi16(_mm_add_epi16(scoreOdd,_mm_and_si128(cmpOdd,mMatchBonus16)),_mm_andnot_si128(cmpOdd,mMismatchPenalty16));
		
		// gapOdd = scoreEven - mGapPenalty16
		__m128i gapOdd = _mm_sub_epi16(scoreEven,mGapPenalty16);
		
		// max
		__m128i nextScoreOdd = _mm_min_epi16(_mm_max_epi16(matchOdd,_mm_max_epi16(_mm_or_si128(_mm_srli_si128(gapOdd,2),mShiftRight16),gapOdd)),mClampCornerScoreOdd16);
		
		revPtr--;
		
		__m128i fwdEven = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i revEven = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmpEven8 = _mm_cmpeq_epi8(fwdEven,revEven);
		__m128i cmpEven = _mm_unpacklo_epi8(cmpEven8,cmpEven8);
		
		// matchEven = scoreEven + cmpEven ? mMatchBonus16 : -mMismatchPenalty16
		__m128i matchEven = _mm_sub_epi16(_mm_add_epi16(scoreEven,_mm_and_si128(cmpEven,mMatchBonus16)),_mm_andnot_si128(cmpEven,mMismatchPenalty16));
		
		// gapEven = scoreOdd - mGapPenalty16
		__m128i gapEven = _mm_sub_epi16(nextScoreOdd,mGapPenalty16);
		
		// max
		__m128i nextScoreEven = _mm_max_epi16(matchEven,_mm_max_epi16(_mm_or_si128(_mm_slli_si128(gapEven,2),mShiftLeft16),gapEven));
		
		fwdPtr++;
		
		scoreOdd = nextScoreOdd;
		scoreEven = nextScoreEven;
		
		__m128i scoreCmpOdd = _mm_cmpgt_epi16(scoreOdd,maxScore);
		refLen = _mm_or_si128(_mm_and_si128(scoreCmpOdd,curLen),_mm_andnot_si128(scoreCmpOdd,refLen));
		
		maxScore = _mm_max_epi16(maxScore,scoreOdd);
		
		maxScore = _mm_slli_si128(maxScore,2);
		maxScore = _mm_or_si128(maxScore,mShiftLeft16);
		
		refLen = _mm_slli_si128(refLen,2);
		
		__m128i scoreCmpEven = _mm_cmpgt_epi16(scoreEven,maxScore);
		refLen = _mm_or_si128(_mm_and_si128(scoreCmpEven,curLen),_mm_andnot_si128(scoreCmpEven,refLen));
		
		maxScore = _mm_max_epi16(maxScore,scoreEven);
		
		curLen = _mm_add_epi16(curLen,mIncrement16);
		
		seqScores[max(0,seqScoreIndex)] = _mm_extract_epi16(maxScore,7);
		refLengths[max(0,seqScoreIndex)] = _mm_extract_epi16(refLen,7);
		seqScoreIndex++;
	}
}

int SimpleAligner::AlignBandedSSE2BW7ScoreRev(const char* refStart, const char* seqStart, const char* seqEnd)
{
	__m128i scoreOdd = mScoreOddInit16;
	__m128i scoreEven = mScoreEvenInit16;
	__m128i maxScore = {0,0};
	
	const char* fwdPtr = seqStart - 3;
	const char* revPtr = refStart - 3;
	
	while (fwdPtr <= seqEnd - 1)
	{
		__m128i fwdOdd = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i revOdd = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmpOdd8 = _mm_cmpeq_epi8(fwdOdd,revOdd);
		__m128i cmpOdd = _mm_unpacklo_epi8(cmpOdd8,cmpOdd8);
		
		// matchOdd = scoreOdd + cmpOdd ? mMatchBonus16 : -mMismatchPenalty16
		__m128i matchOdd = _mm_sub_epi16(_mm_add_epi16(scoreOdd,_mm_and_si128(cmpOdd,mMatchBonus16)),_mm_andnot_si128(cmpOdd,mMismatchPenalty16));
		
		// gapOdd = scoreEven - mGapPenalty16
		__m128i gapOdd = _mm_sub_epi16(scoreEven,mGapPenalty16);
		
		// max
		__m128i nextScoreOdd = _mm_min_epi16(_mm_max_epi16(matchOdd,_mm_max_epi16(_mm_or_si128(_mm_srli_si128(gapOdd,2),mShiftRight16),gapOdd)),mClampCornerScoreOdd16);
		
		revPtr--;
		
		__m128i fwdEven = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i revEven = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmpEven8 = _mm_cmpeq_epi8(fwdEven,revEven);
		__m128i cmpEven = _mm_unpacklo_epi8(cmpEven8,cmpEven8);
		
		// matchEven = scoreEven + cmpEven ? mMatchBonus16 : -mMismatchPenalty16
		__m128i matchEven = _mm_sub_epi16(_mm_add_epi16(scoreEven,_mm_and_si128(cmpEven,mMatchBonus16)),_mm_andnot_si128(cmpEven,mMismatchPenalty16));
		
		// gapEven = scoreOdd - mGapPenalty16
		__m128i gapEven = _mm_sub_epi16(nextScoreOdd,mGapPenalty16);
		
		// max
		__m128i nextScoreEven = _mm_max_epi16(matchEven,_mm_max_epi16(_mm_or_si128(_mm_slli_si128(gapEven,2),mShiftLeft16),gapEven));
		
		fwdPtr++;
		
		scoreOdd = nextScoreOdd;
		scoreEven = nextScoreEven;
		
		maxScore = _mm_max_epi16(maxScore,nextScoreOdd);
		maxScore = _mm_max_epi16(maxScore,nextScoreEven);
	}
	
	maxScore = _mm_max_epi16(_mm_srli_si128(maxScore,8),maxScore);
	maxScore = _mm_max_epi16(_mm_srli_si128(maxScore,4),maxScore);
	maxScore = _mm_max_epi16(_mm_srli_si128(maxScore,2),maxScore);
	
	return _mm_extract_epi16(maxScore,0);
}

void SimpleAligner::AlignBandedSSE2BW7ScoreRev(const char* refStart, const char* seqStart, const char* seqEnd, short int* seqScores)
{
	__m128i scoreOdd = mScoreOddInit16;
	__m128i scoreEven = mScoreEvenInit16;
	__m128i maxScore = mMaxScoreRevInit16;
	
	const char* fwdPtr = seqStart - 3;
	const char* revPtr = refStart - 3;
	
	int seqScoreIndex = -2;
	while (fwdPtr <= seqEnd - 1)
	{
		__m128i fwdOdd = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i revOdd = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmpOdd8 = _mm_cmpeq_epi8(fwdOdd,revOdd);
		__m128i cmpOdd = _mm_unpacklo_epi8(cmpOdd8,cmpOdd8);
		
		// matchOdd = scoreOdd + cmpOdd ? mMatchBonus16 : -mMismatchPenalty16
		__m128i matchOdd = _mm_sub_epi16(_mm_add_epi16(scoreOdd,_mm_and_si128(cmpOdd,mMatchBonus16)),_mm_andnot_si128(cmpOdd,mMismatchPenalty16));
		
		// gapOdd = scoreEven - mGapPenalty16
		__m128i gapOdd = _mm_sub_epi16(scoreEven,mGapPenalty16);
		
		// max
		__m128i nextScoreOdd = _mm_min_epi16(_mm_max_epi16(matchOdd,_mm_max_epi16(_mm_or_si128(_mm_srli_si128(gapOdd,2),mShiftRight16),gapOdd)),mClampCornerScoreOdd16);
		
		revPtr--;
		
		__m128i fwdEven = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i revEven = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmpEven8 = _mm_cmpeq_epi8(fwdEven,revEven);
		__m128i cmpEven = _mm_unpacklo_epi8(cmpEven8,cmpEven8);
		
		// matchEven = scoreEven + cmpEven ? mMatchBonus16 : -mMismatchPenalty16
		__m128i matchEven = _mm_sub_epi16(_mm_add_epi16(scoreEven,_mm_and_si128(cmpEven,mMatchBonus16)),_mm_andnot_si128(cmpEven,mMismatchPenalty16));
		
		// gapEven = scoreOdd - mGapPenalty16
		__m128i gapEven = _mm_sub_epi16(nextScoreOdd,mGapPenalty16);
		
		// max
		__m128i nextScoreEven = _mm_max_epi16(matchEven,_mm_max_epi16(_mm_or_si128(_mm_slli_si128(gapEven,2),mShiftLeft16),gapEven));
		
		fwdPtr++;
		
		scoreOdd = nextScoreOdd;
		scoreEven = nextScoreEven;
		
		maxScore = _mm_srli_si128(maxScore,2);
		maxScore = _mm_or_si128(maxScore,mShiftRight16);
		maxScore = _mm_max_epi16(maxScore,scoreOdd);
		maxScore = _mm_max_epi16(maxScore,scoreEven);
		
		seqScores[max(0,seqScoreIndex)] = _mm_extract_epi16(maxScore,0);
		seqScoreIndex++;
	}
}

void SimpleAligner::AlignBandedSSE2BW7ScoreRev(const char* refStart, const char* seqStart, const char* seqEnd, short int* seqScores, short int* refLengths)
{
	__m128i scoreOdd = mScoreOddInit16;
	__m128i scoreEven = mScoreEvenInit16;
	__m128i maxScore = mMaxScoreRevInit16;
	__m128i curLen = mCurLenRevInit16;
	__m128i refLen = mRefLenRevInit16;
	
	const char* fwdPtr = seqStart - 3;
	const char* revPtr = refStart - 3;
	
	int seqScoreIndex = -2;
	while (fwdPtr <= seqEnd - 1)
	{
		__m128i fwdOdd = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i revOdd = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmpOdd8 = _mm_cmpeq_epi8(fwdOdd,revOdd);
		__m128i cmpOdd = _mm_unpacklo_epi8(cmpOdd8,cmpOdd8);
		
		// matchOdd = scoreOdd + cmpOdd ? mMatchBonus16 : -mMismatchPenalty16
		__m128i matchOdd = _mm_sub_epi16(_mm_add_epi16(scoreOdd,_mm_and_si128(cmpOdd,mMatchBonus16)),_mm_andnot_si128(cmpOdd,mMismatchPenalty16));
		
		// gapOdd = scoreEven - mGapPenalty16
		__m128i gapOdd = _mm_sub_epi16(scoreEven,mGapPenalty16);
		
		// max
		__m128i nextScoreOdd = _mm_min_epi16(_mm_max_epi16(matchOdd,_mm_max_epi16(_mm_or_si128(_mm_srli_si128(gapOdd,2),mShiftRight16),gapOdd)),mClampCornerScoreOdd16);
		
		revPtr--;
		
		__m128i fwdEven = _mm_loadu_si128((__m128i*)fwdPtr);
		__m128i revEven = _mm_loadu_si128((__m128i*)revPtr);
		
		__m128i cmpEven8 = _mm_cmpeq_epi8(fwdEven,revEven);
		__m128i cmpEven = _mm_unpacklo_epi8(cmpEven8,cmpEven8);
		
		// matchEven = scoreEven + cmpEven ? mMatchBonus16 : -mMismatchPenalty16
		__m128i matchEven = _mm_sub_epi16(_mm_add_epi16(scoreEven,_mm_and_si128(cmpEven,mMatchBonus16)),_mm_andnot_si128(cmpEven,mMismatchPenalty16));
		
		// gapEven = scoreOdd - mGapPenalty16
		__m128i gapEven = _mm_sub_epi16(nextScoreOdd,mGapPenalty16);
		
		// max
		__m128i nextScoreEven = _mm_max_epi16(matchEven,_mm_max_epi16(_mm_or_si128(_mm_slli_si128(gapEven,2),mShiftLeft16),gapEven));
		
		fwdPtr++;
		
		scoreOdd = nextScoreOdd;
		scoreEven = nextScoreEven;
		
		maxScore = _mm_srli_si128(maxScore,2);
		maxScore = _mm_or_si128(maxScore,mShiftRight16);
		
		refLen = _mm_srli_si128(refLen,2);
		
		__m128i scoreCmpOdd = _mm_cmpgt_epi16(scoreOdd,maxScore);
		refLen = _mm_or_si128(_mm_and_si128(scoreCmpOdd,curLen),_mm_andnot_si128(scoreCmpOdd,refLen));
		
		maxScore = _mm_max_epi16(maxScore,scoreOdd);
		
		curLen = _mm_add_epi16(curLen,mIncrement16);
		
		__m128i scoreCmpEven = _mm_cmpgt_epi16(scoreEven,maxScore);
		refLen = _mm_or_si128(_mm_and_si128(scoreCmpEven,curLen),_mm_andnot_si128(scoreCmpEven,refLen));
		
		maxScore = _mm_max_epi16(maxScore,scoreEven);
		
		seqScores[max(0,seqScoreIndex)] = _mm_extract_epi16(maxScore,0);
		refLengths[max(0,seqScoreIndex)] = _mm_extract_epi16(refLen,0);
		seqScoreIndex++;
	}
}

void SimpleAligner::AlignBanded(const string& reference, const string& sequence, int bandWidth, int& score, int& refLength, int& seqLength)
{
	int matrixLength = reference.size() + 1;
	int matrixHeight = sequence.size() + 1;
	
	if (matrixLength != mCurrent.size())
	{
		mCurrent.resize(matrixLength);
		mPrevious.resize(matrixLength);
	}
	
	int badScore = numeric_limits<int>::min() / 2;
	
	score = 0;
	refLength = 0;
	seqLength = 0;
	
	for (int j = 0; j < matrixHeight; j++)
	{
		int seqPos = j - 1;
		
		int istart = j - bandWidth;
		int iend = j + bandWidth;
		
		if (istart > 0)
		{
			mCurrent[istart-1] = badScore;
		}
		
		for (int i = max(0,istart); i <= min(matrixLength-1,iend); i++)
		{
			int refPos = i - 1;
			
			if (i == 0 && j == 0)
			{
				mCurrent[i] = 0;
			}
			else if (j == 0)
			{
				mCurrent[i] = mCurrent[i-1] + mGapScore;
			}
			else if (i == 0)
			{
				mCurrent[i] = mPrevious[i] + mGapScore;
			}
			else
			{
				int matchScore = mPrevious[i-1] + ((reference[refPos] == sequence[seqPos]) ? mMatchScore : mMisMatchScore);
				int gapRefScore = mCurrent[i-1] + mGapScore;
				int gapReadScore = mPrevious[i] + mGapScore;
				int maxScore = max(matchScore,max(gapRefScore,gapReadScore));
				
				if (maxScore > score)
				{
					score = maxScore;
					refLength = refPos + 1;
					seqLength = seqPos + 1;
				}
				
				mCurrent[i] = maxScore;
			}
		}
		
		if (iend + 1 < matrixLength)
		{
			mCurrent[iend + 1] = badScore;
		}
		
		swap(mCurrent, mPrevious);
	}
}

void SimpleAligner::AlignBanded(const string& reference, const string& sequence, int bandWidth, IntegerVec& seqLengthScores, IntegerVec& refLengths)
{
	int matrixLength = reference.size() + 1;
	int matrixHeight = sequence.size() + 1;
	
	if (matrixLength != mCurrent.size())
	{
		mCurrent.resize(matrixLength + 2 * bandWidth + 1);
		mPrevious.resize(matrixLength + 2 * bandWidth + 1);
	}
	
	seqLengthScores.resize(sequence.size() + 1);
	refLengths.resize(sequence.size() + 1);
	
	int badScore = numeric_limits<int>::min() / 2;
	
	for (int j = 0; j < matrixHeight; j++)
	{
		int seqPos = j - 1;
		
		int istart = j - bandWidth;
		int iend = j + bandWidth;
		
		if (istart > 0)
		{
			mCurrent[istart-1] = badScore;
		}
		
		int seqLengthScore = badScore;
		int refLength;
		
		for (int i = max(0,istart); i <= iend; i++)
		{
			int refPos = i - 1;
			
			if (i == 0 && j == 0)
			{
				mCurrent[i] = 0;
			}
			else if (j == 0)
			{
				mCurrent[i] = mCurrent[i-1] + mGapScore;
			}
			else if (i == 0)
			{
				mCurrent[i] = mPrevious[i] + mGapScore;
			}
			else
			{
				int matchScore = mPrevious[i-1] + mMisMatchScore;
				if (refPos < reference.size())
				{
					matchScore = mPrevious[i-1] + ((reference[refPos] == sequence[seqPos]) ? mMatchScore : mMisMatchScore);
				}

				int gapRefScore = mCurrent[i-1] + mGapScore;
				int gapReadScore = mPrevious[i] + mGapScore;
				int maxScore = max(matchScore,max(gapRefScore,gapReadScore));
				
				mCurrent[i] = maxScore;
			}
			
			if (mCurrent[i] > seqLengthScore)
			{
				seqLengthScore = mCurrent[i];
				refLength = i;
			}
		}
		
		seqLengthScores[j] = seqLengthScore;
		refLengths[j] = refLength;
		
		mCurrent[iend + 1] = badScore;
		
		swap(mCurrent, mPrevious);
	}
}

void SimpleAligner::AlignBandedFwd(const char* refStart, const char* seqStart, const char* seqEnd, int bandWidth, int& score, int& refLength, int& seqLength)
{
	int matrixLength = seqEnd - seqStart + 1 + bandWidth;
	int matrixHeight = seqEnd - seqStart + 1;
	
	if (matrixLength != mCurrent.size())
	{
		mCurrent.resize(matrixLength);
		mPrevious.resize(matrixLength);
	}
	
	int badScore = numeric_limits<int>::min() / 2;
	
	score = 0;
	refLength = 0;
	seqLength = 0;
	
	for (int j = 0; j < matrixHeight; j++)
	{
		int seqPos = j - 1;
		
		int istart = j - bandWidth;
		int iend = j + bandWidth;
		
		if (istart > 0)
		{
			mCurrent[istart-1] = badScore;
		}
		
		for (int i = max(0,istart); i <= iend; i++) 
		{
			DebugCheck(i < mCurrent.size());
			
			int refPos = i - 1;
			
			if (i == 0 && j == 0)
			{
				mCurrent[i] = 0;
			}
			else if (j == 0)
			{
				mCurrent[i] = mCurrent[i-1] + mGapScore;
			}
			else if (i == 0)
			{
				mCurrent[i] = mPrevious[i] + mGapScore;
			}
			else
			{
				int matchScore = mPrevious[i-1] + ((refStart[refPos] == seqEnd[-seqPos - 1]) ? mMatchScore : mMisMatchScore);
				int gapRefScore = mCurrent[i-1] + mGapScore;
				int gapReadScore = mPrevious[i] + mGapScore;
				int maxScore = max(matchScore,max(gapRefScore,gapReadScore));
				
				if (maxScore > score)
				{
					score = maxScore;
					refLength = refPos + 1;
					seqLength = seqPos + 1;
				}
				
				mCurrent[i] = maxScore;
			}
		}
		
		if (iend + 1 < matrixLength)
		{
			mCurrent[iend + 1] = badScore;
		}
		
		swap(mCurrent, mPrevious);
	}
}

void SimpleAligner::AlignBandedRev(const char* refStart, const char* seqStart, const char* seqEnd, int bandWidth, int& score, int& refLength, int& seqLength)
{
	int matrixLength = seqEnd - seqStart + 1 + bandWidth;
	int matrixHeight = seqEnd - seqStart + 1;
	
	if (matrixLength != mCurrent.size())
	{
		mCurrent.resize(matrixLength);
		mPrevious.resize(matrixLength);
	}
	
	int badScore = numeric_limits<int>::min() / 2;
	
	score = 0;
	refLength = 0;
	seqLength = 0;
	
	for (int j = 0; j < matrixHeight; j++)
	{
		int seqPos = j - 1;
		
		int istart = j - bandWidth;
		int iend = j + bandWidth;
		
		if (istart > 0)
		{
			mCurrent[istart-1] = badScore;
		}
		
		for (int i = max(0,istart); i <= iend; i++) 
		{
			DebugCheck(i < mCurrent.size());
			
			int refPos = i - 1;
			
			if (i == 0 && j == 0)
			{
				mCurrent[i] = 0;
			}
			else if (j == 0)
			{
				mCurrent[i] = mCurrent[i-1] + mGapScore;
			}
			else if (i == 0)
			{
				mCurrent[i] = mPrevious[i] + mGapScore;
			}
			else
			{
				int matchScore = mPrevious[i-1] + ((refStart[-refPos] == seqStart[seqPos]) ? mMatchScore : mMisMatchScore);
				int gapRefScore = mCurrent[i-1] + mGapScore;
				int gapReadScore = mPrevious[i] + mGapScore;
				int maxScore = max(matchScore,max(gapRefScore,gapReadScore));
				
				if (maxScore > score)
				{
					score = maxScore;
					refLength = refPos + 1;
					seqLength = seqPos + 1;
				}
				
				mCurrent[i] = maxScore;
			}
		}
		
		if (iend + 1 < matrixLength)
		{
			mCurrent[iend + 1] = badScore;
		}
		
		swap(mCurrent, mPrevious);
	}
}

void SimpleAligner::AlignPrint(const string& seq1, const string& seq2)
{
	int matrixLength = seq1.size() + 1;
	int matrixHeight = seq2.size() + 1;
	
	Matrix<int> matrix;
	Matrix<Cell> backTrace;
	
	matrix.Resize(matrixLength, matrixHeight);
	backTrace.Resize(matrixLength, matrixHeight);
	
	int bestScore = 0;
	int previousBestScore = 0;
	
	int increaseCount = 0;
	int plateauCount = 0;
	
	int bestChangepointScore = 0;
	int bestChangepointLength = 0;
	Cell bestChangepointCell;
	
	for (int i = 0; i < matrixLength; i++) 
	{
		Cell bestCell;
		
		int seq1Pos = i - 1;
		
		for (int j = 0; j < matrixHeight; j++)
		{
			int seq2Pos = j - 1;
			
			if (i == 0 && j == 0)
			{
				matrix(i,j) = 0;
			}
			else if (j == 0)
			{
				matrix(i,j) = matrix(i-1,j) + mGapScore;
				backTrace(i,j) = Cell(i,j-1);
			}
			else if (i == 0)
			{
				matrix(i,j) = matrix(i,j-1) + mGapScore;
				backTrace(i,j) = Cell(i,j-1);
			}
			else
			{
				int matchScore = matrix(i-1,j-1) + ((seq1[seq1Pos] == seq2[seq2Pos]) ? mMatchScore : mMisMatchScore);
				int gapRefScore = matrix(i-1,j) + mGapScore;
				int gapReadScore = matrix(i,j-1) + mGapScore;
				int maxScore = max(matchScore,max(gapRefScore,gapReadScore));
				
				if (matchScore == maxScore)
				{
					backTrace(i,j) = Cell(i-1,j-1);
				}
				
				if (gapRefScore == maxScore)
				{
					backTrace(i,j) = Cell(i-1,j);
				}
				
				if (gapReadScore == maxScore)
				{
					backTrace(i,j) = Cell(i,j-1);
				}
				
				if (maxScore > bestScore)
				{
					bestScore = maxScore;
					bestCell = Cell(i,j);
				}
				
				matrix(i,j) = maxScore;
			}
		}
		
		if (bestScore > previousBestScore)
		{
			increaseCount++;
		}
		else
		{
			plateauCount++;
		}
		
		int changepointScore = increaseCount - plateauCount;
		
		if (changepointScore > bestChangepointScore)
		{
			bestChangepointScore = changepointScore;
			bestChangepointLength = seq1Pos + 1;
			bestChangepointCell = bestCell;
		}
		
		previousBestScore = bestScore;
	}
	
	IntegerPairVec matches;
	
	Cell cell = bestChangepointCell;
	while (cell.j > 0)
	{
		Cell nextCell = backTrace(cell);
		
		int seq1Pos = cell.i - 1;
		int seq2Pos = cell.j - 1;
		
		if (cell.i - 1 == nextCell.i && cell.j - 1 == nextCell.j)
		{
			matches.push_back(IntegerPair(seq1Pos,seq2Pos));
		}
		
		cell = nextCell;
	}
	
	reverse(matches.begin(), matches.end());
	
	string seq1Gapped;
	string seq2Gapped;
	string seqMatches;
	IntegerPair last(-1,-1);
	for (IntegerPairVecConstIter matchIter = matches.begin(); matchIter != matches.end(); matchIter++)
	{
		for (int seq1Pos = last.first + 1; seq1Pos < matchIter->first; seq1Pos++)
		{
			seq1Gapped.push_back(seq1[seq1Pos]);
			seq2Gapped.append("-");
			seqMatches.append(" ");
		}
		
		for (int seq2Pos = last.second + 1; seq2Pos < matchIter->second; seq2Pos++)
		{
			seq1Gapped.append("-");
			seq2Gapped.push_back(seq2[seq2Pos]);
			seqMatches.append(" ");
		}
		
		seq1Gapped.push_back(seq1[matchIter->first]);
		seq2Gapped.push_back(seq2[matchIter->second]);
		seqMatches.append(seq1[matchIter->first] == seq2[matchIter->second] ? "|" : " ");
		
		last = *matchIter;
	}
	
	cout << seq1Gapped << endl;
	cout << seqMatches << endl;	
	cout << seq2Gapped << endl;
}



