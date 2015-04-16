/*
 *  SimpleAligner.h
 *
 *  Created by Andrew McPherson on 05/09/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SIMPLEALIGNER_H_
#define SIMPLEALIGNER_H_

#include "Common.h"
#include "Matrix.h"

#include <string>
#include <vector>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>

using namespace std;

class SimpleAligner
{
public:
	SimpleAligner(int matchScore, int misMatchScore, int gapScore);
	
	int AlignEndToEnd(const string& reference, const string& sequence);
	void Align16baseSSE2Fwd(const char* refStart, const char* refEnd, const char* seqEnd, short int* refScores);
	void Align16baseSSE2Rev(const char* refStart, const char* refEnd, const char* seqStart, short int* refScores);
	void AlignEndToEndFwd(const char* refStart, const char* refEnd, const char* seqStart, const char* seqEnd, IntegerVec& refLengthScores);
	void AlignEndToEndRev(const char* refStart, const char* refEnd, const char* seqStart, const char* seqEnd, IntegerVec& refLengthScores);
	void AlignPartial(const string& reference, const string& sequence, IntegerVec& seqLengthScores);
	void AlignBanded(const string& reference, const string& sequence, int bandWidth, int& score, int& refLength, int& seqLength);
	void AlignBanded(const string& reference, const string& sequence, int bandWidth, IntegerVec& seqLengthScores, IntegerVec& refLengths);
	void AlignBandedFwd(const char* refStart, const char* seqStart, const char* seqEnd, int bandWidth, int& score, int& refLength, int& seqLength);
	void AlignBandedRev(const char* refStart, const char* seqStart, const char* seqEnd, int bandWidth, int& score, int& refLength, int& seqLength);
	int AlignBandedSSE2BW7ScoreFwd(const char* refStart, const char* seqStart, const char* seqEnd);
	void AlignBandedSSE2BW7ScoreFwd(const char* refStart, const char* seqStart, const char* seqEnd, short int* seqScores);
	void AlignBandedSSE2BW7ScoreFwd(const char* refStart, const char* seqStart, const char* seqEnd, short int* seqScores, short int* refLengths);
	int AlignBandedSSE2BW7ScoreRev(const char* refStart, const char* seqStart, const char* seqEnd);
	void AlignBandedSSE2BW7ScoreRev(const char* refStart, const char* seqStart, const char* seqEnd, short int* seqScores);
	void AlignBandedSSE2BW7ScoreRev(const char* refStart, const char* seqStart, const char* seqEnd, short int* seqScores, short int* refLengths);
	void AlignPrint(const string& seq1, const string& seq2);

private:
	int mMatchScore;
	int mMisMatchScore;
	int mGapScore;
	int mBandWidth;
	
	vector<int> mCurrent;
	vector<int> mPrevious;
	
	vector<int> mD0;
	vector<int> mD1;
	vector<int> mP0;
	vector<int> mP1;
	vector<int> mQ0;
	vector<int> mQ1;
	
	vector<unsigned char> mScoreBuffer8;
	vector<short int> mScoreBuffer16;
	
	__m128i mIncrement16;
	__m128i mMatchBonus16;
	__m128i mMismatchPenalty16;
	__m128i mGapPenalty16;
	__m128i mMaxScoreInit16;
	__m128i mShiftRight16;
	__m128i mShiftLeft16;
	__m128i mScoreOddInit16;
	__m128i mScoreEvenInit16;
	__m128i mClampCornerScoreOdd16;
	
	__m128i mCurLenFwdInit16;
	__m128i mRefLenFwdInit16;
	__m128i mMaxScoreFwdInit16;
	__m128i mCurLenRevInit16;
	__m128i mRefLenRevInit16;
	__m128i mMaxScoreRevInit16;

	uint8_t mOffset8;
	uint8_t mZero8;
	__m128i mMatchBonus8;
	__m128i mMismatchPenalty8;
	__m128i mGapPenalty8;
	
	__m128i mShiftRefGapFwd8;
	__m128i mScoreDiagonalInitFwd8;
	__m128i mScoreOrthogonalInitFwd8;
	
	__m128i mShiftRefGapRev8;
	__m128i mScoreDiagonalInitRev8;
	__m128i mScoreOrthogonalInitRev8;
};

#endif


