/*
 *  testssealign.cpp
 */

#include "Common.h"
#include "DebugCheck.h"
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


void GetRandomSequence(int length, string& sequence)
{
	char ntchars[] = {'A','C','T','G'};
	
	while (sequence.length() < length)
	{
		sequence = sequence + ntchars[rand() % 4];
	}
}

#define CHECKTEST(expr) if(!(expr)) { cout << #expr << endl; cout << seq1 << endl << seq2 << endl; exit(1); }

void RandomSequenceTestsEndToEnd()
{
	cout << "End to end alignment test" << endl;
	
	SimpleAligner aligner1(2, -3, -4);
	
	const int cPadding = 16;
	
	srand(1);
	for (int i = 0; i < 100000; i++)
	{
		int length1 = 200;
		
		string seq1;
		GetRandomSequence(length1, seq1);
		
		int length2 = rand() % 50 + 30;
		
		string seq2;
		GetRandomSequence(length2, seq2);
		
		string reference = seq1;
		string sequence = seq2;
		
		int score = aligner1.AlignEndToEnd(seq1, seq2);
		
		int scoreFwd;
		IntegerVec refLengthScoresFwd;
		{
			string refPad = string(cPadding,'X') + reference + string(cPadding,'X');
			string seqPad = string(cPadding,'Y') + sequence + string(cPadding,'Y');
			reverse(seqPad.begin(), seqPad.end());
		
			const char* refStart = refPad.c_str() + cPadding;
			const char* refEnd = refPad.c_str() + reference.length() + cPadding;
			const char* seqStart = seqPad.c_str() + cPadding;
			const char* seqEnd = seqPad.c_str() + sequence.length() + cPadding;
			
			aligner1.AlignEndToEndFwd(refStart, refEnd, seqStart, seqEnd, refLengthScoresFwd);
			scoreFwd = *max_element(refLengthScoresFwd.begin(), refLengthScoresFwd.end());
		}
		
		int scoreRev;
		IntegerVec refLengthScoresRev;
		{
			string refPad = string(cPadding,'X') + reference + string(cPadding,'X');
			string seqPad = string(cPadding,'Y') + sequence + string(cPadding,'Y');
			reverse(refPad.begin(), refPad.end());
			
			const char* refStart = refPad.c_str() + cPadding;
			const char* refEnd = refPad.c_str() + reference.length() + cPadding;
			const char* seqStart = seqPad.c_str() + cPadding;
			const char* seqEnd = seqPad.c_str() + sequence.length() + cPadding;
			
			aligner1.AlignEndToEndRev(refStart, refEnd, seqStart, seqEnd, refLengthScoresRev);
			scoreRev = *max_element(refLengthScoresRev.begin(), refLengthScoresRev.end());
		}
		
		CHECKTEST(score == scoreFwd)
		CHECKTEST(score == scoreRev)
	}
	
	cout << "Passed" << endl;
}

void RandomSequenceTestsSeed()
{
	cout << "Seed alignment test" << endl;
	
	SimpleAligner aligner1(2, -3, -4);
	
	const int cPadding = 16;
	
	srand(1);
	for (int i = 0; i < 100000; i++)
	{
		int length1 = 200;
		
		string seq1;
		GetRandomSequence(length1, seq1);
		
		int length2 = 16;
		
		string seq2;
		GetRandomSequence(length2, seq2);
		
		string reference = seq1;
		string sequence = seq2;
		
		int score = aligner1.AlignEndToEnd(seq1, seq2);
		
		int scoreFwd;
		IntegerVec refLengthScoresFwd;
		{
			string refPad = string(cPadding,'X') + reference + string(cPadding,'X');
			string seqPad = string(cPadding,'Y') + sequence + string(cPadding,'Y');
			reverse(seqPad.begin(), seqPad.end());
			
			const char* refStart = refPad.c_str() + cPadding;
			const char* refEnd = refPad.c_str() + reference.length() + cPadding;
			const char* seqStart = seqPad.c_str() + cPadding;
			const char* seqEnd = seqPad.c_str() + sequence.length() + cPadding;
			
			aligner1.AlignEndToEndFwd(refStart, refEnd, seqStart, seqEnd, refLengthScoresFwd);
			scoreFwd = *max_element(refLengthScoresFwd.begin(), refLengthScoresFwd.end());
		}
		
		int scoreRev;
		IntegerVec refLengthScoresRev;
		{
			string refPad = string(cPadding,'X') + reference + string(cPadding,'X');
			string seqPad = string(cPadding,'Y') + sequence + string(cPadding,'Y');
			reverse(refPad.begin(), refPad.end());
			
			const char* refStart = refPad.c_str() + cPadding;
			const char* refEnd = refPad.c_str() + reference.length() + cPadding;
			const char* seqStart = seqPad.c_str() + cPadding;
			const char* seqEnd = seqPad.c_str() + sequence.length() + cPadding;
			
			aligner1.AlignEndToEndRev(refStart, refEnd, seqStart, seqEnd, refLengthScoresRev);
			scoreRev = *max_element(refLengthScoresRev.begin(), refLengthScoresRev.end());
		}
		
		IntegerVec scores16baseSSE2Fwd;
		{
			string refPad = string(cPadding,'X') + reference + string(cPadding,'X');
			string seqPad = string(cPadding,'Y') + sequence + string(cPadding,'Y');
			reverse(seqPad.begin(), seqPad.end());
			
			const char* refStart = refPad.c_str() + cPadding;
			const char* refEnd = refPad.c_str() + reference.length() + cPadding;
			const char* seqStart = seqPad.c_str() + cPadding;
			const char* seqEnd = seqPad.c_str() + sequence.length() + cPadding;
			
			vector<short int> scoresBuffer(1000);
			short int* scoresPtr = &scoresBuffer[10];
			aligner1.Align16baseSSE2Fwd(refStart, refEnd, seqStart, scoresPtr);
			scores16baseSSE2Fwd = IntegerVec(scoresPtr, scoresPtr + seq1.size() + 1);
		}
		
		IntegerVec scores16baseSSE2Rev;
		{
			
			string refPad = string(cPadding,'X') + reference + string(cPadding,'X');
			string seqPad = string(cPadding,'Y') + sequence + string(cPadding,'Y');
			reverse(refPad.begin(), refPad.end());
			
			const char* refStart = refPad.c_str() + cPadding;
			const char* refEnd = refPad.c_str() + reference.length() + cPadding;
			const char* seqStart = seqPad.c_str() + cPadding;
			const char* seqEnd = seqPad.c_str() + sequence.length() + cPadding;
			
			vector<short int> scoresBuffer(1000);
			short int* scoresPtr = &scoresBuffer[10];
			aligner1.Align16baseSSE2Rev(refStart, refEnd, seqStart, scoresPtr);
			scores16baseSSE2Rev = IntegerVec(scoresPtr, scoresPtr + seq1.size() + 1);
		}
		
		CHECKTEST(score == scoreFwd)
		CHECKTEST(score == scoreRev)
		CHECKTEST(refLengthScoresFwd.size() == scores16baseSSE2Fwd.size() && equal(refLengthScoresFwd.begin(), refLengthScoresFwd.end(), scores16baseSSE2Fwd.begin()))
		CHECKTEST(refLengthScoresRev.size() == scores16baseSSE2Rev.size() && equal(refLengthScoresRev.begin(), refLengthScoresRev.end(), scores16baseSSE2Rev.begin()))
	}
	
	cout << "Passed" << endl;
}

void RandomSequenceTestsBanded()
{
	cout << "Banded alignment test" << endl;
	
	SimpleAligner aligner1(2, -3, -4);
	
	srand(1);
	for (int i = 0; i < 100000; i++)
	{
		int length1 = 200;
		
		string seq1;
		GetRandomSequence(length1, seq1);
		
		int length2 = rand() % 50 + 30;
		
		string seq2;
		GetRandomSequence(length2, seq2);
		
		string reference = seq1;
		string sequence = seq2;
		
		IntegerVec seqLengthScores;
		IntegerVec refLengths;
		aligner1.AlignBanded(reference, sequence, 7, seqLengthScores, refLengths);
		
		int score;
		int refLength;
		int seqLength;
		aligner1.AlignBanded(seq1, seq2, 7, score, refLength, seqLength);
		
		int scoreFwdSSE;
		{
			string refPad = string(8,'X') + reference + string(8+max(0,(int)(sequence.length()-reference.length())),'X');
			string seqPad = string(8,'Y') + sequence + string(8+max(0,(int)(reference.length()-sequence.length())),'Y');
			reverse(seqPad.begin(), seqPad.end());
			
			const char* refStart = refPad.c_str() + 8;
			const char* seqStart = seqPad.c_str() + seqPad.length() - sequence.length() - 8;
			const char* seqEnd = seqPad.c_str() + seqPad.length() - 8;
			
			scoreFwdSSE = aligner1.AlignBandedSSE2BW7ScoreFwd(refStart, seqStart, seqEnd);
		}
		
		int scoreRevSSE;
		{
			string refPad = string(8,'X') + reference + string(8+max(0,(int)(sequence.length()-reference.length())),'X');
			string seqPad = string(8,'Y') + sequence + string(8+max(0,(int)(reference.length()-sequence.length())),'Y');
			reverse(refPad.begin(), refPad.end());
			
			const char* refStart = refPad.c_str() + refPad.length() - 9;
			const char* seqStart = seqPad.c_str() + 8;
			const char* seqEnd = seqPad.c_str() + sequence.length() + 8;
			
			scoreRevSSE = aligner1.AlignBandedSSE2BW7ScoreRev(refStart, seqStart, seqEnd);
		}
		
		int scoreFwd;
		int refLengthFwd;
		int seqLengthFwd;
		{
			string refPad = string(8,'X') + reference + string(8+max(0,(int)(sequence.length()-reference.length())),'X');
			string seqPad = string(8,'Y') + sequence + string(8+max(0,(int)(reference.length()-sequence.length())),'Y');
			reverse(seqPad.begin(), seqPad.end());
			
			const char* refStart = refPad.c_str() + 8;
			const char* seqStart = seqPad.c_str() + seqPad.length() - sequence.length() - 8;
			const char* seqEnd = seqPad.c_str() + seqPad.length() - 8;
			
			aligner1.AlignBandedFwd(refStart, seqStart, seqEnd, 7, scoreFwd, refLengthFwd, seqLengthFwd);
		}
		
		int scoreRev;
		int refLengthRev;
		int seqLengthRev;
		{
			string refPad = string(8,'X') + reference + string(8+max(0,(int)(sequence.length()-reference.length())),'X');
			string seqPad = string(8,'Y') + sequence + string(8+max(0,(int)(reference.length()-sequence.length())),'Y');
			reverse(refPad.begin(), refPad.end());
			
			const char* refStart = refPad.c_str() + refPad.length() - 9;
			const char* seqStart = seqPad.c_str() + 8;
			const char* seqEnd = seqPad.c_str() + sequence.length() + 8;
			
			aligner1.AlignBandedRev(refStart, seqStart, seqEnd, 7, scoreRev, refLengthRev, seqLengthRev);
		}
		
		IntegerVec scoresSSEFwd1;
		{
			string refPad = string(8,'X') + reference + string(8+max(0,(int)(sequence.length()-reference.length())),'X');
			string seqPad = string(8,'Y') + sequence + string(8+max(0,(int)(reference.length()-sequence.length())),'Y');
			reverse(seqPad.begin(), seqPad.end());
			
			const char* refStart = refPad.c_str() + 8;
			const char* seqStart = seqPad.c_str() + seqPad.length() - sequence.length() - 8;
			const char* seqEnd = seqPad.c_str() + seqPad.length() - 8;
			
			vector<short int> scoresBuffer(100);
			
			short int* scoresPtr = &scoresBuffer[10];
			
			aligner1.AlignBandedSSE2BW7ScoreFwd(refStart, seqStart, seqEnd, scoresPtr);
			
			scoresSSEFwd1 = IntegerVec(scoresPtr, scoresPtr + seq2.size() + 1);
		}
		
		IntegerVec scoresSSEFwd2;
		IntegerVec lengthsSSEFwd;
		{
			string refPad = string(8,'X') + reference + string(8+max(0,(int)(sequence.length()-reference.length())),'X');
			string seqPad = string(8,'Y') + sequence + string(8+max(0,(int)(reference.length()-sequence.length())),'Y');
			reverse(seqPad.begin(), seqPad.end());
			
			const char* refStart = refPad.c_str() + 8;
			const char* seqStart = seqPad.c_str() + seqPad.length() - sequence.length() - 8;
			const char* seqEnd = seqPad.c_str() + seqPad.length() - 8;
			
			vector<short int> scoresBuffer(100);
			vector<short int> lengthsBuffer(100);
			
			short int* scoresPtr = &scoresBuffer[10];
			short int* lengthsPtr = &lengthsBuffer[10];
			
			aligner1.AlignBandedSSE2BW7ScoreFwd(refStart, seqStart, seqEnd, scoresPtr, lengthsPtr);
			
			scoresSSEFwd2 = IntegerVec(scoresPtr, scoresPtr + seq2.size() + 1);
			lengthsSSEFwd = IntegerVec(lengthsPtr, lengthsPtr + seq2.size() + 1);
		}
		
		IntegerVec scoresSSERev1;
		{
			string refPad = string(8,'X') + reference + string(8+max(0,(int)(sequence.length()-reference.length())),'X');
			string seqPad = string(8,'Y') + sequence + string(8+max(0,(int)(reference.length()-sequence.length())),'Y');
			reverse(refPad.begin(), refPad.end());
			
			const char* refStart = refPad.c_str() + refPad.length() - 9;
			const char* seqStart = seqPad.c_str() + 8;
			const char* seqEnd = seqPad.c_str() + sequence.length() + 8;
			
			vector<short int> scoresBuffer(100);
			
			short int* scoresPtr = &scoresBuffer[10];
			
			aligner1.AlignBandedSSE2BW7ScoreRev(refStart, seqStart, seqEnd, scoresPtr);
			
			scoresSSERev1 = IntegerVec(scoresPtr, scoresPtr + seq2.size() + 1);
		}
		
		IntegerVec scoresSSERev2;
		IntegerVec lengthsSSERev;
		{
			string refPad = string(8,'X') + reference + string(8+max(0,(int)(sequence.length()-reference.length())),'X');
			string seqPad = string(8,'Y') + sequence + string(8+max(0,(int)(reference.length()-sequence.length())),'Y');
			reverse(refPad.begin(), refPad.end());
			
			const char* refStart = refPad.c_str() + refPad.length() - 9;
			const char* seqStart = seqPad.c_str() + 8;
			const char* seqEnd = seqPad.c_str() + sequence.length() + 8;
			
			vector<short int> scoresBuffer(100);
			vector<short int> lengthsBuffer(100);
			
			short int* scoresPtr = &scoresBuffer[10];
			short int* lengthsPtr = &lengthsBuffer[10];
			
			aligner1.AlignBandedSSE2BW7ScoreRev(refStart, seqStart, seqEnd, scoresPtr, lengthsPtr);
			
			scoresSSERev2 = IntegerVec(scoresPtr, scoresPtr + seq2.size() + 1);
			lengthsSSERev = IntegerVec(lengthsPtr, lengthsPtr + seq2.size() + 1);
		}
		
		CHECKTEST(scoreFwd == score)
		CHECKTEST(scoreRev == score)
		CHECKTEST(scoreFwdSSE == score)
		CHECKTEST(scoreRevSSE == score)
		CHECKTEST(refLength == refLengthFwd)
		CHECKTEST(refLength == refLengthRev)
		CHECKTEST(seqLength == seqLengthFwd)
		CHECKTEST(seqLength == seqLengthRev)
		CHECKTEST(score == *max_element(seqLengthScores.begin(), seqLengthScores.end()))
		CHECKTEST(seqLengthScores.size() == scoresSSEFwd1.size() && equal(seqLengthScores.begin(), seqLengthScores.end(), scoresSSEFwd1.begin()))
		CHECKTEST(seqLengthScores.size() == scoresSSEFwd2.size() && equal(seqLengthScores.begin(), seqLengthScores.end(), scoresSSEFwd2.begin()))
		CHECKTEST(refLengths.size() == lengthsSSEFwd.size() && equal(refLengths.begin(), refLengths.end(), lengthsSSEFwd.begin()))
		CHECKTEST(seqLengthScores.size() == scoresSSERev1.size() && equal(seqLengthScores.begin(), seqLengthScores.end(), scoresSSERev1.begin()))
		CHECKTEST(seqLengthScores.size() == scoresSSERev2.size() && equal(seqLengthScores.begin(), seqLengthScores.end(), scoresSSERev2.begin()))
		CHECKTEST(refLengths.size() == lengthsSSERev.size() && equal(refLengths.begin(), refLengths.end(), lengthsSSERev.begin()))
	}
	
	cout << "Passed" << endl;
}

int main()
{
	RandomSequenceTestsSeed();
	RandomSequenceTestsEndToEnd();
	RandomSequenceTestsBanded();
}

