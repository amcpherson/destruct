/*
 *  AffineGapDP.cpp
 *
 *  Created by Andrew McPherson on 18/07/12.
 *
 */

#include "AffineGapDP.h"
#include "DebugCheck.h"

#include <iostream>
#include <list>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>

using namespace std;


AffineGapDP::AffineGapDP()
{
}

/*
 --ma <int>
 Sets the match bonus. In --local mode <int> is added to the alignment score for each position where a read character aligns
 to a reference character and the characters match. Not used in --end-to-end mode. Default: 2.
 
 --mp MX,MN
 Sets the maximum (MX) and minimum (MN) mismatch penalties, both integers. A number less than or equal to MX and greater than
 or equal to MN is subtracted from the alignment score for each position where a read character aligns to a reference character,
 the characters do not match, and neither is an N. If --ignore-quals is specified, the number subtracted quals MX. Otherwise, 
 the number subtracted is MN + floor( (MX-MN)(MIN(Q, 40.0)/40.0) ) where Q is the Phred quality value. Default: MX = 6, MN = 2.
 
 --np <int>
 Sets penalty for positions where the read, reference, or both, contain an ambiguous character such as N. Default: 1.
 
 --rdg <int1>,<int2>
 Sets the read gap open (<int1>) and extend (<int2>) penalties. A read gap of length N gets a penalty of <int1> + N * <int2>. 
 Default: 5, 3.
 
 --rfg <int1>,<int2>
 Sets the reference gap open (<int1>) and extend (<int2>) penalties. A reference gap of length N gets a penalty of <int1> + N * <int2>. 
 Default: 5, 3.
 
 --score-min <func>
 Sets a function governing the minimum alignment score needed for an alignment to be considered "valid" (i.e. good enough to report). 
 This is a function of read length. For instance, specifying L,0,-0.6 sets the minimum-score function f to f(x) = 0 + -0.6 * x, where
 x is the read length. See also: setting function options. The default in --end-to-end mode is L,-0.6,-0.6 and the default in --local mode is G,20,8.
*/

int AffineGapDP::EndToEndAlign(const string& reference, const string& sequence, const string& quality)
{
	const int ma = 0;
	const int mx = 6;
	const int mn = 2;
	const int rdgo = 5;
	const int rdge = 3;
	const int rfgo = 5;
	const int rfge = 3;
	
	int minscore = -(mx + rdgo + rdge + rfgo + rfge) * sequence.length();
	
	int matrixLength = reference.size() + 1;
	int matrixHeight = sequence.size() + 1;
	
	if (matrixLength != mD0.size())
	{
		mD0.resize(matrixLength);
		mD1.resize(matrixLength);
		mP0.resize(matrixLength);
		mP1.resize(matrixLength);
		mQ0.resize(matrixLength);
		mQ1.resize(matrixLength);
	}
	
	int score = numeric_limits<int>::min();
	
	for (int j = 0; j < matrixHeight; j++)
	{
		int seqPos = j - 1;
		
		for (int i = 0; i < matrixLength; i++) 
		{
			int refPos = i - 1;
			
			if (j == 0)
			{
				mP1[i] = minscore;
				mQ1[i] = minscore;
				mD1[i] = 0;
			}
			else if (i == 0)
			{
				mP1[i] = minscore;
				mQ1[i] = minscore;
				mD1[i] = -rdgo - j * rdge;
			}
			else
			{
				double q = (double)(quality[seqPos] - 33);
				int mmp = (int)(mn + floor((mx-mn)*(min(q,40.0)/40.0)));
				int ms = (reference[refPos] == sequence[seqPos]) ? ma : -mmp;
				
				mP1[i] = max(mD1[i-1] - rfgo - rfge, mP1[i-1] - rfge);
				mQ1[i] = max(mD0[i] - rdgo - rdge, mQ0[i] - rdge);
				mD1[i] = max(mD0[i-1] + ms, max(mP1[i], mQ1[i]));
				
				if (j + 1 == matrixHeight && mD1[i] > score)
				{
					score = mD1[i];
				}
			}
		}
		
		swap(mP1, mP0);
		swap(mQ1, mQ0);
		swap(mD1, mD0);
	}
	
	return score;
}

void AffineGapDP::PartialAlign(const string& reference, const string& sequence, const string& quality, IntegerVec& scores)
{
	const int ma = 0;
	const int mx = 6;
	const int mn = 2;
	const int rdgo = 5;
	const int rdge = 3;
	const int rfgo = 5;
	const int rfge = 3;
	
	int minscore = -(mx + rdgo + rdge + rfgo + rfge) * sequence.length();
	
	int matrixLength = reference.size() + 1;
	int matrixHeight = sequence.size() + 1;
	
	if (matrixLength != mD0.size())
	{
		mD0.resize(matrixLength);
		mD1.resize(matrixLength);
		mP0.resize(matrixLength);
		mP1.resize(matrixLength);
		mQ0.resize(matrixLength);
		mQ1.resize(matrixLength);
	}
	
	scores.resize(sequence.length());
	
	for (int j = 0; j < matrixHeight; j++)
	{
		int seqPos = j - 1;
		
		int score = numeric_limits<int>::min();
		
		for (int i = 0; i < matrixLength; i++) 
		{
			int refPos = i - 1;
			
			if (j == 0)
			{
				mP1[i] = minscore;
				mQ1[i] = minscore;
				mD1[i] = 0;
			}
			else if (i == 0)
			{
				mP1[i] = minscore;
				mQ1[i] = minscore;
				mD1[i] = -rdgo - j * rdge;
			}
			else
			{
				double q = (double)(quality[seqPos] - 33);
				int mmp = (int)(mn + floor((mx-mn)*(min(q,40.0)/40.0)));
				int ms = (reference[refPos] == sequence[seqPos]) ? ma : -mmp;
				
				mP1[i] = max(mD1[i-1] - rfgo - rfge, mP1[i-1] - rfge);
				mQ1[i] = max(mD0[i] - rdgo - rdge, mQ0[i] - rdge);
				mD1[i] = max(mD0[i-1] + ms, max(mP1[i], mQ1[i]));
				
				if (mD1[i] > score)
				{
					score = mD1[i];
				}
			}
		}
		
		swap(mP1, mP0);
		swap(mQ1, mQ0);
		swap(mD1, mD0);
		
		if (seqPos >= 0)
		{
			scores[seqPos] = score;
		}
	}
}

