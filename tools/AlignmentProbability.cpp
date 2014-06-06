/*
 *  AlignmentProbability.cpp
 *
 */

#include "AlignmentProbability.h"
#include "Parsers.h"
#include "Common.h"

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/negative_binomial.hpp>

using namespace boost;
using namespace std;


void AlignmentProbability::ReadDistributions(const string& filename)
{
	ifstream distFile(filename.c_str());
	
	if (!distFile.good())
	{
		cerr << "Error: Unable to open file " << filename << endl;
		exit(1);
	}
	
	StringVec fields;
	while (ReadTSV(distFile, fields))
	{
		int alignedLength = SAFEPARSE(int, fields[0]);
		
		mNBSizeTrue[alignedLength] = SAFEPARSE(double, fields[2]);
		mNBProbTrue[alignedLength] = SAFEPARSE(double, fields[3]);

		// Simulations produce unrealistic curves
		// As a work around, a score that is better than the mode will get the probability at the mode
		vector<pair<double,int> > probTrues;
		for (int x = 0; x < mMatchScore * alignedLength; x++)
		{
			double probTrue = pdf(math::negative_binomial(mNBSizeTrue[alignedLength], mNBProbTrue[alignedLength]), x);
			probTrues.push_back(make_pair(probTrue, x));
		}
		mProbTrueMode[alignedLength] = *max_element(probTrues.begin(), probTrues.end());;
	}
}

double AlignmentProbability::ProbTrue(int alignedLength, int score) const
{
	// Simulations work around, see above
	if (mMatchScore * alignedLength - score < mProbTrueMode.find(alignedLength)->second.second)
	{
		return mProbTrueMode.find(alignedLength)->second.first;
	}

	return pdf(math::negative_binomial(mNBSizeTrue.find(alignedLength)->second, mNBProbTrue.find(alignedLength)->second), mMatchScore * alignedLength - score);
}

void AlignmentPosterior::Initialize(const AlignmentProbability* alignmentProbability, int alignedLength)
{
	mAlignmentProbability = alignmentProbability;
	mAlignedLength = alignedLength;
	mSumProbTrue = 0.0;
	mMaxScore = 0;
}

void AlignmentPosterior::AddAlignment(int score)
{
	double probTrue = mAlignmentProbability->ProbTrue(mAlignedLength, score);
	
	mSumProbTrue += probTrue;
	mMaxScore = max(mMaxScore, score);
}

double AlignmentPosterior::MaxPosterior()
{
	double probTrue = mAlignmentProbability->ProbTrue(mAlignedLength, mMaxScore);
	
	return probTrue / mSumProbTrue;
}

double AlignmentPosterior::Posterior(int score)
{
	double probTrue = mAlignmentProbability->ProbTrue(mAlignedLength, score);
	
	return probTrue / mSumProbTrue;
}


