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
		
		mMixWeight[alignedLength] = SAFEPARSE(double, fields[1]);
		mNBSizeTrue[alignedLength] = SAFEPARSE(double, fields[2]);
		mNBProbTrue[alignedLength] = SAFEPARSE(double, fields[3]);
		mNBSizeInvalid[alignedLength] = SAFEPARSE(double, fields[4]);
		mNBProbInvalid[alignedLength] = SAFEPARSE(double, fields[5]);
	}
}

double AlignmentProbability::ProbTrue(int alignedLength, int score) const
{
	return pdf(math::negative_binomial(mNBSizeTrue.find(alignedLength)->second, mNBProbTrue.find(alignedLength)->second), mMatchScore * alignedLength - score);
}

double AlignmentProbability::ProbInvalid(int alignedLength, int score) const
{
	return pdf(math::negative_binomial(mNBSizeInvalid.find(alignedLength)->second, mNBProbInvalid.find(alignedLength)->second), mMatchScore * alignedLength - score);
}

double AlignmentProbability::ProbFalse(int alignedLength, int score) const
{
	double probTrue = ProbTrue(alignedLength, score);
	double probInvalid = ProbInvalid(alignedLength, score);
	double mixWeight = mMixWeight.find(alignedLength)->second;
	
	return mixWeight * probTrue / (mixWeight * probTrue + (1.0 - mixWeight) * probInvalid);
}

double AlignmentProbability::Classify(int alignedLength, int score, double prior) const
{
	double probTrue = ProbTrue(alignedLength, score);
	double probInvalid = ProbInvalid(alignedLength, score);
	
	return prior * probTrue / (prior * probTrue + (1.0 - prior) * probInvalid);
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


