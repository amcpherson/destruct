/*
 *  AlignmentProbability.cpp
 *
 */

#include "AlignmentProbability.h"
#include "Common.h"

#include <fstream>

using namespace boost;
using namespace std;


double ScoreLikelihoodCalculator::Calculate(int score) const
{
	return mExponLambda * exp(-mExponLambda * (double)(mMaxScore - score));
}

void AlignmentProbability::ReadDistributions(const string& filename, double cdfThreshold)
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
		double exponLambda = SAFEPARSE(double, fields[1]);

		int maxScore = mMatchScore * alignedLength;

		ScoreLikelihoodCalculator scoreLikelihood(maxScore, exponLambda);
		
		// Calculate a score threshold based on the cdf of the likelihood
		double normalize = 0.0;
		for (int score = mMatchScore * alignedLength; score >= -mMatchScore * alignedLength; score--)
		{
			normalize += scoreLikelihood.Calculate(score);
		}

		double cdf = 0.0;
		for (int score = mMatchScore * alignedLength; score >= -mMatchScore * alignedLength; score--)
		{
			cdf += scoreLikelihood.Calculate(score) / normalize;

			if (cdf > (1.0 - cdfThreshold))
			{
				mScoreThresholds[alignedLength] = score + 1;
				break;
			}
		}

		assert(mScoreThresholds.find(alignedLength) != mScoreThresholds.end());

		mScoreLikelihoods[alignedLength] = scoreLikelihood;
	}
}

bool AlignmentProbability::Good() const
{
	return ((mScoreLikelihoods.size() > 0) && (mScoreThresholds.size() > 0));
}

double AlignmentProbability::Likelihood(int alignedLength, int score) const
{
	return mScoreLikelihoods.find(alignedLength)->second.Calculate(score);
}

bool AlignmentProbability::AboveThreshold(int alignedLength, int score) const
{
	return score >= mScoreThresholds.find(alignedLength)->second;
}

int AlignmentProbability::GetMinAlignedLength() const
{
	return min_element(mScoreThresholds.begin(), mScoreThresholds.end())->first;
}

void AlignmentPosterior::AppendAlignment(int readEnd, int score)
{
	double readLikelihood = mAlignmentProbability.Likelihood(mAlignedLengths[readEnd], score);

	mSumReadLikelihood[readEnd] += readLikelihood;

	mReadEnd.push_back(readEnd);
	mReadLikelihood.push_back(readLikelihood);
	mConcordantLikelihood.push_back(0.0);
}

void AlignmentPosterior::AppendAlignmentWithMate(int readEnd, int score, int mateScore)
{
	int mateEnd = OtherReadEnd(readEnd);

	double readLikelihood = mAlignmentProbability.Likelihood(mAlignedLengths[readEnd], score);
	double mateLikelihood = mAlignmentProbability.Likelihood(mAlignedLengths[mateEnd], mateScore);

	mSumReadLikelihood[readEnd] += readLikelihood;
	mSumReadLikelihood[mateEnd] += mateLikelihood;
	mSumConcordantLikelihood += readLikelihood * mateLikelihood;

	mReadEnd.push_back(readEnd);
	mReadLikelihood.push_back(readLikelihood);
	mConcordantLikelihood.push_back(readLikelihood * mateLikelihood);
}

double AlignmentPosterior::Posterior(int index)
{
	double N = mPriorDiscordant * mReadLikelihood[index] * mSumReadLikelihood[OtherReadEnd(mReadEnd[index])] + 
	           (1.0 - 2.0 * mPriorDiscordant) * mConcordantLikelihood[index];
	double Z = mPriorDiscordant * mSumReadLikelihood[0] * mSumReadLikelihood[1] + 
	           (1.0 - 2.0 * mPriorDiscordant) * mSumConcordantLikelihood;
	
	return N / Z;
}

double AlignmentPosterior::PosteriorConcordant()
{
	double N = (1.0 - mPriorDiscordant) * mSumConcordantLikelihood;
	double Z = mPriorDiscordant * mSumReadLikelihood[0] * mSumReadLikelihood[1] + 
	           (1.0 - 2.0 * mPriorDiscordant) * mSumConcordantLikelihood;
	
	return N / Z;
}


