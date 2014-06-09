/*
 *  AlignmentProbability.h
 *
 */

#ifndef ALIGNMENTPROBABILITY_H_
#define ALIGNMENTPROBABILITY_H_

#include "Common.h"


class ScoreLikelihoodCalculator
{
public:
	ScoreLikelihoodCalculator() : mMaxScore(0), mExponLambda(0.0) {}
	ScoreLikelihoodCalculator(int maxScore, double exponLambda) : mMaxScore(maxScore), mExponLambda(exponLambda) {}

	double Calculate(int score) const;

private:
	int mMaxScore;
	double mExponLambda;
};


class AlignmentProbability
{
public:
	AlignmentProbability(int matchScore) : mMatchScore(matchScore) {}
	
	void ReadDistributions(const string& filename, double cdfThreshold);
	
	double Likelihood(int alignedLength, int score) const;
	bool AboveThreshold(int alignedLength, int score) const;

private:
	int mMatchScore;
	unordered_map<int,ScoreLikelihoodCalculator> mScoreLikelihoods;
	unordered_map<int,int> mScoreThresholds;
};


class AlignmentPosterior
{
public:
	AlignmentPosterior(const AlignmentProbability& alignmentProbability, double priorDiscordant, const vector<int>& alignedLengths)
	: mAlignmentProbability(alignmentProbability), mPriorDiscordant(priorDiscordant), mAlignedLengths(alignedLengths)
	{
		mSumReadLikelihood[0] = 0.0;
		mSumReadLikelihood[1] = 0.0;
		mSumConcordantLikelihood = 0.0;
	}
	
	// Add alignment to the posterior calculation
	void AppendAlignment(int readEnd, int score);

	// Add alignment to the posterior calculation including
	// score of a potentially concordant mate
	void AppendAlignmentWithMate(int readEnd, int score, int mateScore);
	
	// Posterior probability an alignment location is correct
	double Posterior(int index);

	// Posterior probability any concordant paired alignment is correct
	double PosteriorConcordant();
	
private:
	double mPriorDiscordant;
	const AlignmentProbability& mAlignmentProbability;
	vector<int> mAlignedLengths;

	double mSumReadLikelihood[2];
	double mSumConcordantLikelihood;

	vector<double> mReadEnd;
	vector<double> mReadLikelihood;
	vector<double> mConcordantLikelihood;
};


#endif


