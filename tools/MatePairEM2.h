/*
 *  MatePairEM.h
 *
 *  Created by Andrew McPherson.
 *
 */

#ifndef MATEPAIREM_H_
#define MATEPAIREM_H_

#include "Common.h"

#include <vector>

using namespace std;

class MatePairEM
{
public:
	struct MatePair
	{
		int id;
		Region alignment1;
		Region alignment2;
		double fragmentLengthMean;
		double fragmentLengthStdDev;
	};
	
	typedef vector<MatePair> MatePairVec;
	
	MatePairEM();
	
	void CreateMatePair(int id, const Region& alignment1, int strand1, const Region& alignment2, int strand2, double fragmentLengthMean, double fragmentLengthStdDev, MatePair& matePair) const;
	void DoClustering(const MatePairVec& matePairs, IntegerTable& clusters, DoubleTable& likelihoods);
	
private:
	void UpdateCapacity(int n);
	void StrandRemap(const Region& region, int strand, Region& remapped) const;
	
	double PairProbability(double x, double y, double u, double s, double a, double b) const;
	double LogLikelihood() const;
	void UpdateResponsibilities();
	bool MaxLikelihood(const vector<double>& R, const vector<double>& RXO, const vector<double>& RYO, double& a, double& b) const;
	void UpdateMixWeights();
	bool SelectKKZ(int k, vector<double>& A, vector<double>& B);
	bool ExpectationMaximization(double& ll);
	void KExpectationMaximization(double& ll);
	void MStep();

	double mLambda;
	int mKMeansIter;
	double mTolerance;
	int mKMax;
	int mMaxN;
	
	int mN;
	int mK;
	
	vector<double> mX;
	vector<double> mY;
	vector<double> mU;
	vector<double> mS;

	vector<double> mW;
	vector<double> mA;
	vector<double> mB;
	vector<vector<double> > mR;
	
	vector<double> mXO;
	vector<double> mYO;
	vector<vector<double> > mRXO;
	vector<vector<double> > mRYO;
	vector<int> mToXO;
	vector<int> mToYO;
	
	mutable vector<vector<double> > mTempExponents;
	mutable vector<double> mTempSX;
	mutable vector<double> mTempSY;
	mutable vector<double> mTempCX;
	mutable vector<double> mTempCY;
	mutable vector<double> mTempCS;	
	mutable vector<double> mTempPartial;
};

#endif
