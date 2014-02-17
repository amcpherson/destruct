/*
 *  MatePairEM.cpp
 *
 *  Created by Andrew McPherson.
 *
 */

#include "MatePairEM2.h"
#include "DebugCheck.h"

#include "asa136.H"

#include <iostream>
#include <vector>
#include <numeric>
#include <limits>
#include "asa241.H"

using namespace std;


template <typename t>
void printvec(t begin, t end)
{
	cout.precision(30);
	cout << "c(";
	
	t iter = begin;
	while (iter != end)
	{
		cout << *iter;
		iter++;
		
		if (iter != end)
		{
			cout << ", ";
		}
	}
	
	cout << ")" << endl;
}

MatePairEM::MatePairEM()
{
	mKMeansIter = 1000;
	mLambda = 0.1;
	mTolerance = 0.001;
	mKMax = 20;
	mMaxN = 0;
}

void MatePairEM::UpdateCapacity(int n)
{
	if (n < mMaxN)
	{
		return;
	}
	
	mR.resize(mKMax);
	mRXO.resize(mKMax);
	mRYO.resize(mKMax);
	mW.resize(mKMax);
	mA.resize(mKMax);
	mB.resize(mKMax);
	mTempExponents.resize(mKMax);
	
	for (int j = 0; j < mKMax; j++)
	{
		mR[j].resize(n);
		mRXO[j].resize(n);
		mRYO[j].resize(n);
		mTempExponents[j].resize(n);
	}
	
	mTempSX.resize(n);
	mTempSY.resize(n);
	
	mTempCX.resize(4 * n + 1);
	mTempCY.resize(4 * n + 1);
	mTempCS.resize(4 * n + 1);
	mTempPartial.resize(4 * n + 1);
	
	mMaxN = n;
}

void MatePairEM::CreateMatePair(int id, const Region& alignment1, int strand1, const Region& alignment2, int strand2, double fragmentLengthMean, double fragmentLengthStdDev, MatePair& matePair) const
{
	matePair.id = id;
	StrandRemap(alignment1, strand1, matePair.alignment1);
	StrandRemap(alignment2, strand2, matePair.alignment2);
	matePair.fragmentLengthMean = fragmentLengthMean;
	matePair.fragmentLengthStdDev = fragmentLengthStdDev;
}

void MatePairEM::StrandRemap(const Region& region, int strand, Region& remapped) const
{
	Region remappedTemp;
	
	remappedTemp.start = (strand == PlusStrand) ? region.start : -region.end;
	remappedTemp.end = (strand == PlusStrand) ? region.end : -region.start;
	
	remapped = remappedTemp;
}

double MatePairEM::PairProbability(double x, double y, double u, double s, double a, double b) const
{
	return normalpdf(a + b - x - y, u, s) * exp(-mLambda * max(0.0, x - a) - mLambda * max(0.0, y - b));
}

double MatePairEM::LogLikelihood() const
{
	for (int i = 0; i < mN; i++)
	{
		for (int j = 0; j < mK; j++)
		{
			mTempExponents[j][i] = -0.5 * pow((mA[j] + mB[j] - mX[i] - mY[i] - mU[i]) / mS[i], 2.0) - mLambda * max(0.0, mX[i] - mA[j]) - mLambda * max(0.0, mY[i] - mB[j]);
		}
	}
	
	double LL = 0.0;
	for (int i = 0; i < mN; i++)
	{
		double maxexp = mTempExponents[0][i];
		for (int j = 1; j < mK; j++)
		{
			maxexp = max(maxexp, mTempExponents[j][i]);
		}
		
		double sum = 0.0;
		for (int j = 0; j < mK; j++)
		{
			sum += mW[j] * exp(mTempExponents[j][i] - maxexp);
		}
				
		if (sum == 0.0)
		{
			LL = -numeric_limits<double>::max();
			break;
		}
		
		LL = LL + log(sum) + maxexp;
	}
	
	return LL;
}

void MatePairEM::UpdateResponsibilities()
{
	for (int i = 0; i < mN; i++)
	{
		for (int j = 0; j < mK; j++)
		{
			mTempExponents[j][i] = -0.5 * pow((mA[j] + mB[j] - mX[i] - mY[i] - mU[i]) / mS[i], 2.0) - mLambda * max(0.0, mX[i] - mA[j]) - mLambda * max(0.0, mY[i] - mB[j]);
		}
	}
	
	for (int i = 0; i < mN; i++)
	{
		int iXO = mToXO[i];
		int iYO = mToYO[i];
		
		double maxexp = mTempExponents[0][i];
		for (int j = 1; j < mK; j++)
		{
			maxexp = max(maxexp, mTempExponents[j][i]);
		}
		
		double norm = 0.0;
		for (int j = 0; j < mK; j++)
		{
			norm += mW[j] * exp(mTempExponents[j][i] - maxexp);
		}
		
		DebugCheck(norm != 0.0);
		for (int j = 0; j < mK; j++)
		{
			mR[j][i] = mW[j] * exp(mTempExponents[j][i] - maxexp) / norm;
			
			mRXO[j][iXO] = mR[j][i];
			mRYO[j][iYO] = mR[j][i];
		}
	}
}

void MatePairEM::UpdateMixWeights()
{
	for (int j = 0; j < mK; j++)
	{
		double NK = accumulate(mR[j].begin(), mR[j].begin() + mN, 0.0);
		mW[j] = NK / mN;
	}
}

bool MatePairEM::MaxLikelihood(const vector<double>& R, const vector<double>& RXO, const vector<double>& RYO, double& a, double& b) const
{
	partial_sum(RXO.begin(), RXO.begin() + mN, mTempSX.begin());
	partial_sum(RYO.begin(), RYO.begin() + mN, mTempSY.begin());

	int i = 0;
	int j = 0;
	int l = 0;
	
	mTempCX[l] = mXO[0];
	mTempCY[l] = mYO[0];
	mTempCS[l] = 0.0;
	l++;
	
	while (i < mN && j < mN)
	{
		if (i + 1 < mN && mXO[i] == mXO[i+1])
		{
			i++;
			continue;
		}
		
		if (j + 1 < mN && mYO[j] == mYO[j+1])
		{
			j++;
			continue;
		}

		if (mTempSX[i] == mTempSY[j])
		{
			mTempCX[l] = mXO[i];
			mTempCY[l] = mYO[j];
			mTempCS[l] = mTempSX[i];
			l++;
			
			if (i + 1 < mN && j + 1 < mN)
			{
				mTempCX[l] = mXO[i+1];
				mTempCY[l] = mYO[j+1];
				mTempCS[l] = mTempSX[i];
				l++;
			}
			
			i++;
			j++;
		}
		else if (mTempSX[i] < mTempSY[j])
		{
			mTempCX[l] = mXO[i];
			mTempCY[l] = mYO[j];
			mTempCS[l] = mTempSX[i];
			l++;
			
			if (i + 1 < mN)
			{
				mTempCX[l] = mXO[i+1];
				mTempCY[l] = mYO[j];
				mTempCS[l] = mTempSX[i];
				l++;
			}
			
			i++;
		}
		else
		{
			mTempCX[l] = mXO[i];
			mTempCY[l] = mYO[j];
			mTempCS[l] = mTempSY[j];
			l++;
			
			if (j + 1 < mN)
			{
				mTempCX[l] = mXO[i];
				mTempCY[l] = mYO[j+1];
				mTempCS[l] = mTempSY[j];
				l++;
			}
			
			j++;
		}
	}
	
	double RXYU = 0.0;
	double NKS = 0.0;
	for (int i = 0; i < mN; i++)
	{
		RXYU += R[i] * (mX[i] + mY[i] + mU[i]) / pow(mS[i], 2.0);
		NKS += R[i] / pow(mS[i], 2.0);
	}
	
	if (NKS == 0.0)
	{
		return false;
	}
	
	for (int i = 0; i < l; i++)
	{
		mTempPartial[i] = RXYU - NKS * (mTempCX[i] + mTempCY[i]) + mLambda * mTempCS[i];
	}
	
	int minindex = 0;
	while (minindex < l)
	{
		if (mTempPartial[minindex] > 0)
		{
			break;
		}
		
		minindex++;
	}
	
	double aplusb = (RXYU + mLambda * mTempCS[minindex]) / NKS;
	
	if (minindex == 0)
	{
		double min_a = mTempCX[minindex];
		double max_a = aplusb - mTempCY[minindex];
		a = 0.5 * (min_a + max_a);
		b = aplusb - a;
	}
	else if (mTempCS[minindex] != mTempCS[minindex-1])
	{
		a = mTempCX[minindex];
		b = mTempCY[minindex];
	}
	else
	{
		double min_a = max(mTempCX[minindex], aplusb - mTempCY[minindex-1]);
		double max_a = min(mTempCX[minindex-1], aplusb - mTempCY[minindex]);
		a = 0.5 * (min_a + max_a);
		b = aplusb - a;
	}
	
	return true;
}

bool MatePairEM::SelectKKZ(int k, vector<double>& A, vector<double>& B)
{
	DebugCheck(k <= mN);
	
	A.clear();
	B.clear();
	
	double L2max = mX[0] * mY[0];
	int iL2max = 0;
	for (int i = 1; i < mN; i++)
	{
		double L2 = mX[i] * mY[i];
		
		if (L2 > L2max)
		{
			iL2max = i;
			L2max = L2;
		}
	}
	
	A.push_back(mX[iL2max]);
	B.push_back(mY[iL2max]);
		
	while ((int)A.size() < k)
	{
		vector<double> DistMin(mN);
		for (int i = 0; i < mN; i++)
		{
			double minDist = pow(mX[i] - A[0], 2.0) + pow(mY[i] - B[0], 2.0);
			for (int j = 1; j < (int)A.size(); j++)
			{
				double dist = pow(mX[i] - A[j], 2.0) + pow(mY[i] - B[j], 2.0);
				minDist = min(minDist, dist);
			}
			
			DistMin[i] = minDist;
		}
		
		double DistsMax = DistMin[0];
		int iDistsMax = 0;
		for (int i = 0; i < mN; i++)
		{
			if (DistMin[i] > DistsMax)
			{
				DistsMax = DistMin[i];
				iDistsMax = i;
			}
		}
		
		if (DistsMax == 0.0)
		{
			return false;
		}
		
		A.push_back(mX[iDistsMax]);
		B.push_back(mY[iDistsMax]);
	}
	
	return true;
}

void MatePairEM::MStep()
{
	for (int j = 0; j < mK; j++)
	{
		double a;
		double b;
		if (MaxLikelihood(mR[j], mRXO[j], mRYO[j], a, b))
		{
			mA[j] = a;
			mB[j] = b;
		}
	}
}

bool MatePairEM::ExpectationMaximization(double& logLikelihood)
{
	DebugCheck(mK <= mN);
	
	double lastLikelihood;
	bool lastLikelihoodValid = false;
	while (true)
	{
		MStep();
		
		UpdateMixWeights();
		
		double likelihood = LogLikelihood();
		if (lastLikelihoodValid && abs(likelihood - lastLikelihood) < mTolerance)
		{
			break;
		}
		
		if (lastLikelihoodValid && likelihood == -numeric_limits<double>::max())
		{
			return false;
		}
		
		DebugCheck(!lastLikelihoodValid || (likelihood / lastLikelihood < 1.0000001));
		
		lastLikelihood = likelihood;
		lastLikelihoodValid = true;

		UpdateResponsibilities();
	}
	
	logLikelihood = lastLikelihood;
	
	return true;
}

void MatePairEM::KExpectationMaximization(double& logLikelihood)
{
	mK = 1;
	
	for (int i = 0; i < mN; i++)
	{
		mR[0][i] = 1.0;
		mRXO[0][i] = 1.0;
		mRYO[0][i] = 1.0;
	}
	
	double previousBIC = numeric_limits<double>::max();
	double previousLogLikelihood;
	vector<double> previousW;
	vector<double> previousA;
	vector<double> previousB;	
	while (true)
	{
		ExpectationMaximization(logLikelihood);
		
		double BIC = -2.0 * logLikelihood + mK * 2.0 * log(mN);
		
		if (BIC > previousBIC)
		{
			mK--;
			
			mW = previousW;
			mA = previousA;
			mB = previousB;
			
			UpdateResponsibilities();
			
			logLikelihood = previousLogLikelihood;
			
			return;
		}
		
		if (false)//mK == mKMax)
		{
			cout << "X = ";
			printvec(mX.begin(),mX.begin()+mN);
			cout << "Y = ";
			printvec(mY.begin(),mY.begin()+mN);
			cout << "U = ";
			printvec(mU.begin(),mU.begin()+mN);
			cout << "S = ";
			printvec(mS.begin(),mS.begin()+mN);
			cout << endl;
		}
		
		if (mK >= mKMax || mK >= mN)
		{
			return;
		}
		
		previousBIC = BIC;
		previousLogLikelihood = logLikelihood;
		previousW = mW;
		previousA = mA;
		previousB = mB;
		
		double minModelProb;
		int minModelProbIndex;
		for (int i = 0; i < mN; i++)
		{
			double modelProb = 0.0;
			for (int j = 1; j < mK; j++)
			{
				modelProb += mW[j] * PairProbability(mX[i], mY[i], mU[i], mS[i], mA[j], mB[j]);
			}
			
			if (i == 0 || modelProb < minModelProb)
			{
				minModelProb = modelProb;
				minModelProbIndex = i;
			}
		}
		
		int minModelProbXO = mToXO[minModelProbIndex];
		int minModelProbYO = mToYO[minModelProbIndex];
		
		for (int j = 0; j < mK; j++)
		{
			mR[j][minModelProbIndex] = 0.0;
			mRXO[j][minModelProbXO] = 0.0;
			mRYO[j][minModelProbYO] = 0.0;
		}
		
		for (int i = 0; i < mN; i++)
		{
			mR[mK][i] = 0.0;
			mRXO[mK][i] = 0.0;
			mRYO[mK][i] = 0.0;
		}
		
		mR[mK][minModelProbIndex] = 1.0;
		mRXO[mK][minModelProbXO] = 1.0;
		mRYO[mK][minModelProbYO] = 1.0;
		
		for (int j = 0; j < mK; j++)
		{
			mW[j] *= (1 - 1/((double)mK + 1));
		}
		mW[mK] = 1/((double)mK + 1);
		
		mK++;
		
		MStep();
		
		UpdateResponsibilities();
	}
}

struct MatePairInfo
{
	int i;
	double x;
	double y;
	double u;
	double s;
};

typedef vector<MatePairInfo> MatePairInfoVec;

bool XGreaterThan(const MatePairInfo& m1, const MatePairInfo& m2)
{
	return (m1.x > m2.x);
}

bool YGreaterThan(const MatePairInfo& m1, const MatePairInfo& m2)
{
	return (m1.y > m2.y);
}





//int xlist[]={-993917, -972421, -960630, -993910, -993918, -994491, -993907, -993915, -993907, -960618, -993906, -993901, -988339, -960618, -987217, -993906, -993918, -993899, -993906, -988317, -993932, -993906, -993910, -993915, -993910, -960618, -993918, -994133, -993925, -991084, -988321, -972421, -991084, -988321, -993920, -993917, -993903, -993906};
//int ylist[]={991037, 975193, 947665, 991035, 991049, 992004, 991035, 991036, 991035, 947653, 991036, 991043, 985018, 947653, 989891, 991036, 991034, 991026, 991036, 985599, 991054, 991036, 991026, 991035, 991035, 947653, 991035, 992071, 991035, 987200, 985012, 975193, 987200, 985012, 991039, 991035, 991039, 991036};

//void MatePairEM::DoClustering(const MatePairVec& matePairs, IntegerTable& clusters)
//{
//	if (matePairs.size() < mMinClusterSize)
//	{
//		return;
//	}
//	
//	mN = sizeof(xlist) / sizeof(int);
//	
//	MatePairInfoVec info(mN);
//	for (int matePairIndex = 0; matePairIndex < mN; matePairIndex++)
//	{
//		info[matePairIndex].i = matePairIndex;
//		info[matePairIndex].x = xlist[matePairIndex];
//		info[matePairIndex].y = ylist[matePairIndex];
//		info[matePairIndex].u = mFragmentMean - 100;
//	}

void MatePairEM::DoClustering(const MatePairVec& matePairs, IntegerTable& clusters, DoubleTable& likelihoods)
{
	mN = matePairs.size();
	
	UpdateCapacity(mN);
	
	MatePairInfoVec info(mN);
	for (int matePairIndex = 0; matePairIndex < (int)matePairs.size(); matePairIndex++)
	{
		info[matePairIndex].i = matePairIndex;
		info[matePairIndex].x = matePairs[matePairIndex].alignment1.end;
		info[matePairIndex].y = matePairs[matePairIndex].alignment2.end;
		info[matePairIndex].u = matePairs[matePairIndex].fragmentLengthMean - Length(matePairs[matePairIndex].alignment1) - Length(matePairs[matePairIndex].alignment2);
		info[matePairIndex].s = matePairs[matePairIndex].fragmentLengthStdDev;
	}
	
	mX.resize(mN);
	mY.resize(mN);
	mU.resize(mN);
	mS.resize(mN);
	for (int matePairIndex = 0; matePairIndex < (int)info.size(); matePairIndex++)
	{
		mX[matePairIndex] = info[matePairIndex].x;
		mY[matePairIndex] = info[matePairIndex].y;
		mU[matePairIndex] = info[matePairIndex].u;
		mS[matePairIndex] = info[matePairIndex].s;
	}
	
	sort(info.begin(), info.end(), XGreaterThan);
	
	mXO.resize(mN);
	mToXO.resize(mN);
	for (int sortIndex = 0; sortIndex < (int)info.size(); sortIndex++)
	{
		mXO[sortIndex] = info[sortIndex].x;
		mToXO[info[sortIndex].i] = sortIndex;
	}
	
	sort(info.begin(), info.end(), YGreaterThan);
	
	mYO.resize(mN);
	mToYO.resize(mN);
	for (int sortIndex = 0; sortIndex < (int)info.size(); sortIndex++)
	{
		mYO[sortIndex] = info[sortIndex].y;
		mToYO[info[sortIndex].i] = sortIndex;
	}
	
	double logLikelihood;
	KExpectationMaximization(logLikelihood);
	
	IntegerTable newClusters(mK);
	DoubleTable newLikelihoods(mK);
	for (int i = 0; i < mN; i++)
	{
		double maxResp = -1.0;
		int jMaxResp = 0;
		for (int j = 0; j < mK; j++)
		{
			if (mR[j][i] > maxResp)
			{
				maxResp = mR[j][i];
				jMaxResp = j;
			}
		}
		newClusters[jMaxResp].push_back(i);
		newLikelihoods[jMaxResp].push_back(PairProbability(mX[i], mY[i], mU[i], mS[i], mA[jMaxResp], mB[jMaxResp]));
	}
	
	for (int j = 0; j < mK; j++)
	{
		//if (newClusters[j].empty())
		if (newClusters[j].size() > 100)
		{
/*
			cerr << "Warning: EM Error" << endl;
			printvec(mX.begin(), mX.end());
			printvec(mY.begin(), mY.end());
			printvec(mU.begin(), mU.end());
			printvec(mS.begin(), mS.end());
			cerr << mK << endl;
			printvec(mW.begin(), mW.end());
			exit(1);
			continue;
*/		}
		
		clusters.push_back(IntegerVec());
		swap(clusters.back(), newClusters[j]);
		
		likelihoods.push_back(DoubleVec());
		swap(likelihoods.back(), newLikelihoods[j]);
	}
}

