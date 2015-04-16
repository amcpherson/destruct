/*
 *  MatePairGibbs.cpp
 *
 *  Created by Andrew McPherson.
 *
 */

#include "MatePairGibbs.h"
#include "DebugCheck.h"

#include "asa136.H"

#include <iostream>
#include <vector>
#include <list>
#include <numeric>
#include <limits>
#include "asa241.H"
#include "BinaryMinHeap.h"
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace boost;
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

void PrintMatePairs(const MatePairVec& matePairs)
{
	DoubleVec x;
	DoubleVec y;
	DoubleVec u;
	DoubleVec s;
	for (int i = 0; i < matePairs.size(); i++)
	{
		x.push_back(matePairs[i].x);
		y.push_back(matePairs[i].y);
		u.push_back(matePairs[i].u);
		s.push_back(matePairs[i].s);
	}
	
	cout << "X="; printvec(x.begin(), x.end());
	cout << "Y="; printvec(y.begin(), y.end());
	cout << "U="; printvec(u.begin(), u.end());
	cout << "S="; printvec(s.begin(), s.end());
}

double phi(double x)
{
	return cdf(math::normal(),x);
	
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
	
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
	
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
	
    return 0.5*(1.0 + sign*y);
}

MatePairGibbs::MatePairGibbs()
{
}

class MatePairGrid
{
public:
	explicit MatePairGrid(const MatePairVec& matePairs, int spacing) : mSpacing(spacing)
	{
		for (int matePairIndex = 0; matePairIndex < matePairs.size(); matePairIndex++)
		{
			int gridX = (int)matePairs[matePairIndex].x / mSpacing;
			int gridY = (int)matePairs[matePairIndex].y / mSpacing;
			
			for (int i = -1; i <= 1; i++)
			{
				for (int j = -1; j <= 1; j++)
				{
					int k = (i + 3 * j + 9) % 9;
					
					list<int>& gridList = mGrid.insert(make_pair(pair<int,int>(gridX+i,gridY+j),list<int>())).first->second;
					
					mAdjacentList[k].push_back(&gridList);
					mAdjacentActive[k].push_back(gridList.end());
				}
			}
		}
	}
	
	void Activate(int matePairIndex)
	{
		for (int i = 0; i < 9; i++)
		{
			DebugCheck(mAdjacentActive[i][matePairIndex] == mAdjacentList[i][matePairIndex]->end());
			
			mAdjacentList[i][matePairIndex]->push_front(matePairIndex);
			mAdjacentActive[i][matePairIndex] = mAdjacentList[i][matePairIndex]->begin();
		}
	}
	
	void DeActivate(int matePairIndex)
	{
		for (int i = 0; i < 9; i++)
		{
			DebugCheck(mAdjacentActive[i][matePairIndex] != mAdjacentList[i][matePairIndex]->end());
			
			mAdjacentList[i][matePairIndex]->erase(mAdjacentActive[i][matePairIndex]);
			mAdjacentActive[i][matePairIndex] = mAdjacentList[i][matePairIndex]->end();
		}
	}
	
	const list<int>& AdjacentActivated(int matePairIndex)
	{
		return *mAdjacentList[0][matePairIndex];
	}
	
private:
	typedef unordered_map<pair<int,int>,list<int> > Grid;
	typedef vector<list<int>*> ListPtrVec;
	typedef vector<list<int>::iterator> ListIterVec;
	
	int mSpacing;
	Grid mGrid;
	ListPtrVec mAdjacentList[9];
	ListIterVec mAdjacentActive[9];
};

class ClusterSet
{
public:
	ClusterSet(int numElements, MatePairGrid& matePairGrid) : mClusters(numElements), mIterators(numElements), mSize(0), mMatePairGrid(matePairGrid)
	{
	}
	
	void Add(int matePairIndex, int clusterIndex)
	{
		if (matePairIndex == clusterIndex)
		{
			DebugCheck(mClusters[matePairIndex] == 0);
			
			mClusters[matePairIndex] = new list<int>();
			mSize++;
			mIterators[matePairIndex] = mClusters[matePairIndex]->insert(mClusters[matePairIndex]->end(), matePairIndex);
			mMatePairGrid.Activate(clusterIndex);
		}
		else
		{
			mClusters[matePairIndex] = mClusters[clusterIndex];
			mIterators[matePairIndex] = mClusters[matePairIndex]->insert(mClusters[matePairIndex]->end(), matePairIndex);
		}
	}
	
	void Remove(int matePairIndex)
	{
		bool deactivated = false;
		if (mIterators[matePairIndex] == mClusters[matePairIndex]->begin())
		{
			mMatePairGrid.DeActivate(matePairIndex);
			deactivated = true;
		}
		
		mClusters[matePairIndex]->erase(mIterators[matePairIndex]);
		
		if (mClusters[matePairIndex]->empty())
		{
			delete mClusters[matePairIndex];
			mSize--;
		}
		else if (deactivated)
		{
			mMatePairGrid.Activate(mClusters[matePairIndex]->front());
		}
		
		mClusters[matePairIndex] = 0;
	}
	
	int Size()
	{
		return mSize;
	}
	
	void Save(vector<int>& membership) const
	{
		membership.resize(mClusters.size());
		for (int matePairIndex = 0; matePairIndex < mClusters.size(); matePairIndex++)
		{
			membership[matePairIndex] = mClusters[matePairIndex]->front();
		}
	}
	
private:
	typedef vector<list<int>*> ListPtrVec;
	typedef vector<list<int>::iterator> ListIterVec;
	
	ListPtrVec mClusters;
	ListIterVec mIterators;
	int mSize;
	MatePairGrid& mMatePairGrid;
};

/*
 pe_cluster_prob <- function(X, Y, U, S)
 {
 Norm = (S * exp(-U^2/(2*S^2)) + U * sqrt(2*pi) * pnorm(U/S)) / sqrt(2*pi)
 W = U - max(X) + X - max(Y) + Y
 alpha = sum(1/S^2)
 beta = sum(W/S^2)
 gamma = sum(W^2/S^2)
 a1 = prod(1/(Norm*sqrt(2*pi)*S))
 a2 = exp(-(1/2)*(gamma-beta^2/alpha))/alpha
 a3 = exp(-beta^2/(2*alpha)) + beta * sqrt(2*pi) * pnorm(beta/sqrt(alpha)) / sqrt(alpha)
 pA = a1*a2*a3
 return(pA)
 }*/

double ClusterProb(const DoubleVec& X, const DoubleVec& Y, const DoubleVec& U, const DoubleVec& S)
{
	double maxX = X[0];
	double maxY = Y[0];
	for (int i = 0; i < X.size(); i++)
	{
		maxX = max(maxX, X[i]);
		maxY = max(maxY, Y[i]);
	}
	
	double alpha = 0.0;
	double beta = 0.0;
	double gamma = 0.0;
	double zprod = 1.0;
	for (int i = 0; i < X.size(); i++)
	{
		double norm = (S[i] * exp(-pow(U[i],2.0)/(2.0*pow(S[i],2.0))) + U[i] * sqrt(2.0*M_PI) * phi(U[i]/S[i])) / sqrt(2.0*M_PI);
		double w = U[i] - maxX + X[i] - maxY + Y[i];
		alpha += 1.0/pow(S[i],2.0);
		beta += w/pow(S[i],2.0);
		gamma += pow(w,2.0)/pow(S[i],2.0);
		zprod *= 1.0/(norm*sqrt(2.0*M_PI)*S[i]);
	}
	
	double a2 = exp(-0.5*(gamma-pow(beta,2.0)/alpha))/alpha;
	double a3 = exp(-pow(beta,2.0)/(2.0*alpha)) + beta * sqrt(2.0*M_PI) * phi(beta/sqrt(alpha)) / sqrt(alpha);
	
	return zprod * a2 * a3;
}

double LogClusterProb(const MatePairVec& matePairs)
{
	double maxX = matePairs[0].x;
	double maxY = matePairs[0].y;
	for (int i = 0; i < matePairs.size(); i++)
	{
		maxX = max(maxX, matePairs[i].x);
		maxY = max(maxY, matePairs[i].y);
	}
	
	double alpha = 0.0;
	double beta = 0.0;
	double gamma = 0.0;
	double zprod = 0.0;
	for (int i = 0; i < matePairs.size(); i++)
	{
		double norm = (matePairs[i].s * exp(-pow(matePairs[i].u,2.0)/(2.0*pow(matePairs[i].s,2.0))) + matePairs[i].u * sqrt(2.0*M_PI) * phi(matePairs[i].u/matePairs[i].s)) / sqrt(2.0*M_PI);
		double w = matePairs[i].u - maxX + matePairs[i].x - maxY + matePairs[i].y;
		alpha += 1.0/pow(matePairs[i].s,2.0);
		beta += w/pow(matePairs[i].s,2.0);
		gamma += pow(w,2.0)/pow(matePairs[i].s,2.0);
		zprod -= log(norm*sqrt(2.0*M_PI)*matePairs[i].s);
	}
	
	double a2 = -0.5*(gamma-pow(beta,2.0)/alpha) - log(alpha);
	double a3 = log(exp(-pow(beta,2.0)/(2.0*alpha)) + beta * sqrt(2.0*M_PI) * phi(beta/sqrt(alpha)) / sqrt(alpha));
	
	return zprod + a2 + a3;
}

double dpprior(const IntegerVec& memberships)
{
	double alpha = 1.0;
	double dpprior = 1.0;
	unordered_map<int,int> membercounts;
	for (IntegerVecConstIter iter = memberships.begin(); iter != memberships.end(); iter++)
	{
		pair<unordered_map<int,int>::iterator,bool> result = membercounts.insert(make_pair(*iter,0));
		if (result.second)
		{
			dpprior *= alpha;
		}
		else
		{
			dpprior *= result.first->second;
		}
		result.first->second++;
	}
	return dpprior;
}

double calclogf(double x)
{
	double step = 0.2;
	double minx = -30.0;
	double maxx = 30.0;
	
	if (x < minx)
	{
		return -2.0*log(-x);
	}
	else if (x > maxx)
	{
		return log(x) + log(sqrt(2.0*M_PI)) + 0.5*x*x;
	}
	
	static DoubleVec xs;
	static DoubleVec fs;
	
	if (xs.empty())
	{
		for (long double x = minx; x < maxx + 2.0*step; x+= step)
		{
			double f = log(1.0 + x/math::hazard(math::normal_distribution<long double>(),-x));
			
			if (f < -10.0 || f > 5000.0)
			{
				cerr << "Fatal error: Could not correctly calculate log f(x)" << endl;
				cerr << "Value calculated was " << f << endl;
				exit(1);
			}
			
			xs.push_back(x);
			fs.push_back(f);
		}
	}
	
	int lookup = (int)((x - minx) / step);
	double remainder = (x - minx) - (double)lookup * step;
	
	DebugCheck(lookup >= 0);
	DebugCheck(lookup+1 < fs.size());
	
	double f = fs[lookup] + remainder * (fs[lookup+1] - fs[lookup]) / step;
	
	return f;
}

class ClusterStats
{
public:
	ClusterStats(vector<int>& indX, vector<int>& indY) : mLogZProd(0.0), mAlpha(0.0), mBeta(0.0), mGamma(0.0), mLogProb(0.0), mXs(indX), mYs(indY) {}
	
	double X() const
	{
		return mXs.MaxPriority();
	}
	
	double Y() const
	{
		return mYs.MaxPriority();
	}
	
	double Size() const
	{
		return (double)mXs.Size();
	}
	
	double LogProb() const
	{
		return mLogProb;
	}
	
	double LogPosteriorPredictive(double norm, const MatePair& mp) const
	{
		DebugCheck(mXs.Size() > 0);
		
		double q = mXs.MaxPriority() - max(mXs.MaxPriority(), mp.x);
		double r = mYs.MaxPriority() - max(mYs.MaxPriority(), mp.y);
		
		double w = mp.u - max(mXs.MaxPriority(),mp.x) + mp.x - max(mYs.MaxPriority(),mp.y) + mp.y;
		
		double newLogZProd = mLogZProd - log(norm*sqrt(2.0*M_PI)*mp.s);
		double newGamma = mGamma + 2.0*(q+r)*mBeta + pow(q+r,2.0)*mAlpha + (w*w/(mp.s*mp.s));
		double newBeta = mBeta + (q+r)*mAlpha + (w/(mp.s*mp.s));
		double newAlpha = mAlpha + (1.0/(mp.s*mp.s));
		
		return LogClusterProbability(newLogZProd,newAlpha,newBeta,newGamma) - mLogProb;
	}
	
	void Add(int id, double norm, const MatePair& mp)
	{
		double q = 0.0;
		double r = 0.0;
		if (mXs.Size() > 0)
		{
			q = mXs.MaxPriority() - max(mXs.MaxPriority(), mp.x);
			r = mYs.MaxPriority() - max(mYs.MaxPriority(), mp.y);
		}
		
		mXs.Push(id,mp.x);
		mYs.Push(id,mp.y);
		
		double w = mp.u - mXs.MaxPriority() + mp.x - mYs.MaxPriority() + mp.y;
		
		mLogZProd = mLogZProd - log(norm*sqrt(2.0*M_PI)*mp.s);
		mGamma = mGamma + 2.0*(q+r)*mBeta + pow(q+r,2.0)*mAlpha + (w*w/(mp.s*mp.s));
		mBeta = mBeta + (q+r)*mAlpha + (w/(mp.s*mp.s));
		mAlpha = mAlpha + (1.0/(mp.s*mp.s));
		
		mLogProb = LogClusterProbability(mLogZProd,mAlpha,mBeta,mGamma);
	}
	
	void Remove(int id, double norm, const MatePair& mp)
	{
		DebugCheck(mXs.Size() > 1);
		
		double w = mp.u - mXs.MaxPriority() + mp.x - mYs.MaxPriority() + mp.y;
		
		mXs.Remove(id);
		mYs.Remove(id);
		
		double q = mXs.MaxPriority() - max(mXs.MaxPriority(), mp.x);
		double r = mYs.MaxPriority() - max(mYs.MaxPriority(), mp.y);
		
		mLogZProd = mLogZProd + log(norm*sqrt(2.0*M_PI)*mp.s);
		mAlpha = mAlpha - (1.0/(mp.s*mp.s));
		mBeta = mBeta - (q+r)*mAlpha - (w/(mp.s*mp.s));
		mGamma = mGamma - 2.0*(q+r)*mBeta - pow(q+r,2.0)*mAlpha - (w*w/(mp.s*mp.s));
		
		mLogProb = LogClusterProbability(mLogZProd,mAlpha,mBeta,mGamma);
	}
	
private:
	double LogClusterProbability(double logzprod, double alpha, double beta, double gamma) const
	{
		double a = -0.5*gamma - log(alpha);
		double b = calclogf(beta/sqrt(alpha));
		
		return logzprod + a + b;
	}
	
	double mLogZProd;
	double mAlpha;
	double mBeta;
	double mGamma;
	double mLogProb;
	BinaryMaxHeapSharedIndex mXs;
	BinaryMaxHeapSharedIndex mYs;
};


double NormConstant(double u, double s)
{
	return (s * exp(-pow(u,2.0)/(2.0*pow(s,2.0))) + u * sqrt(2.0*M_PI) * phi(u/s)) / sqrt(2.0*M_PI);
}

void MatePairGibbs::DoClustering(const MatePairVec& matePairs, IntegerTable& clusters)
{
	if (matePairs.size() == 0)
	{
		return;
	}
	else if (matePairs.size() == 1)
	{
		clusters.push_back(IntegerVec(1, 0));
		return;
	}
	
	double numStdDev = 6.0;
	double prior = 1e-10;
	double alpha = 1.0;
	int numIterations = 100;
	int numStable = 5;
	int burnin = 2;
	
	double maxRange = 0.0;
	DoubleVec norms;
	IntegerVec shuffledIndices;
	for (int matePairIndex = 0; matePairIndex < matePairs.size(); matePairIndex++)
	{
		maxRange = max(maxRange, matePairs[matePairIndex].u + numStdDev * matePairs[matePairIndex].s);
		norms.push_back(NormConstant(matePairs[matePairIndex].u, matePairs[matePairIndex].s));
		shuffledIndices.push_back(matePairIndex);
	}
	
	vector<ClusterStats*> clusterStats((int)matePairs.size(), (ClusterStats*)0);
	vector<int> sharedHeapIndexX((int)matePairs.size(), -1);
	vector<int> sharedHeapIndexY((int)matePairs.size(), -1);
	MatePairGrid matePairGrid(matePairs, (int)maxRange);
	ClusterSet clusterSet((int)matePairs.size(), matePairGrid);
	
	double partitionLog = 0.0;
	double clustersLogProb = 0.0;
	double bestClusteringScore = -numeric_limits<double>::max();
	IntegerVec bestMembership;
	int stableCount = 0;
	for (int iteration = 0; iteration < numIterations; iteration++)
	{
		random_shuffle(shuffledIndices.begin(), shuffledIndices.end());
		
		for (IntegerVecConstIter shuffleIter = shuffledIndices.begin(); shuffleIter != shuffledIndices.end(); shuffleIter++)
		{
			int matePairIndex = *shuffleIter;
			
			// Remove from previous cluster
			if (clusterStats[matePairIndex])
			{
				// Remove previous probability from total 
				clustersLogProb -= clusterStats[matePairIndex]->LogProb();
				
				// Decrement previous contribution to partition probability
				if (clusterStats[matePairIndex]->Size() > 1)
				{
					partitionLog -= log(clusterStats[matePairIndex]->Size() - 1);
				}
				
				if (clusterStats[matePairIndex]->Size() == 1)
				{
					// Delete soon to be empty cluster
					delete clusterStats[matePairIndex];
				}
				else
				{
					// Remove from cluster
					clusterStats[matePairIndex]->Remove(matePairIndex, norms[matePairIndex], matePairs[matePairIndex]);
					
					// Add current probability of new cluster to total
					clustersLogProb += clusterStats[matePairIndex]->LogProb();
				}
				
				// This mate pair no longer has a cluster
				clusterSet.Remove(matePairIndex);
				clusterStats[matePairIndex] = 0;
			}
			
			// Calculate membership probabilities for adjacent clusters
			IntegerVec clusterIDs;
			DoubleVec logClusterProb;
			double maxLogClusterProb = -numeric_limits<double>::max();
			const list<int>& adjacent = matePairGrid.AdjacentActivated(matePairIndex);
			for (list<int>::const_iterator clusterIDIter = adjacent.begin(); clusterIDIter != adjacent.end(); clusterIDIter++)
			{
				ClusterStats* cluster = clusterStats[*clusterIDIter];
				clusterIDs.push_back(*clusterIDIter);
				logClusterProb.push_back(log(cluster->Size()) + cluster->LogPosteriorPredictive(norms[matePairIndex], matePairs[matePairIndex]));
				maxLogClusterProb = max(maxLogClusterProb, logClusterProb.back());
			}
			logClusterProb.push_back(log(alpha * prior));
			maxLogClusterProb = max(maxLogClusterProb, logClusterProb.back());
			
			// Normalization
			double clusterProbNorm = 0.0;
			for (int probIndex = 0; probIndex < logClusterProb.size(); probIndex++)
			{
				clusterProbNorm += exp(logClusterProb[probIndex] - maxLogClusterProb);
			}
			
			// Randomly select cluster 
			double probAccum = 0.0;
			double randProb = (double)rand() / (double)RAND_MAX;
			int sampleIndex = 0;
			for (int probIndex = 0; probIndex < logClusterProb.size(); probIndex++)
			{
				double clusterProb = exp(logClusterProb[probIndex] - maxLogClusterProb) / clusterProbNorm;
				probAccum += clusterProb;
				if (randProb <= probAccum)
				{
					sampleIndex = probIndex;
					break;
				}
			}
			
			// Assign to next cluster
			if (sampleIndex < clusterIDs.size())
			{
				// Assign to existing cluster
				clusterSet.Add(matePairIndex, clusterIDs[sampleIndex]);
				clusterStats[matePairIndex] = clusterStats[clusterIDs[sampleIndex]];
			}
			else
			{
				// Create new cluster
				clusterSet.Add(matePairIndex, matePairIndex);
				clusterStats[matePairIndex] = new ClusterStats(sharedHeapIndexX, sharedHeapIndexY);
			}
			
			// Remove previous probability from total 
			clustersLogProb -= clusterStats[matePairIndex]->LogProb();
			
			// Add the mate pair to its new cluster
			clusterStats[matePairIndex]->Add(matePairIndex, norms[matePairIndex], matePairs[matePairIndex]);
			
			// Add current probability to total
			clustersLogProb += clusterStats[matePairIndex]->LogProb();
			
			// Increment contribution to partition probability
			if (clusterStats[matePairIndex]->Size() > 1)
			{
				partitionLog += log(clusterStats[matePairIndex]->Size() - 1);
			}
			
			// Allow a burnin period
			if (iteration > burnin)
			{
				// Check clustering score
				double k = (double)clusterSet.Size();
				double clusteringScore = k * log(alpha*prior) + partitionLog + clustersLogProb;
				
				if (clusteringScore > bestClusteringScore)
				{
					bestClusteringScore = clusteringScore;
					clusterSet.Save(bestMembership);
					stableCount = 0;
				}
				else
				{
					stableCount++;
				}
				
				if (stableCount >= numStable)
				{
					break;
				}
			}
		}
	}
	
	unordered_map<int,int> clusterIDToCluster;
	for (int matePairIndex = 0; matePairIndex < matePairs.size(); matePairIndex++)
	{
		pair<unordered_map<int,int>::const_iterator,bool> insertResult = clusterIDToCluster.insert(make_pair(bestMembership[matePairIndex],clusters.size()));
		
		if (insertResult.second)
		{
			clusters.push_back(IntegerVec());
		}
		
		clusters[insertResult.first->second].push_back(matePairIndex);
	}
}

void RandomizedClusterStatsTests()
{
	srand(20);
	vector<int> sharedHeapIndexX(100000, -1);
	vector<int> sharedHeapIndexY(100000, -1);
	ClusterStats teststats(sharedHeapIndexX, sharedHeapIndexY);
	vector<int> ids;
	vector<MatePair> matepairs;
	vector<double> norms;
	int nextID = 0;
	for (int i = 0; i < 100000; i++)
	{
		double removetest = (double)rand() / (double)RAND_MAX;
		if (removetest < min(1.0,(double)matepairs.size()/20.0))
		{
			cout << "remove" << endl;
			
			int toremove = rand() % matepairs.size();
			
			teststats.Remove(ids[toremove], norms[toremove], matepairs[toremove]);
			
			ids.erase(ids.begin()+toremove);
			matepairs.erase(matepairs.begin()+toremove);
			norms.erase(norms.begin()+toremove);
		}
		else if (matepairs.size() > 0)
		{
			cout << "add" << endl;
			
			MatePair mp;
			mp.x = 200.0 * (double)rand() / (double)RAND_MAX;
			mp.y = 200.0 * (double)rand() / (double)RAND_MAX;
			mp.u = 100.0 + 50.0 * (double)rand() / (double)RAND_MAX;
			mp.s = 25.0 + 25.0 * (double)rand() / (double)RAND_MAX;
			
			double beforeProb = LogClusterProb(matepairs);
			
			ids.push_back(nextID);
			matepairs.push_back(mp);
			norms.push_back(NormConstant(mp.u, mp.s));
			
			double afterProb = LogClusterProb(matepairs);
			double postPredTrue = afterProb - beforeProb;
			
			double postPred = teststats.LogPosteriorPredictive(norms.back(), matepairs.back());
			
			cout << postPredTrue << "\t" << postPred << "\t" << postPredTrue - postPred << endl;
			
			teststats.Add(ids.back(), norms.back(), matepairs.back());
			
			nextID++;
		}
		else
		{
			cout << "addfirst" << endl;
			
			MatePair mp;
			mp.x = 200.0 * (double)rand() / (double)RAND_MAX;
			mp.y = 200.0 * (double)rand() / (double)RAND_MAX;
			mp.u = 100.0 + 50.0 * (double)rand() / (double)RAND_MAX;
			mp.s = 25.0 + 25.0 * (double)rand() / (double)RAND_MAX;
			
			ids.push_back(nextID);
			matepairs.push_back(mp);
			norms.push_back(NormConstant(mp.u, mp.s));
			
			teststats.Add(ids.back(), norms.back(), matepairs.back());
			
			nextID++;
		}
		
		cout << LogClusterProb(matepairs) << "\t" << teststats.LogProb() << "\t" << LogClusterProb(matepairs) - teststats.LogProb() << endl;
		
		if ((LogClusterProb(matepairs) - teststats.LogProb()) > 0.000001)
		{
			cout << "failed" << endl;
			exit(1);
		}
	}
}

double runif(double a, double b)
{
	double r = (double)rand() / (double)RAND_MAX;
	return r * (b - a) + a;
}

int rdiscrete(const DoubleVec& p)
{
	double r = runif(0.0,1.0);
	double accum = 0.0;
	for (int i = 0; i < p.size(); i++)
	{
		accum += p[i];
		if (r <= accum)
		{
			return i;
		}
	}
	return p.size() - 1;
}

double rnorm(double u, double s)
{
	mt19937 rng; 
	rng.seed(rand());
	normal_distribution<double> normaldist(u,s);		
	boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > vnorm(rng, normaldist);
	return vnorm();
}

void GenerateRandomMatePairs(int seed, int K, int N, double range, MatePairVec& matePairs, IntegerTable& clusters)
{
	srand(seed);
	
	DoubleVec clusterX;
	DoubleVec clusterY;
	DoubleVec clusterW;
	double weightSum = 0.0;
	for (int i = 0; i < K; i++)
	{
		clusterX.push_back(runif(0.0,range));
		clusterY.push_back(runif(0.0,range));
		clusterW.push_back(runif(0.0,1.0));
		weightSum += clusterW.back();
	}
	
	for (int i = 0; i < K; i++)
	{
		clusterW[i] /= weightSum;
	}
	
	clusters.clear();
	clusters.resize(K);
	for (int i = 0; i < N; i++)
	{
		int clusterIndex = rdiscrete(clusterW);
		clusters[clusterIndex].push_back(i);
		matePairs.push_back(MatePair());
		matePairs.back().u = runif(200.0,500.0);
		matePairs.back().s = runif(25.0,75.0);
		double length = rnorm(matePairs.back().u,matePairs.back().s);
		while (length <= 20.0)
		{
			length = rnorm(matePairs.back().u,matePairs.back().s);
		}
		double offset = runif(0.0,length);
		matePairs.back().x = clusterX[clusterIndex] - offset;
		matePairs.back().y = clusterY[clusterIndex] - (length - offset);
	}
}

double Evaluate(const IntegerTable& known, const IntegerTable& predicted)
{
	unordered_set<IntegerPair> pairs;
	for (IntegerTableConstIter clusterIter = known.begin(); clusterIter != known.end(); clusterIter++)
	{
		for (IntegerVecConstIter element1Iter = clusterIter->begin(); element1Iter != clusterIter->end(); element1Iter++)
		{
			for (IntegerVecConstIter element2Iter = element1Iter + 1; element2Iter != clusterIter->end(); element2Iter++)
			{
				pairs.insert(IntegerPair(*element1Iter,*element2Iter));
				pairs.insert(IntegerPair(*element2Iter,*element1Iter));
			}
		}
	}
	
	double totalPairCount = pairs.size() / 2;
	
	double correctPairCount = 0.0;
	for (IntegerTableConstIter clusterIter = predicted.begin(); clusterIter != predicted.end(); clusterIter++)
	{
		for (IntegerVecConstIter element1Iter = clusterIter->begin(); element1Iter != clusterIter->end(); element1Iter++)
		{
			for (IntegerVecConstIter element2Iter = element1Iter + 1; element2Iter != clusterIter->end(); element2Iter++)
			{
				if (pairs.find(IntegerPair(*element1Iter,*element2Iter)) != pairs.end())
				{
					correctPairCount++;
				}
			}
		}
	}
	
	return correctPairCount / totalPairCount;
}

void RandomizedClusteringTest()
{
	MatePairGibbs matePairGibbs;
	
	MatePairVec matePairs;
	IntegerTable simulatedClusters;
	GenerateRandomMatePairs(3, 20, 1000, 10000, matePairs, simulatedClusters);
	
	IntegerTable clusters;
	matePairGibbs.DoClustering(matePairs, clusters);
	
	PrintMatePairs(matePairs);
	
	cout << Evaluate(simulatedClusters, clusters) << endl;
}

#if 0

int main()
{
	MatePairGibbs matePairGibbs;
	
	/*
	MatePairVec matePairs;
	for (int i = 0; i < sizeof(X)/sizeof(double); i++)
	{
		matePairs.push_back(MatePair());
		matePairs.back().x = X[i];
		matePairs.back().y = Y[i];
		matePairs.back().u = U[i];
		matePairs.back().s = S[i];
	}*/
	
	MatePairVec matePairs;
	IntegerTable simulatedClusters;
	GenerateRandomMatePairs(3, 20, 50000, 5000, matePairs, simulatedClusters);
	
	IntegerTable clusters;
	matePairGibbs.DoClustering(matePairs, clusters);
	
	Print(simulatedClusters);
	Print(clusters);
	
	cout << Evaluate(simulatedClusters, clusters) << endl;
}
#endif


