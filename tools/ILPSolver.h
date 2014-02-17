/*
 *  ILPSolver.h
 *
 *  Created by Andrew McPherson on 10-09-07.
 *
 */

#ifndef ILPSOLVER_H_
#define ILPSOLVER_H_

#include "Common.h"

#include <vector>
#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

class ILPSolver
{
public:
	ILPSolver(const IntegerVecMap& dnaRnaOverlap, const IntegerVecMap& rnaDnaOverlap, double w1, double w2, double w3, double w4, int timeout);
	double Solve(const IntegerVecMap& dnaClusters, const IntegerVecMap& rnaClusters, IntegerVecMap& dnaSolution, IntegerVecMap& rnaSolution, bool reduce);

private:
	IntegerVecMap mDnaRnaOverlap;
	IntegerVecMap mRnaDnaOverlap;
	
	double mW1;
	double mW2;
	double mW3;
	double mW4;
	int mTimeout;
};

#endif

