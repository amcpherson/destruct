/*
 *  MatePairGibbs.h
 *
 *  Created by Andrew McPherson.
 *
 */

#ifndef MATEPAIRGIBBS_H_
#define MATEPAIRGIBBS_H_

#include "Common.h"

#include <vector>

using namespace std;

class MatePairGibbs
{
public:
	MatePairGibbs();
	
	void DoClustering(const MatePairVec& matePairs, IntegerTable& clusters);
	
private:
};

#endif
