/*
 *  AffineGapDP.h
 *
 *  Created by Andrew McPherson on 18/07/12.
 *
 */

#ifndef AFFINEGAPDP_H_
#define AFFINEGAPDP_H_

#include "Common.h"
#include "Matrix.h"

#include <string>
#include <vector>

using namespace std;

class AffineGapDP
{
public:
	AffineGapDP();
	
	int EndToEndAlign(const string& reference, const string& sequence, const string& quality);
	void PartialAlign(const string& reference, const string& sequence, const string& quality, IntegerVec& scores);

private:
	vector<int> mD0;
	vector<int> mD1;
	vector<int> mP0;
	vector<int> mP1;
	vector<int> mQ0;
	vector<int> mQ1;
};

#endif

