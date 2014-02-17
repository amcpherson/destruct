/*
 *  Algorithms.h
 *
 *  Created by Andrew McPherson on 11/09/12.
 *
 */

#ifndef ALGORITHMS_H_
#define ALGORITHMS_H_

#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "Common.h"

using namespace std;
using namespace boost;

void SetCover(const IntegerVecMap& sets, const DoubleMap& weights, IntegerVec& solution);
void AssignInOrder(const IntegerVecMap& sets, const IntegerVec& order, IntegerVecMap& result);

#endif
