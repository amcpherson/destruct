/*
 *  Algorithms.cpp
 *
 *  Created by Andrew McPherson on 11/09/12.
 *
 */

#include "Algorithms.h"
#include "BinaryMinHeap.h"

#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>


using namespace std;
using namespace boost;


void SetCover(const IntegerVecMap& sets, const DoubleMap& weights, IntegerVec& solution)
{
	BinaryMinHeap minHeap;
	
	IntegerVecMap elementSets;
	IntegerMap setSizes;
	for (IntegerVecMapConstIter setIter = sets.begin(); setIter != sets.end(); setIter++)
	{
		setSizes[setIter->first] = (int)setIter->second.size();
		
		for (IntegerVecConstIter elementIter = setIter->second.begin(); elementIter != setIter->second.end(); elementIter++)
		{
			elementSets[*elementIter].push_back(setIter->first);
			
			minHeap.Push(setIter->first, weights.find(setIter->first)->second / (double)setIter->second.size());
		}
	}
	
	IntegerSet assigned;
	while (!minHeap.Empty())
	{
		int nextSetID = minHeap.MinID();
		
		solution.push_back(nextSetID);
		
		const IntegerVec& set = sets.find(nextSetID)->second;
		
		IntegerSet alteredSets;
		for (IntegerVecConstIter elementIter = set.begin(); elementIter != set.end(); elementIter++)
		{
			if (assigned.insert(*elementIter).second)
			{
				for (IntegerVecConstIter setIter = elementSets[*elementIter].begin(); setIter != elementSets[*elementIter].end(); setIter++)
				{
					setSizes[*setIter]--;
					alteredSets.insert(*setIter);
				}
			}
		}
		
		for (IntegerSetConstIter setIter = alteredSets.begin(); setIter != alteredSets.end(); setIter++)
		{
			DebugCheck(setSizes[*setIter] >= 0);
			
			if (setSizes[*setIter] == 0)
			{
				minHeap.Remove(*setIter);
			}
			else
			{
				minHeap.IncreaseKey(*setIter, weights.find(*setIter)->second / (double)setSizes[*setIter]);
			}
		}
	}
}

void AssignInOrder(const IntegerVecMap& sets, const IntegerVec& order, IntegerVecMap& result)
{
	IntegerSet assigned;
	for (IntegerVecConstIter idIter = order.begin(); idIter != order.end(); idIter++)
	{
		const IntegerVec& set = sets.find(*idIter)->second;
		for (IntegerVecConstIter elementIter = set.begin(); elementIter != set.end(); elementIter++)
		{
			if (assigned.find(*elementIter) == assigned.end())
			{
				result[*idIter].push_back(*elementIter);
				assigned.insert(*elementIter);
			}
		}
	}
}
