/*
 *  Algorithms.h
 *
 *  Created by Andrew McPherson on 11/09/12.
 *
 */

#ifndef ALGORITHMS_H_
#define ALGORITHMS_H_

#include "BinaryMinHeap.h"

#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>


using namespace std;
using namespace boost;


template<typename TElemType>
void SetCover(const vector<vector<TElemType> >& sets, const vector<double>& weights, vector<int>& solution)
{
	BinaryMinHeap minHeap;
	
	unordered_map<TElemType,vector<int> > elementSets;

	vector<int> setSizes(sets.size(), 0);

	for (int setIdx = 0; setIdx < (size_t)sets.size(); setIdx++)
	{
		const vector<TElemType>& set = sets[setIdx];

		setSizes[setIdx] = (int)set.size();
		
		for (typename vector<TElemType>::const_iterator elementIter = set.begin(); elementIter != set.end(); elementIter++)
		{
			elementSets[*elementIter].push_back(setIdx);
			
			minHeap.Push(setIdx, weights[setIdx] / (double)set.size());
		}
	}
	
	unordered_set<TElemType> assigned;
	while (!minHeap.Empty())
	{
		int setIdx = minHeap.MinID();
		
		solution.push_back(setIdx);
		
		const vector<TElemType>& set = sets[setIdx];
		
		unordered_set<int> alteredSets;

		for (typename vector<TElemType>::const_iterator elementIter = set.begin(); elementIter != set.end(); elementIter++)
		{
			if (assigned.insert(*elementIter).second)
			{
				for (vector<int>::const_iterator setIter = elementSets[*elementIter].begin(); setIter != elementSets[*elementIter].end(); setIter++)
				{
					setSizes[*setIter]--;
					alteredSets.insert(*setIter);
				}
			}
		}
		
		for (unordered_set<int>::const_iterator setIter = alteredSets.begin(); setIter != alteredSets.end(); setIter++)
		{
			DebugCheck(setSizes[*setIter] >= 0);
			
			if (setSizes[*setIter] == 0)
			{
				minHeap.Remove(*setIter);
			}
			else
			{
				minHeap.IncreaseKey(*setIter, weights[*setIter] / (double)setSizes[*setIter]);
			}
		}
	}
}

template<typename TElemType>
void AssignInOrder(const vector<vector<TElemType> >& sets, const vector<int>& order, vector<vector<TElemType> >& result)
{
	result = vector<vector<TElemType> >(sets.size(), vector<TElemType>());

	unordered_set<TElemType> assigned;
	for (vector<int>::const_iterator idxIter = order.begin(); idxIter != order.end(); idxIter++)
	{
		const vector<TElemType>& set = sets[*idxIter];

		for (typename vector<TElemType>::const_iterator elementIter = set.begin(); elementIter != set.end(); elementIter++)
		{
			if (assigned.find(*elementIter) == assigned.end())
			{
				result[*idxIter].push_back(*elementIter);
				assigned.insert(*elementIter);
			}
		}
	}
}

#endif
