/*
 *  ShortestPath.cpp
 */

#include "ShortestPath.h"
#include "BinaryMinHeap.h"

#include <vector>
#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

bool ShortestPath(const IBreakpointGraph* breakpointGraph, unsigned int startVertex, unsigned int endVertex, int visitLimit, double scoreLimit, vector<unsigned int>& path, double& pathScore, int& numVisited)
{
	BinaryMinHeap toVisit;
	unordered_map<unsigned int,double> toScore;
	unordered_map<unsigned int,unsigned int> fromIDs;
	unordered_map<unsigned int,vector<unsigned int> > parents;
	unordered_map<unsigned int,IChildIterator*> childIterators;
	unordered_set<unsigned int> visited;
	
	toVisit.Push(startVertex, 0.0);
	toScore[startVertex] = 0.0;
	
	numVisited = 0;
	while (!toVisit.Empty())
	{
		numVisited++;
		if (numVisited > visitLimit)
		{
			break;
		}
		
		unsigned int currentVertex = toVisit.Pop();
		
		if (currentVertex == endVertex)
		{
			break;
		}
		
		if (toScore[currentVertex] > scoreLimit)
		{
			break;
		}
		
		// Start iteration for children
		childIterators[currentVertex] = breakpointGraph->BeginChildIteration(currentVertex, toScore[currentVertex]);
		
		// Create vector to increment
		vector<unsigned int> toIncrement;
		swap(parents[currentVertex], toIncrement);
		parents.erase(currentVertex);
		
		// Add first child of currentVertex
		toIncrement.push_back(currentVertex);
		
		// Set as visited
		visited.insert(currentVertex);
		
		// Iterate over vertices requiring increment
		for (vector<unsigned int>::const_iterator vertexIter = toIncrement.begin(); vertexIter != toIncrement.end(); vertexIter++)
		{
			if (childIterators.find(*vertexIter) == childIterators.end())
			{
				continue;
			}
			
			IChildIterator* childIter = childIterators.find(*vertexIter)->second;
			
			// Increment parents until an unvisited vertex found or finished
			while (childIter->Valid() && visited.find(childIter->ChildID()) != visited.end())
			{
				childIter->Next();
			}
			
			// Check for end of iteration
			if (childIter->Valid())
			{
				// Update best score for new child
				if (toScore.find(childIter->ChildID()) == toScore.end() || childIter->ToChildScore() < toScore[childIter->ChildID()])
				{
					toVisit.Push(childIter->ChildID(), childIter->ToChildScore());
					toScore[childIter->ChildID()] = childIter->ToChildScore();
					fromIDs[childIter->ChildID()] = *vertexIter;
				}
				
				// Update list of parents, regardless of score
				parents[childIter->ChildID()].push_back(*vertexIter);
				
				// Increment for next update
				childIter->Next();
			}
			
			if (!childIter->Valid())
			{
				// Remove finished iterators
				breakpointGraph->EndChildIteration(*vertexIter);
				childIterators.erase(*vertexIter);
			}
		}
	}
	
	// Create output path
	if (fromIDs.find(endVertex) != fromIDs.end())
	{
		unsigned int traverseVertex = endVertex;
		while (fromIDs.find(traverseVertex) != fromIDs.end() && traverseVertex != startVertex)
		{
			path.push_back(traverseVertex);
			traverseVertex = fromIDs[traverseVertex];
		}
		path.push_back(traverseVertex);
		
		reverse(path.begin(), path.end());
		
		pathScore = toScore[endVertex];
		
		return true;
	}
	
	return false;
}

