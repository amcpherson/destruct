/*
 *  ShortestPath.h
 */

#ifndef SHORTESTPATH_H_
#define SHORTESTPATH_H_

#include <vector>

using namespace std;

struct IChildIterator
{
	virtual ~IChildIterator() {}
	
	virtual bool Valid() = 0;
	virtual void Next() = 0;
	virtual unsigned int ChildID() = 0;
	virtual double ToChildScore() = 0;
};

struct IBreakpointGraph
{
	virtual ~IBreakpointGraph() {}
	
	virtual IChildIterator* BeginChildIteration(unsigned int vertex, double toScore) const = 0;
	virtual void EndChildIteration(unsigned int vertex) const = 0;
};

bool ShortestPath(const IBreakpointGraph* breakpointGraph, unsigned int startVertex, unsigned int endVertex, int visitLimit, double scoreLimit, vector<unsigned int>& path, double& pathScore, int& numVisited);

#endif
