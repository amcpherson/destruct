/*
 *  MatePairDelauny.h
 */

#ifndef MATEPAIRDELAUNY_H_
#define MATEPAIRDELAUNY_H_

#include "Common.h"
#include "DebugCheck.h"

#include <string>
#include <map>
#include <iostream>
#include <boost/unordered_map.hpp>

extern "C" {
#include "../external/Triangle/triangle.h"
}

using namespace boost;
using namespace std;


class MatePairDelauny
{
public:
	MatePairDelauny()
	{}
	
	void DoClustering(const MatePairVec& matePairs, IntegerTable& clusters)
	{
		vector<double> points;
		unordered_map<pair<double,double>,int> transformedToPointIndex;
		unordered_set<int> unvisited;
		double maxFragmentMean = 0;
		double maxFragmentStdDev = 0;
		for (int matePairIndex = 0; matePairIndex < matePairs.size(); matePairIndex++)
		{
			maxFragmentMean = max(maxFragmentMean, matePairs[matePairIndex].u);
			maxFragmentStdDev = max(maxFragmentStdDev, matePairs[matePairIndex].s);
			
			pair<double,double> transformed = GetTransformedPoint(matePairs[matePairIndex]);
			
			pair<unordered_map<pair<double,double>,int>::const_iterator,bool> insertResult = transformedToPointIndex.insert(make_pair(transformed,points.size()/2));
			
			if (insertResult.second)
			{
				unvisited.insert(insertResult.first->second);
				
				points.push_back(transformed.first);
				points.push_back(transformed.second);
			}
		}
		
		vector<int> edgelist;
		Triangulate(points, edgelist);
		
		double maxEdge = maxFragmentMean + 6 * maxFragmentStdDev;
		double maxEdgeSq = pow(maxEdge, 2.0);
		
		unordered_map<int,vector<int> > edges;
		for (int edgeIndex = 0; edgeIndex < edgelist.size() / 2; edgeIndex++)
		{
			int vertex1 = edgelist[edgeIndex * 2];
			int vertex2 = edgelist[edgeIndex * 2 + 1];
			
			double x1 = points[2 * vertex1];
			double y1 = points[2 * vertex1 + 1];
			double x2 = points[2 * vertex2];
			double y2 = points[2 * vertex2 + 1];
			
			double distSq = pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0);
			
			if (distSq < maxEdgeSq)
			{
				edges[vertex1].push_back(vertex2);
				edges[vertex2].push_back(vertex1);
			}
		}
		
		IntegerVec membership(points.size() / 2);
		int currentCluster = 0;
		
		vector<int> stack;
		stack.push_back(*unvisited.begin());
		
		while (!unvisited.empty())
		{
			while (!stack.empty())
			{
				int current = stack.back();
				stack.pop_back();
				
				if (unvisited.find(current) == unvisited.end())
				{
					continue;
				}
				
				unvisited.erase(current);
				
				membership[current] = currentCluster;
				
				for (vector<int>::const_iterator iter = edges[current].begin(); iter != edges[current].end(); iter++)
				{
					if (unvisited.find(*iter) != unvisited.end())
					{
						stack.push_back(*iter);
					}
				}
			}
			
			if (!unvisited.empty())
			{
				stack.push_back(*unvisited.begin());
			}
			
			currentCluster++;
		}
		
		IntegerTable newClusters(currentCluster);
		for (int matePairIndex = 0; matePairIndex < matePairs.size(); matePairIndex++)
		{
			pair<double,double> transformed = GetTransformedPoint(matePairs[matePairIndex]);
			int pointIndex = transformedToPointIndex[transformed];
			newClusters[membership[pointIndex]].push_back(matePairIndex);
		}
		
		for (IntegerTableIter clusterIter = newClusters.begin(); clusterIter != newClusters.end(); clusterIter++)
		{
			clusters.push_back(IntegerVec());
			swap(clusters.back(), *clusterIter);
		}
	}
	
private:
	pair<double,double> GetTransformedPoint(const MatePair& matePair) const
	{
		double remap1 = matePair.x + matePair.y + matePair.u;
		double remap2 = matePair.x - matePair.y;
		
		return pair<double,double>(remap1,remap2);
	}
	
	void Triangulate(const vector<double>& points, vector<int>& edges) const
	{
		if (points.size() == 0)
		{
			return;
		}
		
		if (points.size() <= 6)
		{
			for (int pointIndex1 = 0; pointIndex1 < points.size() / 2; pointIndex1++)
			{
				for (int pointIndex2 = pointIndex1 + 1; pointIndex2 < points.size() / 2; pointIndex2++)
				{
					edges.push_back(pointIndex1);
					edges.push_back(pointIndex2);
				}
			}
			
			return;
		}
		
		triangulateio triangulateIn;
		triangulateIn.numberofpoints = points.size() / 2;
		triangulateIn.numberofpointattributes = 0;
		triangulateIn.pointlist = const_cast<double*>(&points.front());
		triangulateIn.pointmarkerlist = NULL;
		triangulateIn.numberofsegments = 0;
		triangulateIn.numberofholes = 0;
		triangulateIn.numberofregions = 0;
		
		triangulateio triangulateOut;
		triangulateOut.pointlist = NULL;
		triangulateOut.pointattributelist = NULL;
		triangulateOut.pointmarkerlist = NULL;
		triangulateOut.trianglelist = NULL;
		triangulateOut.triangleattributelist = NULL;
		triangulateOut.neighborlist = NULL;
		triangulateOut.segmentlist = NULL;
		triangulateOut.segmentmarkerlist = NULL;
		triangulateOut.edgelist = NULL;
		triangulateOut.edgemarkerlist = NULL;
		
		const char* options = "pczeQ";
		triangulate(const_cast<char*>(options), &triangulateIn, &triangulateOut, NULL);
		
		edges.resize(triangulateOut.numberofedges * 2);
		copy(triangulateOut.edgelist, triangulateOut.edgelist + 2 * triangulateOut.numberofedges, edges.begin());
		
		free(triangulateOut.pointlist);
		free(triangulateOut.pointattributelist);
		free(triangulateOut.pointmarkerlist);
		free(triangulateOut.trianglelist);
		free(triangulateOut.triangleattributelist);
		free(triangulateOut.neighborlist);
		free(triangulateOut.segmentlist);
		free(triangulateOut.segmentmarkerlist);
		free(triangulateOut.edgelist);
		free(triangulateOut.edgemarkerlist);
	}
};


#endif

