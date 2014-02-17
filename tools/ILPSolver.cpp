/*
 *  ILPSolver.cpp
 *
 *  Created by Andrew McPherson on 10-09-07.
 *
 */

#include "ILPSolver.h"
#include "Common.h"
#include "DebugCheck.h"

#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>

#include <algorithm>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;
using namespace boost;

typedef unordered_set<int> IntegerSet;
typedef unordered_set<int>::const_iterator IntegerSetIter;
typedef unordered_set<IntegerPair> IntPairSet;
typedef unordered_set<IntegerPair>::const_iterator IntPairSetIter;

int hash32shiftmult(int key)
{
	int c2=0x27d4eb2d; // a prime or an odd constant
	key = (key ^ 61) ^ (key >> 16);
	key = key + (key << 3);
	key = key ^ (key >> 4);
	key = key * c2;
	key = key ^ (key >> 15);
	return key;
}

void Print(const IntegerMap& weights)
{
	for (IntegerMapConstIter weightIter = weights.begin(); weightIter != weights.end(); weightIter++)
	{
		cout << weightIter->first << "\t" << weightIter->second << endl;
	}
}

void Print(const IntegerPairVec& overlaps)
{
	for (IntegerPairVecConstIter overlapIter = overlaps.begin(); overlapIter != overlaps.end(); overlapIter++)
	{
		cout << overlapIter->first << "\t" << overlapIter->second << endl;
	}
}

void Sort(IntegerVecMap& clusters)
{
	for (IntegerVecMapIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		sort(clusterIter->second.begin(), clusterIter->second.end());
	}
}

int Count(IntegerVecMap& clusters)
{
	int numFragments = 0;
	for (IntegerVecMapIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		numFragments += clusterIter->second.size();
	}
	return numFragments;
}

void Erase(IntegerSet& redundant, IntegerVecMap& clusters)
{
	IntegerVecMap filteredClusters;
	
	for (IntegerVecMapIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		IntegerVec filtered;
		for (IntegerVecIter fragmentIter = clusterIter->second.begin(); fragmentIter != clusterIter->second.end(); fragmentIter++)
		{
			if (redundant.find(*fragmentIter) == redundant.end())
			{
				filtered.push_back(*fragmentIter);
			}
		}
		
		if (filtered.size() > 0)
		{
			filteredClusters[clusterIter->first] = filtered;
		}
		
		clusterIter->second.clear();
	}
	
	swap(clusters, filteredClusters);
}

void RemoveRedundantClusters(IntegerVecMap& clusters, IntegerVecMap& overlaps)
{
	IntegerVecMap clusterSum;
	for (IntegerVecMapIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		int sum = 0;
		for (IntegerVecIter fragmentIter = clusterIter->second.begin(); fragmentIter != clusterIter->second.end(); fragmentIter++)
		{
			sum += hash32shiftmult(*fragmentIter);
		}
		clusterSum[sum].push_back(clusterIter->first);
	}
	
	IntegerSet redundant;
	for (IntegerVecMapIter clustersIter = clusterSum.begin(); clustersIter != clusterSum.end(); clustersIter++)
	{
		for (IntegerVecIter cluster1Iter = clustersIter->second.begin(); cluster1Iter != clustersIter->second.end(); cluster1Iter++)
		{
			int clusterID1 = *cluster1Iter;
			const IntegerVec& cluster1 = clusters[clusterID1];
			const IntegerVec& overlap1 = overlaps[clusterID1];
			
			for (IntegerVecIter cluster2Iter = clustersIter->second.begin(); cluster2Iter != clustersIter->second.end(); cluster2Iter++)
			{
				int clusterID2 = *cluster2Iter;
				const IntegerVec& cluster2 = clusters[clusterID2];
				const IntegerVec& overlap2 = overlaps[clusterID2];
				
				bool clustersEqual = cluster1.size() == cluster2.size() && equal(cluster1.begin(), cluster1.end(), cluster2.begin());
				bool overlapsEqual = overlap1.size() == overlap2.size() && equal(overlap1.begin(), overlap1.end(), overlap2.begin());
				
				if (clustersEqual && overlapsEqual && clusterID1 < clusterID2)
				{
					redundant.insert(clusterID1);
				}
			}	
		}
	}
	
	for (IntegerSetIter redundantIter = redundant.begin(); redundantIter != redundant.end(); redundantIter++)
	{
		clusters.erase(*redundantIter);
	}
}

void RemoveRedundantFragments(IntegerVecMap& fragments, IntegerVecMap& clusters)
{
	IntegerVecMap fragmentSum;
	for (IntegerVecMapIter fragmentIter = fragments.begin(); fragmentIter != fragments.end(); fragmentIter++)
	{
		int sum = 0;
		for (IntegerVecIter clusterIter = fragmentIter->second.begin(); clusterIter != fragmentIter->second.end(); clusterIter++)
		{
			sum += hash32shiftmult(*clusterIter);
		}
		fragmentSum[sum].push_back(fragmentIter->first);
	}
	
	IntegerVecMap equivalent;
	IntegerSet foundEquiv;
	for (IntegerVecMapIter fragmentsIter = fragmentSum.begin(); fragmentsIter != fragmentSum.end(); fragmentsIter++)
	{
		for (IntegerVecIter fragment1Iter = fragmentsIter->second.begin(); fragment1Iter != fragmentsIter->second.end(); fragment1Iter++)
		{
			int fragmentID1 = *fragment1Iter;
			const IntegerVec& clusters1 = fragments[fragmentID1];
			
			if (foundEquiv.find(fragmentID1) != foundEquiv.end())
			{
				continue;
			}
			
			for (IntegerVecIter fragment2Iter = fragmentsIter->second.begin(); fragment2Iter != fragmentsIter->second.end(); fragment2Iter++)
			{
				int fragmentID2 = *fragment2Iter;
				const IntegerVec& clusters2 = fragments[fragmentID2];
				
				if (fragmentID1 == fragmentID2)
				{
					continue;
				}
				
				bool fragmentsEqual = clusters1.size() == clusters2.size() && equal(clusters1.begin(), clusters1.end(), clusters2.begin());
				
				if (fragmentsEqual)
				{
					equivalent[fragmentID1].push_back(fragmentID2);
					foundEquiv.insert(fragmentID2);
				}
			}	
		}
	}

	IntegerSet redundant;
	for (IntegerVecMapIter equivIter = equivalent.begin(); equivIter != equivalent.end(); equivIter++)
	{
		int fragmentID = equivIter->first;
		int numClusters = fragments[fragmentID].size();
		int numFragments = equivIter->second.size() + 1;
		int numRedundant = numFragments - numClusters;
		
		int madeRedundant = 0;
		for (IntegerVecIter fragmentIter = equivIter->second.begin(); fragmentIter != equivIter->second.end(); fragmentIter++)
		{
			if (madeRedundant == numRedundant)
			{
				break;
			}
			
			redundant.insert(*fragmentIter);

			madeRedundant++;
		}
	}
	
	Erase(redundant, clusters);
}

ILPSolver::ILPSolver(const IntegerVecMap& dnaRnaOverlap, const IntegerVecMap& rnaDnaOverlap, double w1, double w2, double w3, double w4, int timeout)
: mDnaRnaOverlap(dnaRnaOverlap), mRnaDnaOverlap(rnaDnaOverlap), mW1(w1), mW2(w2), mW3(w3), mW4(w4), mTimeout(timeout)
{
	Sort(mDnaRnaOverlap);
	Sort(mRnaDnaOverlap);
}

double ILPSolver::Solve(const IntegerVecMap& dnaClusters, const IntegerVecMap& rnaClusters, IntegerVecMap& dnaSolution, IntegerVecMap& rnaSolution, bool reduce)
{
	IntegerVecMap dnaClustersR = dnaClusters;
	IntegerVecMap rnaClustersR = rnaClusters;
	
	Sort(dnaClustersR);
	Sort(rnaClustersR);
		
	IntegerVecMap dnaFragments;
	IntegerVecMap rnaFragments;
	Transpose(dnaClustersR, dnaFragments);
	Transpose(rnaClustersR, rnaFragments);
	Sort(dnaFragments);
	Sort(rnaFragments);

	if (reduce)
	{	
		RemoveRedundantClusters(dnaClustersR, mDnaRnaOverlap);
		RemoveRedundantClusters(rnaClustersR, mRnaDnaOverlap);
		
		RemoveRedundantFragments(dnaFragments, dnaClustersR);
		RemoveRedundantFragments(rnaFragments, rnaClustersR);
	}	
		
	dnaFragments.clear();
	rnaFragments.clear();	
	Transpose(dnaClustersR, dnaFragments);
	Transpose(rnaClustersR, rnaFragments);
	Sort(dnaFragments);
	Sort(rnaFragments);
	
	glp_prob* lp = glp_create_prob();
	glp_set_prob_name(lp, "sample");
	glp_set_obj_dir(lp, GLP_MIN);
	
	// Count number of columns required
	int numColumns = dnaClustersR.size() + rnaClustersR.size();
	for (IntegerVecMapIter clusterIter = dnaClustersR.begin(); clusterIter != dnaClustersR.end(); clusterIter++)
	{
		IntegerVecMapIter overlapsIter = mDnaRnaOverlap.find(clusterIter->first);
		if (overlapsIter != mDnaRnaOverlap.end() && overlapsIter->second.size() > 0)
		{
			numColumns++;
		}
	}
	for (IntegerVecMapIter clusterIter = rnaClustersR.begin(); clusterIter != rnaClustersR.end(); clusterIter++)
	{
		IntegerVecMapIter overlapsIter = mRnaDnaOverlap.find(clusterIter->first);
		if (overlapsIter != mRnaDnaOverlap.end() && overlapsIter->second.size() > 0)
		{
			numColumns++;
		}
	}
	numColumns += Count(dnaClustersR);
	numColumns += Count(rnaClustersR);
	
	// Add columns
	glp_add_cols(lp, numColumns);
	
	// Column index
	int column = 1;
	
	// Dna columns
	IntegerMap dnaCluster2QColumn;
	IntegerMap dnaCluster2RColumn;
	for (IntegerVecMapIter clusterIter = dnaClustersR.begin(); clusterIter != dnaClustersR.end(); clusterIter++)
	{
		stringstream sstr;
		
		sstr.str("");
		sstr << "dna_cluster_indicator_" << clusterIter->first;
		
		glp_set_col_name(lp, column, sstr.str().c_str());
		glp_set_col_bnds(lp, column, GLP_DB, 0.0, 1.0);
		glp_set_col_kind(lp, column, GLP_BV);
		glp_set_obj_coef(lp, column, mW2);
		dnaCluster2QColumn[clusterIter->first] = column;
		column++;
		
		IntegerVecMapIter overlapsIter = mDnaRnaOverlap.find(clusterIter->first);
		if (overlapsIter == mDnaRnaOverlap.end() || overlapsIter->second.size() == 0)
		{
			continue;
		}
		
		sstr.str("");
		sstr << "dna_cluster_supported_" << clusterIter->first;
		
		glp_set_col_name(lp, column, sstr.str().c_str());
		glp_set_col_bnds(lp, column, GLP_DB, 0.0, 1.0);
		glp_set_col_kind(lp, column, GLP_BV);
		glp_set_obj_coef(lp, column, mW1 - mW2);
		dnaCluster2RColumn[clusterIter->first] = column;
		column++;
	}
	
	// Rna Columns
	IntegerMap rnaCluster2SColumn;
	IntegerMap rnaCluster2TColumn;
	for (IntegerVecMapIter clusterIter = rnaClustersR.begin(); clusterIter != rnaClustersR.end(); clusterIter++)
	{
		stringstream sstr;
		
		sstr.str("");
		sstr << "rna_cluster_indicator_" << clusterIter->first;
		
		glp_set_col_name(lp, column, sstr.str().c_str());
		glp_set_col_bnds(lp, column, GLP_DB, 0.0, 1.0);
		glp_set_col_kind(lp, column, GLP_BV);
		glp_set_obj_coef(lp, column, mW4);
		rnaCluster2SColumn[clusterIter->first] = column;
		column++;
		
		IntegerVecMapIter overlapsIter = mRnaDnaOverlap.find(clusterIter->first);
		if (overlapsIter == mRnaDnaOverlap.end() || overlapsIter->second.size() == 0)
		{
			continue;
		}
		
		sstr.str("");
		sstr << "rna_cluster_supported_" << clusterIter->first;
		
		glp_set_col_name(lp, column, sstr.str().c_str());
		glp_set_col_bnds(lp, column, GLP_DB, 0.0, 1.0);
		glp_set_col_kind(lp, column, GLP_BV);
		glp_set_obj_coef(lp, column, mW3 - mW4);
		rnaCluster2TColumn[clusterIter->first] = column;
		column++;
	}
	
	// Dna assignment columns
	IntegerPairMap dnaAssignmentColumn;
	for (IntegerVecMapIter clusterIter = dnaClustersR.begin(); clusterIter != dnaClustersR.end(); clusterIter++)
	{
		for (IntegerVecIter fragmentIter = clusterIter->second.begin(); fragmentIter != clusterIter->second.end(); fragmentIter++)
		{
			stringstream sstr;
			
			sstr.str("");
			sstr << "dna_cluster_fragment_assignment_" << clusterIter->first << "-" << *fragmentIter;
			
			glp_set_col_name(lp, column, sstr.str().c_str());
			glp_set_col_bnds(lp, column, GLP_DB, 0.0, 1.0);
			glp_set_col_kind(lp, column, GLP_BV);
			glp_set_obj_coef(lp, column, 0.0);
			dnaAssignmentColumn[IntegerPair(clusterIter->first,*fragmentIter)] = column;
			column++;
		}
	}
	
	// Rna assignment columns
	IntegerPairMap rnaAssignmentColumn;
	for (IntegerVecMapIter clusterIter = rnaClustersR.begin(); clusterIter != rnaClustersR.end(); clusterIter++)
	{
		for (IntegerVecIter fragmentIter = clusterIter->second.begin(); fragmentIter != clusterIter->second.end(); fragmentIter++)
		{
			stringstream sstr;
			
			sstr.str("");
			sstr << "rna_cluster_fragment_assignment_" << clusterIter->first << "-" << *fragmentIter;
			
			glp_set_col_name(lp, column, sstr.str().c_str());
			glp_set_col_bnds(lp, column, GLP_DB, 0.0, 1.0);
			glp_set_col_kind(lp, column, GLP_BV);
			glp_set_obj_coef(lp, column, 0.0);
			rnaAssignmentColumn[IntegerPair(clusterIter->first,*fragmentIter)] = column;
			column++;
		}
	}
	
	DebugCheck(column - 1 == numColumns);
	
	// Count number of rows required
	int numRows = dnaFragments.size() + rnaFragments.size();
	for (IntegerVecMapIter clusterIter = dnaClustersR.begin(); clusterIter != dnaClustersR.end(); clusterIter++)
	{
		IntegerVecMapIter overlapsIter = mDnaRnaOverlap.find(clusterIter->first);
		if (overlapsIter != mDnaRnaOverlap.end() && overlapsIter->second.size() > 0)
		{
			numRows += 2;
		}
	}
	for (IntegerVecMapIter clusterIter = rnaClustersR.begin(); clusterIter != rnaClustersR.end(); clusterIter++)
	{
		IntegerVecMapIter overlapsIter = mRnaDnaOverlap.find(clusterIter->first);
		if (overlapsIter != mRnaDnaOverlap.end() && overlapsIter->second.size() > 0)
		{
			numRows += 2;
		}
	}
	numRows += Count(dnaClustersR);
	numRows += Count(rnaClustersR);
	
	// Add rows
	glp_add_rows(lp, numRows);
	
	// Row index
	int row = 1;
	
	// Add dna ifelse rows
	for (IntegerVecMapIter clusterIter = dnaClustersR.begin(); clusterIter != dnaClustersR.end(); clusterIter++)
	{
		IntegerVecMapIter overlapsIter = mDnaRnaOverlap.find(clusterIter->first);
		if (overlapsIter == mDnaRnaOverlap.end() || overlapsIter->second.size() == 0)
		{
			continue;
		}
		
		stringstream sstr;
		sstr << "dna_cluster_indicated_supported_" << clusterIter->first;
		
		glp_set_row_name(lp, row, sstr.str().c_str());
		glp_set_row_bnds(lp, row, GLP_LO, 0.0, 0.0);
		row++;
	}
	
	// Add rna ifelse rows
	for (IntegerVecMapIter clusterIter = rnaClustersR.begin(); clusterIter != rnaClustersR.end(); clusterIter++)
	{
		IntegerVecMapIter overlapsIter = mRnaDnaOverlap.find(clusterIter->first);
		if (overlapsIter == mRnaDnaOverlap.end() || overlapsIter->second.size() == 0)
		{
			continue;
		}
		
		stringstream sstr;
		sstr << "rna_cluster_indicated_supported_" << clusterIter->first;
		
		glp_set_row_name(lp, row, sstr.str().c_str());
		glp_set_row_bnds(lp, row, GLP_LO, 0.0, 0.0);
		row++;
	}
	
	// Add dna overlap rows
	for (IntegerVecMapIter clusterIter = dnaClustersR.begin(); clusterIter != dnaClustersR.end(); clusterIter++)
	{
		IntegerVecMapIter overlapsIter = mDnaRnaOverlap.find(clusterIter->first);
		if (overlapsIter == mDnaRnaOverlap.end() || overlapsIter->second.size() == 0)
		{
			continue;
		}
		
		stringstream sstr;
		sstr << "dna_cluster_supported_by_rna_clusters_" << clusterIter->first;
		
		glp_set_row_name(lp, row, sstr.str().c_str());
		glp_set_row_bnds(lp, row, GLP_LO, 0.0, 0.0);
		row++;
	}
	
	// Add rna overlap rows
	for (IntegerVecMapIter clusterIter = rnaClustersR.begin(); clusterIter != rnaClustersR.end(); clusterIter++)
	{
		IntegerVecMapIter overlapsIter = mRnaDnaOverlap.find(clusterIter->first);
		if (overlapsIter == mRnaDnaOverlap.end() || overlapsIter->second.size() == 0)
		{
			continue;
		}
		
		stringstream sstr;
		sstr << "rna_cluster_supported_by_dna_clusters_" << clusterIter->first;
		
		glp_set_row_name(lp, row, sstr.str().c_str());
		glp_set_row_bnds(lp, row, GLP_LO, 0.0, 0.0);
		row++;
	}
	
	// Add dna fragment assignment rows
	for (IntegerVecMapIter fragmentIter = dnaFragments.begin(); fragmentIter != dnaFragments.end(); fragmentIter++)
	{
		stringstream sstr;
		sstr << "dna_fragment_assigned_once_" << fragmentIter->first;
		
		glp_set_row_name(lp, row, sstr.str().c_str());
		glp_set_row_bnds(lp, row, GLP_FX, 1.0, 1.0);
		row++;
	}
	
	// Add rna fragment assignment rows
	for (IntegerVecMapIter fragmentIter = rnaFragments.begin(); fragmentIter != rnaFragments.end(); fragmentIter++)
	{
		stringstream sstr;
		sstr << "rna_fragment_assigned_once_" << fragmentIter->first;
		
		glp_set_row_name(lp, row, sstr.str().c_str());
		glp_set_row_bnds(lp, row, GLP_FX, 1.0, 1.0);
		row++;
	}
	
	// Add dna cluster support rows
	for (IntegerVecMapIter clusterIter = dnaClustersR.begin(); clusterIter != dnaClustersR.end(); clusterIter++)
	{
		for (IntegerVecIter fragmentIter = clusterIter->second.begin(); fragmentIter != clusterIter->second.end(); fragmentIter++)
		{
			stringstream sstr;
			sstr << "dna_cluster_supported_by_fragment_" << clusterIter->first << "-" << *fragmentIter;
			
			glp_set_row_name(lp, row, sstr.str().c_str());
			glp_set_row_bnds(lp, row, GLP_LO, 0.0, 0.0);
			row++;
		}
	}
	
	// Add rna cluster support rows
	for (IntegerVecMapIter clusterIter = rnaClustersR.begin(); clusterIter != rnaClustersR.end(); clusterIter++)
	{
		for (IntegerVecIter fragmentIter = clusterIter->second.begin(); fragmentIter != clusterIter->second.end(); fragmentIter++)
		{
			stringstream sstr;
			sstr << "rna_cluster_supported_by_fragment_" << clusterIter->first << "-" << *fragmentIter;
		
			glp_set_row_name(lp, row, sstr.str().c_str());
			glp_set_row_bnds(lp, row, GLP_LO, 0.0, 0.0);
			row++;
		}
	}
		
	DebugCheck(row - 1 == numRows);
	
	// Create matrix
	IntegerVec ia;
	IntegerVec ja;
	DoubleVec ar;
	ia.push_back(0);
	ja.push_back(0);
	ar.push_back(0.0);
	
	// Row index
	row = 1;
	
	// Add dna ifelse constraints
	for (IntegerVecMapIter clusterIter = dnaClustersR.begin(); clusterIter != dnaClustersR.end(); clusterIter++)
	{
		IntegerVecMapIter overlapsIter = mDnaRnaOverlap.find(clusterIter->first);
		if (overlapsIter == mDnaRnaOverlap.end() || overlapsIter->second.size() == 0)
		{
			continue;
		}
		
		ia.push_back(row);
		ja.push_back(dnaCluster2QColumn[clusterIter->first]);
		ar.push_back(1.0);
		
		ia.push_back(row);
		ja.push_back(dnaCluster2RColumn[clusterIter->first]);
		ar.push_back(-1.0);
		row++;
	}
	
	// Add rna ifelse constraints
	for (IntegerVecMapIter clusterIter = rnaClustersR.begin(); clusterIter != rnaClustersR.end(); clusterIter++)
	{
		IntegerVecMapIter overlapsIter = mRnaDnaOverlap.find(clusterIter->first);
		if (overlapsIter == mRnaDnaOverlap.end() || overlapsIter->second.size() == 0)
		{
			continue;
		}
		
		ia.push_back(row);
		ja.push_back(rnaCluster2SColumn[clusterIter->first]);
		ar.push_back(1.0);
		
		ia.push_back(row);
		ja.push_back(rnaCluster2TColumn[clusterIter->first]);
		ar.push_back(-1.0);
		row++;
	}
	
	// Add dna overlap constraints
	for (IntegerVecMapIter clusterIter = dnaClustersR.begin(); clusterIter != dnaClustersR.end(); clusterIter++)
	{
		IntegerVecMapIter overlapsIter = mDnaRnaOverlap.find(clusterIter->first);
		if (overlapsIter == mDnaRnaOverlap.end() || overlapsIter->second.size() == 0)
		{
			continue;
		}
		
		ia.push_back(row);
		ja.push_back(dnaCluster2RColumn[clusterIter->first]);
		ar.push_back(-1.0);
		
		for (IntegerVecConstIter overlapIter = overlapsIter->second.begin(); overlapIter != overlapsIter->second.end(); overlapIter++)
		{
			ia.push_back(row);
			ja.push_back(rnaCluster2SColumn[*overlapIter]);
			ar.push_back(1.0);
		}
		
		row++;
	}
	
	// Add rna overlap constraints
	for (IntegerVecMapIter clusterIter = rnaClustersR.begin(); clusterIter != rnaClustersR.end(); clusterIter++)
	{
		IntegerVecMapIter overlapsIter = mRnaDnaOverlap.find(clusterIter->first);
		if (overlapsIter == mRnaDnaOverlap.end() || overlapsIter->second.size() == 0)
		{
			continue;
		}
		
		ia.push_back(row);
		ja.push_back(rnaCluster2TColumn[clusterIter->first]);
		ar.push_back(-1.0);
		
		for (IntegerVecConstIter overlapIter = overlapsIter->second.begin(); overlapIter != overlapsIter->second.end(); overlapIter++)
		{
			ia.push_back(row);
			ja.push_back(dnaCluster2QColumn[*overlapIter]);
			ar.push_back(1.0);
		}
		
		row++;
	}
	
	// Add dna fragment assignment constraints
	for (IntegerVecMapIter fragmentIter = dnaFragments.begin(); fragmentIter != dnaFragments.end(); fragmentIter++)
	{
		for (IntegerVecIter clusterIter = fragmentIter->second.begin(); clusterIter != fragmentIter->second.end(); clusterIter++)
		{
			ia.push_back(row);
			ja.push_back(dnaAssignmentColumn[IntegerPair(*clusterIter,fragmentIter->first)]);
			ar.push_back(1.0);
		}
		
		row++;
	}
	
	// Add rna fragment assignment constraints
	for (IntegerVecMapIter fragmentIter = rnaFragments.begin(); fragmentIter != rnaFragments.end(); fragmentIter++)
	{
		for (IntegerVecIter clusterIter = fragmentIter->second.begin(); clusterIter != fragmentIter->second.end(); clusterIter++)
		{
			ia.push_back(row);
			ja.push_back(rnaAssignmentColumn[IntegerPair(*clusterIter,fragmentIter->first)]);
			ar.push_back(1.0);
		}
		
		row++;
	}
	
	// Add dna cluster support constraints
	for (IntegerVecMapIter clusterIter = dnaClustersR.begin(); clusterIter != dnaClustersR.end(); clusterIter++)
	{
		for (IntegerVecIter fragmentIter = clusterIter->second.begin(); fragmentIter != clusterIter->second.end(); fragmentIter++)
		{
			ia.push_back(row);
			ja.push_back(dnaCluster2QColumn[clusterIter->first]);
			ar.push_back(1.0);
			
			ia.push_back(row);
			ja.push_back(dnaAssignmentColumn[IntegerPair(clusterIter->first,*fragmentIter)]);
			ar.push_back(-1.0);
			
			row++;
		}
	}
	
	// Add rna cluster support constraints
	for (IntegerVecMapIter clusterIter = rnaClustersR.begin(); clusterIter != rnaClustersR.end(); clusterIter++)
	{
		for (IntegerVecIter fragmentIter = clusterIter->second.begin(); fragmentIter != clusterIter->second.end(); fragmentIter++)
		{
			ia.push_back(row);
			ja.push_back(rnaCluster2SColumn[clusterIter->first]);
			ar.push_back(1.0);
			
			ia.push_back(row);
			ja.push_back(rnaAssignmentColumn[IntegerPair(clusterIter->first,*fragmentIter)]);
			ar.push_back(-1.0);
			
			row++;
		}
	}
	
	glp_load_matrix(lp, ia.size() - 1, &ia.front(), &ja.front(), &ar.front());
	
	glp_smcp smcp;
	glp_init_smcp(&smcp);

	glp_iocp iocp;
	glp_init_iocp(&iocp);

	if (ia.size() - 1 < 100000)
	{
//		smcp.msg_lev = GLP_MSG_OFF;
//		iocp.msg_lev = GLP_MSG_OFF;
	}
	
	iocp.tm_lim = mTimeout;
	
	glp_simplex(lp, &smcp);
	glp_intopt(lp, &iocp);

	glp_write_prob(lp, 0, "prob.txt");
	glp_print_sol(lp, "sol.txt");
	
	for (IntegerVecMapIter clusterIter = dnaClustersR.begin(); clusterIter != dnaClustersR.end(); clusterIter++)
	{
		for (IntegerVecIter fragmentIter = clusterIter->second.begin(); fragmentIter != clusterIter->second.end(); fragmentIter++)
		{
			double value = glp_mip_col_val(lp, dnaAssignmentColumn[IntegerPair(clusterIter->first,*fragmentIter)]);
		
			DebugCheck(value == 0.0 || value == 1.0);

			if (value > 0.0)
			{
				dnaSolution[clusterIter->first].push_back(*fragmentIter);
			}
		}
	}
	
	for (IntegerVecMapIter clusterIter = rnaClustersR.begin(); clusterIter != rnaClustersR.end(); clusterIter++)
	{
		for (IntegerVecIter fragmentIter = clusterIter->second.begin(); fragmentIter != clusterIter->second.end(); fragmentIter++)
		{
			double value = glp_mip_col_val(lp, rnaAssignmentColumn[IntegerPair(clusterIter->first,*fragmentIter)]);
			
			DebugCheck(value == 0.0 || value == 1.0);
			
			if (value > 0.0)
			{
				rnaSolution[clusterIter->first].push_back(*fragmentIter);
			}
		}
	}
	
	return glp_mip_obj_val(lp);
}
	

