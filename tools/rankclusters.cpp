/*
 *  rankclusters.cpp
 *
 *  Created by Andrew McPherson
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Sequences.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <numeric>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


class MapWithDefault
{
public:
	MapWithDefault(double defaultValue, DoubleMap& map) : mDefaultValue(defaultValue), mMap(map) {}
	MapWithDefault(double defaultValue) : mDefaultValue(defaultValue) {}
	
	double operator [](int key) const
	{
		DoubleMapConstIter mapIter = mMap.find(key);
		if (mapIter == mMap.end())
		{
			return mDefaultValue;
		}
		else
		{
			return mapIter->second;
		}
	}
	
private:
	double mDefaultValue;
	DoubleMap mMap;
};

double CalculateAlignProbability(const DoubleTable& probabilities, const MapWithDefault& depthProbabilities)
{
	DoubleVec previousRow;
	previousRow.push_back(1.0 - probabilities[0][0] * probabilities[0][1]);
	previousRow.push_back(probabilities[0][0] * probabilities[0][1]);
	
	for (int l = 1; l < probabilities.size(); l++)
	{
		DoubleVec nextRow;
		
		double prob0 = previousRow[0] * (1.0 - probabilities[l][0] * probabilities[l][1]);
		nextRow.push_back(prob0);
		
		for (int k = 1; k <= l; k++)
		{
			double prob = previousRow[k-1] * probabilities[l][0] * probabilities[l][1] + previousRow[k] * (1.0 - probabilities[l][0] * probabilities[l][1]);
			nextRow.push_back(prob);
		}
		
		double prob = previousRow[l] * probabilities[l][0] * probabilities[l][1];
		nextRow.push_back(prob);
		
		swap(nextRow, previousRow);
	}
	
	double clusterProb = 0.0;
	for (int n = 1; n <= probabilities.size(); n++)
	{
		clusterProb += previousRow[n] * depthProbabilities[n];
	}
	
	return clusterProb;
}

double CalculateChimericProbability(const DoubleTable& probabilities)
{
	double sum = 0.0;
	
	for (int i = 0; i < probabilities.size(); i++)
	{
		sum += probabilities[i][0];
	}
	
	return sum / (double)probabilities.size();
}

double CalculateValidProbability(const DoubleTable& probabilities)
{
	double sum = 0.0;
	
	for (int i = 0; i < probabilities.size(); i++)
	{
		sum += probabilities[i][0] * probabilities[i][1];
	}
	
	return sum / (double)probabilities.size();
}

void OutputClusterProbabilities(const string& clustersFilename, const MapWithDefault& depthProbabilities)
{
	ifstream clustersFile(clustersFilename.c_str());
	CheckFile(clustersFile, clustersFilename);
	
	ClusterReader clusterReader(clustersFile);
	
	while (clusterReader.Next())
	{
		int clusterID = clusterReader.FetchClusterID();
		DoubleTable alignProbabilities = clusterReader.FetchAlignProbabilities();
		DoubleTable chimericProbabilities = clusterReader.FetchChimericProbabilities();
		DoubleTable validProbabilities = clusterReader.FetchValidProbabilities();
		
		double clusterProbability = CalculateAlignProbability(alignProbabilities, depthProbabilities);
		double chimericProbability = CalculateChimericProbability(chimericProbabilities);
		double validProbability = CalculateValidProbability(validProbabilities);
		
		cout << clusterID << "\t" << clusterProbability << "\t" << chimericProbability << "\t" << validProbability << endl;
	}
}

int main(int argc, char* argv[])
{
	string clustersFilename;
	string depthProbFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Cluster ranking tool");
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> depthProbFilenameArg("d","depthprob","Depth Probability Filename",false,"","string",cmd);
		cmd.parse(argc,argv);
		
		clustersFilename = clustersFilenameArg.getValue();
		depthProbFilename = depthProbFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	cerr << "Calculating cluster probabilities" << endl;
	
	if (depthProbFilename.empty())
	{
		OutputClusterProbabilities(clustersFilename, MapWithDefault(1.0));
	}
	else
	{
		DoubleMap depthProbability;
		ReadDoubleMap(depthProbFilename, depthProbability);
		OutputClusterProbabilities(clustersFilename, MapWithDefault(1.0,depthProbability));
	}
}


