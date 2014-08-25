/*
 *  setcover_dnarna.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-08-24.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Algorithms.h"
#include "BinaryMinHeap.h"

#include <fstream>
#include <iostream>
#include <string>
#include <tclap/CmdLine.h>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;


void FindCommonRearrangements(const IntegerVecPairVec& corroboration, IntegerTable& rearrangementBreakpoints, IntegerTable& rearrangementFusionTranscripts, IntegerVecMap& dnaToRearrangement, IntegerVecMap& rnaToRearrangements)
{
	unordered_map<IntegerVec,int> uniqueRearrangements;
	for (IntegerVecPairVecConstIter corrobIter = corroboration.begin(); corrobIter != corroboration.end(); corrobIter++)
	{
		pair<unordered_map<IntegerVec,int>::iterator,bool> insertResult = uniqueRearrangements.insert(make_pair(corrobIter->second,rearrangementBreakpoints.size()));
		
		if (insertResult.second)
		{
			rearrangementBreakpoints.push_back(corrobIter->second);
			rearrangementFusionTranscripts.push_back(IntegerVec());
		}
		
		for (IntegerVecConstIter rnaClusterIter = corrobIter->first.begin(); rnaClusterIter != corrobIter->first.end(); rnaClusterIter++)
		{
			rearrangementFusionTranscripts[insertResult.first->second].push_back(*rnaClusterIter);
			rnaToRearrangements[*rnaClusterIter].push_back(insertResult.first->second);
		}
	}
	uniqueRearrangements.clear();
	
	for (int rearrangementIndex = 0; rearrangementIndex < rearrangementBreakpoints.size(); rearrangementIndex++)
	{
		for (IntegerVecConstIter dnaClusterIter = rearrangementBreakpoints[rearrangementIndex].begin(); dnaClusterIter != rearrangementBreakpoints[rearrangementIndex].end(); dnaClusterIter++)
		{
			dnaToRearrangement[*dnaClusterIter].push_back(rearrangementIndex);
		}
	}
}

void SetCoverRearrangments(const IntegerVecMap& sets, const IntegerVecMap& setToRearrangements, const IntegerTable& rearrangements, const DoubleVec& rearrangementBonus, IntegerVec& solution)
{
	BinaryMinHeap minHeap;
	
	IntegerMap rearrangementSizes;
	for (int rearrangementIndex = 0; rearrangementIndex < (int)rearrangements.size(); rearrangementIndex++)
	{
		rearrangementSizes[rearrangementIndex] = (int)rearrangements[rearrangementIndex].size();
	}
	
	IntegerMap setSizes;
	DoubleMap weights;
	IntegerVecMap elementSets;
	for (IntegerVecMapConstIter setIter = sets.begin(); setIter != sets.end(); setIter++)
	{
		setSizes[setIter->first] = (int)setIter->second.size();
		
		double bonus = 0.0;
		IntegerVecMapConstIter rearrangementsIter = setToRearrangements.find(setIter->first);
		if (rearrangementsIter != setToRearrangements.end())
		{
			for (IntegerVecConstIter rearrangementIter = rearrangementsIter->second.begin(); rearrangementIter != rearrangementsIter->second.end(); rearrangementIter++)
			{
				bonus += rearrangementBonus[*rearrangementIter] / (double)rearrangementSizes[*rearrangementIter];
			}
		}
		
		double weight = 1.0 - bonus;
		
		minHeap.Push(setIter->first, weight / (double)setIter->second.size());
		
		for (IntegerVecConstIter elementIter = setIter->second.begin(); elementIter != setIter->second.end(); elementIter++)
		{
			elementSets[*elementIter].push_back(setIter->first);
		}
	}
	
	IntegerSet assigned;
	IntegerSet zeroedRearrangements;
	while (!minHeap.Empty())
	{
		int nextSetID = minHeap.MinID();
		
		solution.push_back(nextSetID);
		
		const IntegerVec& set = sets.find(nextSetID)->second;
		
		IntegerSet alteredSets;
		IntegerSet rearrangementAlteredSets;
		for (IntegerVecConstIter elementIter = set.begin(); elementIter != set.end(); elementIter++)
		{
			if (assigned.insert(*elementIter).second)
			{
				for (IntegerVecConstIter setIter = elementSets[*elementIter].begin(); setIter != elementSets[*elementIter].end(); setIter++)
				{
					setSizes[*setIter]--;
					alteredSets.insert(*setIter);
					
					if (setSizes[*setIter] == 0)
					{
						rearrangementAlteredSets.insert(*setIter);
					}
				}
			}
		}
		
		rearrangementAlteredSets.insert(nextSetID);
		
		for (IntegerSetConstIter setIter = rearrangementAlteredSets.begin(); setIter != rearrangementAlteredSets.end(); setIter++)
		{
			IntegerVecMapConstIter rearrangementsIter = setToRearrangements.find(*setIter);
			if (rearrangementsIter != setToRearrangements.end())
			{
				for (IntegerVecConstIter rearrangementIter = rearrangementsIter->second.begin(); rearrangementIter != rearrangementsIter->second.end(); rearrangementIter++)
				{
					if (*setIter == nextSetID)
					{
						rearrangementSizes[*rearrangementIter]--;
					}
					else
					{
						zeroedRearrangements.insert(*rearrangementIter);
					}
					
					const IntegerVec& rearrangement = rearrangements[*rearrangementIter];
					
					for (IntegerVecConstIter rearrSetIter = rearrangement.begin(); rearrSetIter != rearrangement.end(); rearrSetIter++)
					{
						if (setSizes[*rearrSetIter] > 0)
						{
							alteredSets.insert(*rearrSetIter);
						}
					}
				}
			}
		}
		
		for (IntegerSetConstIter setIter = alteredSets.begin(); setIter != alteredSets.end(); setIter++)
		{
			DebugCheck(setSizes[*setIter] >= 0);
			
			if (setSizes[*setIter] > 0)
			{
				double bonus = 0.0;
				IntegerVecMapConstIter rearrangementsIter = setToRearrangements.find(*setIter);
				if (rearrangementsIter != setToRearrangements.end())
				{
					for (IntegerVecConstIter rearrangementIter = rearrangementsIter->second.begin(); rearrangementIter != rearrangementsIter->second.end(); rearrangementIter++)
					{
						if (zeroedRearrangements.find(*rearrangementIter) == zeroedRearrangements.end())
						{
							bonus += rearrangementBonus[*rearrangementIter] / (double)rearrangementSizes[*rearrangementIter];
						}
					}
				}
				
				double weight = 1.0 - bonus;
				
				minHeap.ReplaceKey(*setIter, weight / (double)setSizes[*setIter]);
			}
			else
			{
				minHeap.Remove(*setIter);
			}
		}
	}
}

void EstimateSpliceVariantCounts(const IntegerVecMap& rnaClusters, const IntegerTable& rearrangementFusionTranscripts, int multiplier, IntegerVec& spliceVariantCounts)
{
	spliceVariantCounts.resize(rearrangementFusionTranscripts.size());
	for (int rearrangementIndex = 0; rearrangementIndex < rearrangementFusionTranscripts.size(); rearrangementIndex++)
	{
		IntegerVecMap transcripts;
		DoubleMap weights;
		for (IntegerVecConstIter rnaClusterIter = rearrangementFusionTranscripts[rearrangementIndex].begin(); rnaClusterIter != rearrangementFusionTranscripts[rearrangementIndex].end(); rnaClusterIter++)
		{
			transcripts[*rnaClusterIter] = rnaClusters.find(*rnaClusterIter)->second;
			weights[*rnaClusterIter] = 1.0;
		}
		
		IntegerVec solutionTranscripts;
		SetCover(transcripts, weights, solutionTranscripts);
		
		spliceVariantCounts[rearrangementIndex] = (int)solutionTranscripts.size() * multiplier;
	}
}

void FlipFlopSetCover(const IntegerVecMap& dnaClusters, const IntegerVecMap& rnaClusters, const IntegerVecPairVec& corroboration, double epsilon, IntegerVec& dnaSolution, IntegerVec& rnaSolution)
{
	int initialSamples = 10;
	int numSamples = 10;
	
	IntegerVecMap dnaToRearrangements;
	IntegerVecMap rnaToRearrangements;
	IntegerTable rearrangementBreakpoints;
	IntegerTable rearrangementFusionTranscripts;
	FindCommonRearrangements(corroboration, rearrangementBreakpoints, rearrangementFusionTranscripts, dnaToRearrangements, rnaToRearrangements);
	
	IntegerVec spliceVariantCounts;
	EstimateSpliceVariantCounts(rnaClusters, rearrangementFusionTranscripts, initialSamples, spliceVariantCounts);
	
	IntegerVec rearrangementSuccess(rearrangementBreakpoints.size(), initialSamples);
	
	int iteration = initialSamples;
	while (iteration < initialSamples + numSamples)
	{
		// Calculate rearrangement weights
		DoubleVec rearrangementBonus(rearrangementBreakpoints.size(), 1.0);
		for (int rearrangementIndex = 0; rearrangementIndex < rearrangementFusionTranscripts.size(); rearrangementIndex++)
		{
			if (rearrangementFusionTranscripts[rearrangementIndex].size() != 0)
			{
				double spliceVariantCount = (double)spliceVariantCounts[rearrangementIndex] / (double)iteration;
				double pathLength = (double)rearrangementBreakpoints[rearrangementIndex].size();
				
				rearrangementBonus[rearrangementIndex] = epsilon * spliceVariantCount / (1.0 + pathLength);
			}
		}
		
		// Set cover for rearrangements
		dnaSolution.clear();
		SetCoverRearrangments(dnaClusters, dnaToRearrangements, rearrangementBreakpoints, rearrangementBonus, dnaSolution);
		
		// Update rearrangement success count
		IntegerSet dnaSolutionLookup(dnaSolution.begin(), dnaSolution.end());
		for (int rearrangementIndex = 0; rearrangementIndex < rearrangementBreakpoints.size(); rearrangementIndex++)
		{
			bool success = true;
			for (IntegerVecConstIter dnaClusterIter = rearrangementBreakpoints[rearrangementIndex].begin(); dnaClusterIter != rearrangementBreakpoints[rearrangementIndex].end(); dnaClusterIter++)
			{
				if (dnaSolutionLookup.find(*dnaClusterIter) == dnaSolutionLookup.end())
				{
					success = false;
					break;
				}
			}
			
			if (success)
			{
				rearrangementSuccess[rearrangementIndex]++;
			}
		}
		
		// Calculate fusion transcript weights
		DoubleMap rnaWeights;
		for (IntegerVecMapConstIter rnaClusterIter = rnaClusters.begin(); rnaClusterIter != rnaClusters.end(); rnaClusterIter++)
		{
			rnaWeights[rnaClusterIter->first] = 1.0;
		}
		for (int rearrangementIndex = 0; rearrangementIndex < rearrangementBreakpoints.size(); rearrangementIndex++)
		{
			double rearrangementIndicator = (double)rearrangementSuccess[rearrangementIndex] / ((double)iteration + 1.0);
			double pathLength = (double)rearrangementBreakpoints[rearrangementIndex].size();
			
			double bonus = epsilon * rearrangementIndicator / (1.0 + pathLength);
			
			for (IntegerVecConstIter rnaClusterIter = rearrangementFusionTranscripts[rearrangementIndex].begin(); rnaClusterIter != rearrangementFusionTranscripts[rearrangementIndex].end(); rnaClusterIter++)
			{
				rnaWeights[*rnaClusterIter] -= bonus;
			}
		}
		
		// Set cover for fusion transcripts
		rnaSolution.clear();
		SetCover(rnaClusters, rnaWeights, rnaSolution);
		
		// Update splice variant count estimate
		IntegerSet rnaSolutionLookup(rnaSolution.begin(), rnaSolution.end());
		for (int rearrangementIndex = 0; rearrangementIndex < rearrangementFusionTranscripts.size(); rearrangementIndex++)
		{
			for (IntegerVecConstIter rnaClusterIter = rearrangementFusionTranscripts[rearrangementIndex].begin(); rnaClusterIter != rearrangementFusionTranscripts[rearrangementIndex].end(); rnaClusterIter++)
			{
				if (rnaSolutionLookup.find(*rnaClusterIter) != rnaSolutionLookup.end())
				{
					spliceVariantCounts[rearrangementIndex]++;
				}
			}
		}
		
		iteration++;
	}
}

int main(int argc, char* argv[])
{
/*
	IntegerVecMap sets;
	IntegerVecMap setToRearrangements;
	IntegerTable rearrangements1;
	IntegerSet solution1;
	
	sets[10].push_back(1);
	sets[10].push_back(2);
	sets[10].push_back(3);
	
	sets[11].push_back(4);
	sets[11].push_back(5);
	sets[11].push_back(6);
	
	sets[12].push_back(1);
	sets[12].push_back(2);
	sets[12].push_back(3);
	
	sets[13].push_back(1);
	sets[13].push_back(2);
	sets[13].push_back(3);
	
	sets[14].push_back(4);
	sets[14].push_back(5);
	sets[14].push_back(6);
	
	sets[15].push_back(7);
	sets[15].push_back(8);
	sets[15].push_back(9);
	
	setToRearrangements[11].push_back(0);
	setToRearrangements[12].push_back(0);
	
	setToRearrangements[13].push_back(1);
	setToRearrangements[14].push_back(1);
	setToRearrangements[15].push_back(1);
	
	
	rearrangements1.push_back(IntegerVec());
	rearrangements1.back().push_back(11);
	rearrangements1.back().push_back(12);
	rearrangements1.push_back(IntegerVec());
	rearrangements1.back().push_back(13);
	rearrangements1.back().push_back(14);
	rearrangements1.back().push_back(15);
	
	SetCoverRearrangments(sets, setToRearrangements, rearrangements1, 0.001, solution1);
	
	exit(1);
	
*/
	string dnaClustersFilename;
	string rnaClustersFilename;
	string corroborationFilename;
	double epsilon;
	string dnaOutClustersFilename;
	string rnaOutClustersFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Set cover for rna/dna maximum parsimony");
		TCLAP::ValueArg<string> dnaClustersFilenameArg("d","dna","DNA clusters filename",true,"","string",cmd);
		TCLAP::ValueArg<string> rnaClustersFilenameArg("r","rna","RNA clusters filename",true,"","string",cmd);
		TCLAP::ValueArg<string> corroborationFilenameArg("c","corr","Corroboration filename",true,"","string",cmd);
		TCLAP::ValueArg<double> epsilonArg("e","epsilon","Bonus Weight for Corroboration",true,-1,"string",cmd);
		TCLAP::ValueArg<string> dnaOutClustersFilenameArg("a","outdna","Output DNA clusters filename",true,"","string",cmd);
		TCLAP::ValueArg<string> rnaOutClustersFilenameArg("b","outrna","Output RNA clusters filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		dnaClustersFilename = dnaClustersFilenameArg.getValue();
		rnaClustersFilename = rnaClustersFilenameArg.getValue();
		corroborationFilename = corroborationFilenameArg.getValue();
		epsilon = epsilonArg.getValue();
		dnaOutClustersFilename = dnaOutClustersFilenameArg.getValue();
		rnaOutClustersFilename = rnaOutClustersFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	cout << "Reading DNA clusters" << endl;
	
	ClusterMembership dnaMembership(dnaClustersFilename);
	
	IntegerVecMap dnaClusters;
	dnaMembership.Read(dnaClusters);
	
	cout << "Reading RNA clusters" << endl;
	
	ClusterMembership rnaMembership(rnaClustersFilename);
	
	IntegerVecMap rnaClusters;	
	rnaMembership.Read(rnaClusters);
	
	cout << "Reading corroboration" << endl;
	
	IntegerVecPairVec corroboration;
	ReadCorroboration(corroborationFilename, corroboration);
	
	cout << "Set cover for rearrangements" << endl;
	
	IntegerVec dnaSolution;
	IntegerVec rnaSolution;
	FlipFlopSetCover(dnaClusters, rnaClusters, corroboration, epsilon, dnaSolution, rnaSolution);
	
	cout << "Assigning fragments to solution sets" << endl;
	
	IntegerVecMap dnaAssignment;
	AssignInOrder(dnaClusters, dnaSolution, dnaAssignment);
	
	IntegerVecMap rnaAssignment;
	AssignInOrder(rnaClusters, rnaSolution, rnaAssignment);
	
	cout << "Writing out clusters" << endl;

	dnaMembership.Write(dnaOutClustersFilename, dnaAssignment);
	rnaMembership.Write(rnaOutClustersFilename, rnaAssignment);
}

