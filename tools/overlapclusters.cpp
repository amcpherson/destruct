/*
 *  overlapclusters.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-08-11.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "CompactBreakRegion.h"
#include "DebugCheck.h"
#include "Matrix.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;

typedef unordered_map<IntegerPair,BrRegVec> RefPairBreakRegionsMap;
typedef unordered_map<IntegerPair,BrRegVec>::iterator RefPairBreakRegionsMapIter;

void CreateRefPairBreakRegionsMap(BrRegVec& breakRegions, RefPairBreakRegionsMap& refPairMap)
{
	while (breakRegions.size() != 0)
	{
		CompactBreakRegion breakRegionEnd1 = breakRegions.back(); breakRegions.pop_back();
		CompactBreakRegion breakRegionEnd2 = breakRegions.back(); breakRegions.pop_back();
		
		if (breakRegionEnd1.clustEnd.clusterID != breakRegionEnd2.clustEnd.clusterID)
		{
			cerr << "Error: break regions should be ordered by cluster ID" << endl;		
			exit(1);		
		}
		
		if (breakRegionEnd1.refStrand.referenceIndex > breakRegionEnd2.refStrand.referenceIndex)
		{
			swap(breakRegionEnd1,breakRegionEnd2);
		}
		
		refPairMap[IntegerPair(breakRegionEnd1.refStrand.referenceIndex,breakRegionEnd2.refStrand.referenceIndex)].push_back(breakRegionEnd1);
		refPairMap[IntegerPair(breakRegionEnd1.refStrand.referenceIndex,breakRegionEnd2.refStrand.referenceIndex)].push_back(breakRegionEnd2);
	}
}

int main(int argc, char* argv[])
{
	string dnaBreakRegionsFilename;
	string rnaBreakRegionsFilename;
	string overlapFilename;
	bool matchReciprocal;

	try
	{
		TCLAP::CmdLine cmd("Overlap between break regions tool");
		TCLAP::ValueArg<string> dnaBreakRegionsFilenameArg("a","dnabreaks","Input DNA Breakpoint Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> rnaBreakRegionsFilenameArg("b","rnabreaks","Input RNA Breakpoint Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> overlapFilenameArg("o","overlap","Output Overlap Filename",true,"","string",cmd);
		TCLAP::ValueArg<bool> matchReciprocalArg("r","reciprocal","Match Reciprocal Rearrangement",false,false,"boolean",cmd);
		cmd.parse(argc,argv);
		
		dnaBreakRegionsFilename = dnaBreakRegionsFilenameArg.getValue();
		rnaBreakRegionsFilename = rnaBreakRegionsFilenameArg.getValue();
		overlapFilename = overlapFilenameArg.getValue();
		matchReciprocal = matchReciprocalArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	if (sizeof(long) != 8)
	{
		cerr << "Error: sizeof(long) != 8" << endl;		
		exit(1);
	}
	
	cout << "Reading break regions" << endl;
	
	BrRegVec dnaBreakRegions;
	BrRegVec rnaBreakRegions;
	NameIndex referenceNames;

	ReadBreakRegions(dnaBreakRegionsFilename, referenceNames, dnaBreakRegions);
	ReadBreakRegions(rnaBreakRegionsFilename, referenceNames, rnaBreakRegions);
	
	if ((dnaBreakRegions.size() % 2) == 1)
	{
		cerr << "Error: dna break regions has uneven number of break regions" << endl;		
		exit(1);		
	}
	
	if ((rnaBreakRegions.size() % 2) == 1)
	{
		cerr << "Error: rna break regions has uneven number of break regions" << endl;		
		exit(1);		
	}
	
	cout << "Creating reference pair break regions map" << endl;

	RefPairBreakRegionsMap refPairDnaBreakRegionsMap;
	RefPairBreakRegionsMap refPairRnaBreakRegionsMap;
	
	CreateRefPairBreakRegionsMap(dnaBreakRegions, refPairDnaBreakRegionsMap);
	CreateRefPairBreakRegionsMap(rnaBreakRegions, refPairRnaBreakRegionsMap);
	
	// Open cluster pairs file
	ofstream overlapFile(overlapFilename.c_str());
	if (!overlapFile)
	{
		cerr << "Error: unable to write to overlap file" << endl;		
		exit(1);
	}
	
	cout << "Iterating reference pairs" << endl;
	
	// Iterate through break regions for each ref pair
	int numOverlaps = 0;
	for (RefPairBreakRegionsMapIter refPairIter = refPairDnaBreakRegionsMap.begin(); refPairIter != refPairDnaBreakRegionsMap.end(); refPairIter++)
	{
		const IntegerPair& refPair = refPairIter->first;
		if (refPairRnaBreakRegionsMap.find(refPair) == refPairRnaBreakRegionsMap.end())
		{
			continue;
		}
		
		BrRegVec& refPairDnaBreakRegions = refPairDnaBreakRegionsMap[refPair];
		BrRegVec& refPairRnaBreakRegions = refPairRnaBreakRegionsMap[refPair];

		sort(refPairDnaBreakRegions.begin(), refPairDnaBreakRegions.end(), BreakRegionStartLessThan);
		sort(refPairRnaBreakRegions.begin(), refPairRnaBreakRegions.end(), BreakRegionStartLessThan);
	
		typedef unordered_map<IntegerPair,IntegerPairVec> IntegerPairMapVec;
		typedef unordered_map<IntegerPair,IntegerPairVec>::const_iterator IntegerPairMapVecIter;		
		IntegerPairMapVec clusterPairs;

		BrRegVecConstIter dnaBrRegIter = refPairDnaBreakRegions.begin();
		BrRegVecConstIter rnaBrRegIter = refPairRnaBreakRegions.begin();
		while (dnaBrRegIter != refPairDnaBreakRegions.end())
		{
			while (rnaBrRegIter != refPairRnaBreakRegions.end() && (rnaBrRegIter->refStrand.referenceIndex < dnaBrRegIter->refStrand.referenceIndex || (rnaBrRegIter->refStrand.referenceIndex == dnaBrRegIter->refStrand.referenceIndex && rnaBrRegIter->end < dnaBrRegIter->start)))
			{
				rnaBrRegIter++;
			}

			BrRegVecConstIter brRegBOverlapIter = rnaBrRegIter;
			while (brRegBOverlapIter != refPairRnaBreakRegions.end() && (brRegBOverlapIter->refStrand.referenceIndex == dnaBrRegIter->refStrand.referenceIndex && brRegBOverlapIter->start <= dnaBrRegIter->end))
			{
				if (brRegBOverlapIter->end >= dnaBrRegIter->start)
				{
					IntegerPair clusterIDs(dnaBrRegIter->clustEnd.clusterID,brRegBOverlapIter->clustEnd.clusterID);
					IntegerPair clusterBreakIndices(dnaBrRegIter - refPairDnaBreakRegions.begin(), brRegBOverlapIter - refPairRnaBreakRegions.begin());
					
					clusterPairs[clusterIDs].push_back(clusterBreakIndices);
				}
				
				brRegBOverlapIter++;
			}
			
			dnaBrRegIter++;
		}
		
		for (IntegerPairMapVecIter clusterPairIter = clusterPairs.begin(); clusterPairIter != clusterPairs.end(); clusterPairIter++)
		{
			Matrix<int> found(2,2);
			found.Clear(0);
			
			Matrix<int> strandEqual(2,2);
			
			Matrix<Region> combined(2,2);

			unordered_set<int> referenceIndices;

			for (IntegerPairVecConstIter breakPairIter = clusterPairIter->second.begin(); breakPairIter != clusterPairIter->second.end(); breakPairIter++)
			{
				const CompactBreakRegion& dnaBreakRegion = refPairDnaBreakRegions[breakPairIter->first];
				const CompactBreakRegion& rnaBreakRegion = refPairRnaBreakRegions[breakPairIter->second];

				DebugCheck(dnaBreakRegion.refStrand.referenceIndex == rnaBreakRegion.refStrand.referenceIndex);
				
				referenceIndices.insert(dnaBreakRegion.refStrand.referenceIndex);
				
				int combinedStart;
				int combinedEnd;
				if (rnaBreakRegion.refStrand.strand == PlusStrand)
				{
					combinedStart = rnaBreakRegion.start;
					combinedEnd = dnaBreakRegion.end;
				}
				else
				{
					combinedStart = dnaBreakRegion.start;
					combinedEnd = rnaBreakRegion.end;
				}
				
				found(dnaBreakRegion.clustEnd.clusterEnd, rnaBreakRegion.clustEnd.clusterEnd) = 1;
				strandEqual(dnaBreakRegion.clustEnd.clusterEnd, rnaBreakRegion.clustEnd.clusterEnd) = (dnaBreakRegion.refStrand.strand == rnaBreakRegion.refStrand.strand);
				combined(dnaBreakRegion.clustEnd.clusterEnd, rnaBreakRegion.clustEnd.clusterEnd).start = combinedStart;
				combined(dnaBreakRegion.clustEnd.clusterEnd, rnaBreakRegion.clustEnd.clusterEnd).end = combinedEnd;
			}

			bool sameChromosome = (referenceIndices.size() == 1);
			
			typedef vector<Cell> CellVec;
			typedef vector<CellVec> CellTable;
			
			CellTable possibleCellPairs;
			possibleCellPairs.push_back(CellVec()); possibleCellPairs.back().push_back(Cell(0,0)); possibleCellPairs.back().push_back(Cell(1,1));
			possibleCellPairs.push_back(CellVec()); possibleCellPairs.back().push_back(Cell(0,1)); possibleCellPairs.back().push_back(Cell(1,0));
			
			bool foundOverlap = false;
			for (int possibilityIndex = 0; possibilityIndex < 2; possibilityIndex++)
			{
				const CellVec& possibleCellPair = possibleCellPairs[possibilityIndex];
				
				bool valid = found(possibleCellPair[0]) && found(possibleCellPair[1]);
				bool conflict = Overlap(combined(possibleCellPair[0]), combined(possibleCellPair[1])) && sameChromosome;
				bool normal = strandEqual(possibleCellPair[0]) && strandEqual(possibleCellPair[1]);
				bool reciprocal = !strandEqual(possibleCellPair[0]) && !strandEqual(possibleCellPair[1]);
				
				if (valid && !conflict && (normal || (matchReciprocal && reciprocal)))
				{
					foundOverlap = true;
				}
			}
			
			if (foundOverlap)
			{
				overlapFile << clusterPairIter->first.first << "\t" << clusterPairIter->first.second << endl;
				numOverlaps++;
			}
		}
	}
	
	overlapFile.close();

	cout << "Output " << numOverlaps << " overlaps" << endl;
}

