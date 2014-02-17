/*
 *  calcspandist.cpp
 *
 *  Created by Andrew McPherson on 17/08/11.
 *
 */

#include "AlignmentStream.h"
#include "Common.h"
#include "DebugCheck.h"
#include "Parsers.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/unordered_set.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;

class BinnedSamplePositions
{
public:
	BinnedSamplePositions(int binSpacing) : mBinSpacing(binSpacing) {}
	
	void Load( IntegerPairVec& samplePositions)
	{
		for (int samplePosIndex = 0; samplePosIndex < samplePositions.size(); samplePosIndex++)
		{
			int bin = samplePositions[samplePosIndex].second / mBinSpacing;
			mBinned[IntegerPair(samplePositions[samplePosIndex].first,bin)].push_back(samplePosIndex);
		}
	}
	
	void ApproxContained(int ref, int start, int end, unordered_set<int>& indices) const
	{
		int startBin = start / mBinSpacing;
		int endBin = end / mBinSpacing;
		
		for (int bin = startBin; bin <= endBin; bin++)
		{
			unordered_map<IntegerPair,IntegerVec>::const_iterator findIter = mBinned.find(IntegerPair(ref,bin));
			if (findIter != mBinned.end())
			{
				indices.insert(findIter->second.begin(),findIter->second.end());
			}
		}
	}
	
private:	
	int mBinSpacing;
	unordered_map<IntegerPair,IntegerVec> mBinned;
};

int main(int argc, char* argv[])
{
	string alignmentsFilename;
	int numSamples;
	int trimLength;
	string faiFilename;
	string ignoreList;
	
	try
	{
		TCLAP::CmdLine cmd("Calculate spanning distribution");
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Concordant Alignments Sam Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> numSamplesArg("n","numsamples","Number of Samples",true,0,"integer",cmd);
		TCLAP::ValueArg<int> trimLengthArg("t","trim","Trim Length for Spanning Alignments",true,0,"integer",cmd);		
		TCLAP::ValueArg<string> faiFilenameArg("f","fai","Samtools Fasta Index",false,"","string",cmd);
		TCLAP::ValueArg<string> ignoreListArg("i","ignore","Reference Sequences to Ignore (comma separated)",false,"","string",cmd);
		cmd.parse(argc,argv);

		alignmentsFilename = alignmentsFilenameArg.getValue();
		numSamples = numSamplesArg.getValue();
		trimLength = trimLengthArg.getValue();
		faiFilename = faiFilenameArg.getValue();
		ignoreList = ignoreListArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	vector<string> referenceNames;
	vector<long> referenceLengths;
	ReadFAI(faiFilename, referenceNames, referenceLengths);
	
	unordered_map<string,int> referenceNamesLookup;
	for (int refNameIndex = 0; refNameIndex < referenceNames.size(); refNameIndex++)
	{
		referenceNamesLookup[referenceNames[refNameIndex]] = refNameIndex;
	}
	
	vector<string> ignoreListFields;
	split(ignoreListFields, ignoreList, is_any_of(","));
	unordered_set<string> ignoreReferences(ignoreListFields.begin(),ignoreListFields.end());
	
	IntegerPairVec samplePositions;
	
	RandomGenomicPositionGenerator randomPosition(referenceLengths);
	
	int samples = 0;
	while (samples < numSamples)
	{
		int refNameIndex;
		long position;
		randomPosition.Next(refNameIndex, position);
		
		if (ignoreReferences.find(referenceNames[refNameIndex]) != ignoreReferences.end())
		{
			continue;
		}
		
		samplePositions.push_back(IntegerPair(refNameIndex,position));
		
		samples++;
	}
	
	IntegerVec sampleCounts(samplePositions.size(), 0);
	
	BinnedSamplePositions binnedSamplePositions(10000);
	binnedSamplePositions.Load(samplePositions);
	
	AlignmentStream* alignmentStream = new SamAlignmentStream(alignmentsFilename);
	FragmentAlignmentStream fragmentAlignments(alignmentStream);
	
	RawAlignmentVec alignments;
	while (fragmentAlignments.GetNextAlignments(alignments))
	{
		DebugCheck(alignments.size() > 0);
		
		if (alignments.size() != 2)
		{
			cerr << "Error: expected 2 alignments per fragment" << endl;
			cerr << "retrieved " << alignments.size() << " alignments for " << alignments[0].fragment << endl;		
			exit(1);
		}
		
		if (referenceNamesLookup.find(alignments[0].reference) == referenceNamesLookup.end())
		{
			continue;
		}
		
		int refNameIndex = referenceNamesLookup.find(alignments[0].reference)->second;
		int start = min(alignments[0].region.start + trimLength, alignments[1].region.start + trimLength);
		int end = max(alignments[0].region.end - trimLength, alignments[1].region.end - trimLength);
		
		unordered_set<int> containedSamplePositions;
		binnedSamplePositions.ApproxContained(refNameIndex, start, end, containedSamplePositions);
		
		for (unordered_set<int>::const_iterator samplePosIndexIter = containedSamplePositions.begin(); samplePosIndexIter != containedSamplePositions.end(); samplePosIndexIter++)
		{
			int samplePos = samplePositions[*samplePosIndexIter].second;
			if (samplePos >= start && samplePos <= end)
			{
				sampleCounts[*samplePosIndexIter]++;
			}
		}
	}
	
	for (int samplePosIndex = 0; samplePosIndex < sampleCounts.size(); samplePosIndex++)
	{
		cout << sampleCounts[samplePosIndex] << endl;
	}
}

