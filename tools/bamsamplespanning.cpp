/*
 *  bamsamplespanning.cpp
 *
 *  Created by Andrew McPherson on 17/08/11.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Parsers.h"
#include "api/BamReader.h"

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

using namespace BamTools;


int main(int argc, char* argv[])
{
	string bamFilename;
	int numSamples;
	int maxFragmentLength;
	string ignoreList;
	int minPartialAlignLength;
	
	try
	{
		TCLAP::CmdLine cmd("Sample spanning distribution");
		TCLAP::ValueArg<string> bamFilenameArg("b","bam","Bam Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> numSamplesArg("n","numsamples","Number of Samples",true,0,"integer",cmd);
		TCLAP::ValueArg<int> maxFragmentLengthArg("m","maxfrag","Maximum Fragment Length",true,0,"integer",cmd);
		TCLAP::ValueArg<string> ignoreListArg("i","ignore","Reference Sequences to Ignore (comma separated)",false,"","string",cmd);
		TCLAP::ValueArg<int> minPartialAlignLengthArg("p","partial","Minimum Length of a Partial Alignment",true,0,"integer",cmd);
		cmd.parse(argc,argv);

		bamFilename = bamFilenameArg.getValue();
		numSamples = numSamplesArg.getValue();
		maxFragmentLength = maxFragmentLengthArg.getValue();
		ignoreList = ignoreListArg.getValue();
		minPartialAlignLength = minPartialAlignLengthArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	unordered_set<string> ignoreReferences;
	if (!ignoreList.empty())
	{
		vector<string> ignoreListFields;
		split(ignoreListFields, ignoreList, is_any_of(","));
		
		for (vector<string>::const_iterator fieldIter = ignoreListFields.begin(); fieldIter != ignoreListFields.end(); fieldIter++)
		{
			ignoreReferences.insert(*fieldIter);
		}
	}
	
	BamReader bamReader;
	if (!bamReader.Open(bamFilename))
	{
		cerr << "Error: Unable to open bam file " << bamFilename << endl;
		exit(1);
	}
	
	if (!bamReader.LocateIndex())
	{
		cerr << "Error: Unable to find index for bam file " << bamFilename << endl;
		exit(1);
	}
	
	vector<long> chromosomeLengths;
	for (int chrIdx = 0; chrIdx < bamReader.GetReferenceCount(); chrIdx++)
	{
		chromosomeLengths.push_back(bamReader.GetReferenceData()[chrIdx].RefLength);
	}
	
	RandomGenomicPositionGenerator randomPosition(chromosomeLengths);
	
	int samples = 0;
	while (samples < numSamples)
	{
		int chrIdx;
		long position;
		randomPosition.Next(chrIdx, position);
		
		const string& chromosome = bamReader.GetReferenceData()[chrIdx].RefName;
		
		if (ignoreReferences.find(chromosome) != ignoreReferences.end())
		{
			continue;
		}
		
		int spanningCount = 0;
		
		BamAlignment alignment;
		bamReader.SetRegion(chrIdx, position - maxFragmentLength, chrIdx, position + maxFragmentLength);
		while (bamReader.GetNextAlignmentCore(alignment))
		{
			if (!alignment.IsProperPair())
			{
				continue;
			}
			
			int leftPosition = alignment.Position + minPartialAlignLength - 1;
			int rightPosition = alignment.MatePosition + alignment.Length - minPartialAlignLength;
			
			if (leftPosition >= rightPosition)
			{
				continue;
			}
			
			if (leftPosition < position && rightPosition >= position)
			{
				spanningCount++;
			}
		}
		
		cout << spanningCount << endl;
		
		samples++;
	}
}

