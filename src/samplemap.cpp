/*
 *  samplegc.cpp
 *
 */

#include "DebugCheck.h"
#include "Sequences.h"
#include "api/BamReader.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;
using namespace BamTools;


double CalculateMappability(const Sequences& mappability, const string& chromosome, long start, int fragmentLength, int alignLength)
{
	vector<uint8_t> fragmentMappability;
	mappability.Get(chromosome, start, start + fragmentLength - alignLength - 1, fragmentMappability);
	
	double mappabilitySum = 0.0;
	for (int i = 0; i < fragmentMappability.size(); i++)
	{
		if (fragmentMappability[i] == 0)
		{
			return -1.0;
		}
		
		mappabilitySum += 1.0 / (double)fragmentMappability[i];
	}
	
	return mappabilitySum / (double)fragmentMappability.size();
}

int main(int argc, char* argv[])
{
	int numSamples;
	int fragmentLength;
	int alignLength;
	string mappabilityFilename;
	string bamFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Sample Positions and Calculate Mappability");
		TCLAP::ValueArg<int> numSamplesArg("n","num","Number of Samples",true,0,"integer",cmd);
		TCLAP::ValueArg<int> fragmentLengthArg("f","frlen","Mean Fragment Length",true,0,"integer",cmd);
		TCLAP::ValueArg<int> alignLengthArg("a","alignlen","Aligned Length for Mappability",true,0,"integer",cmd);
		TCLAP::ValueArg<string> mappabilityFilenameArg("m","map","Mappability BedGraph Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> bamFilenameArg("b","bam","Bam Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		numSamples = numSamplesArg.getValue();
		fragmentLength = fragmentLengthArg.getValue();
		alignLength = alignLengthArg.getValue();
		mappabilityFilename = mappabilityFilenameArg.getValue();
		bamFilename = bamFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	BamReader bamReader;
	if (!bamReader.Open(bamFilename))
	{
		cerr << "Error: Unable to open bam file " << bamFilename << endl;
		exit(1);
	}
	
	vector<long> chromosomeLengths;
	for (int chrIdx = 0; chrIdx < bamReader.GetReferenceCount(); chrIdx++)
	{
		chromosomeLengths.push_back(bamReader.GetReferenceData()[chrIdx].RefLength);
	}
	
	cerr << "Reading Mappability" << endl;
	
	Sequences mappability;
	mappability.ReadMappabilityBedGraph(mappabilityFilename);
	
	cerr << "Sampling Positions" << endl;
	
	RandomGenomicPositionGenerator randomPosition(chromosomeLengths);
	
	typedef pair<int,long> SamplePosition;
	typedef pair<double,int> MappabilityCount;
	
	unordered_map<SamplePosition,MappabilityCount> samples;
	
	while (samples.size() < numSamples)
	{
		// Sample 1 based position
		int chrIdx;
		long position;
		randomPosition.Next(chrIdx, position);
		
		const string& chromosome = bamReader.GetReferenceData()[chrIdx].RefName;
		
		double positionMappability = CalculateMappability(mappability, chromosome, position, fragmentLength, alignLength);
		
		if (positionMappability < 0.0)
		{
			continue;
		}
		
		samples[SamplePosition(chrIdx, position)] = MappabilityCount(positionMappability, 0);
	}
	
	cerr << "Counting reads from bam" << endl;
	
	int rowCount = 0;
	
	BamAlignment alignment;
	while (bamReader.GetNextAlignmentCore(alignment))
	{
		if (++rowCount % 2000000 == 0)
		{
			cerr << ".";
			cerr.flush();
		}
		
		if (alignment.IsProperPair() && !alignment.IsReverseStrand() && alignment.MapQuality > 0 && alignment.InsertSize < 2*fragmentLength)
		{
			unordered_map<SamplePosition,MappabilityCount>::iterator sampleIter = samples.find(SamplePosition(alignment.RefID, alignment.Position + 1));
			if (sampleIter != samples.end())
			{
				sampleIter->second.second++;
			}
		}
	}
	
	for (unordered_map<SamplePosition,MappabilityCount>::const_iterator sampleIter = samples.begin(); sampleIter != samples.end(); sampleIter++)
	{
		cout << bamReader.GetReferenceData()[sampleIter->first.first].RefName << "\t";
		cout << sampleIter->first.second << "\t";
		cout << sampleIter->second.first << "\t";
		cout << sampleIter->second.second << endl;
	}
}

