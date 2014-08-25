/*
 *  samplemappable.cpp
 *
 */

#include "DebugCheck.h"
#include "Sequences.h"
#include "RegionDB.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;


void ReadMappabilityWig(const string& mappabilityWigFilename, unordered_map<string,vector<uint8_t> >& mappability)
{
	ifstream mappabilityWigFile(mappabilityWigFilename.c_str());
	CheckFile(mappabilityWigFile, mappabilityWigFilename);
	
	StringVec fields;
	int line = 1;
	while (ReadTSV(mappabilityWigFile, fields))
	{
		if (fields.size() < 4)
		{
			cerr << "Error: line " << line << " has too few fields" << endl;
			exit(1);
		}
		
		string chromosome = fields[0];
		int start = SAFEPARSE(int, fields[1]) + 1;
		int end = SAFEPARSE(int, fields[2]);
		double mapRegion = SAFEPARSE(float, fields[3]);
		uint8_t mapRegionApprox = (uint8_t)(mapRegion * 255.5);
		
		if (chromosome.substr(0, 3) == "chr")
		{
			chromosome = chromosome.substr(3);
		}
		
		if (mappability[chromosome].size() < end + 1)
		{
			mappability[chromosome].resize(end + 1, 0);
		}
		
		for (int pos = start; pos <= end; pos++)
		{
			mappability[chromosome][pos] = mapRegionApprox;
		}
		
		line++;
	}
}

double CalculateMappability(const unordered_map<string,vector<uint8_t> >& mappability, const string& chromosome, long position, int neighbourhood)
{
	unordered_map<string,vector<uint8_t> >::const_iterator mappIter = mappability.find(chromosome);
	
	if (mappIter == mappability.end())
	{
		return 0.0;
	}
	
	double mappabilitySum = 0.0;
	for (long pos = position - neighbourhood; pos <= position + neighbourhood; pos++)
	{
		if (position >= 0 && position < mappIter->second.size())
		{
			mappabilitySum += (double)mappIter->second[position];
		}
	}
	
	return mappabilitySum / (255.0 * (double)(2 * neighbourhood + 1));
}

int main(int argc, char* argv[])
{
	int numSamples;
	string faiFilename;
	string mappabilityFilename;
	string regionsFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Sample Positions and Calculate Mappability");
		TCLAP::ValueArg<int> numSamplesArg("n","num","Number of Samples",true,0,"integer",cmd);
		TCLAP::ValueArg<string> faiFilenameArg("f","fai","Samtools Fasta Index",false,"","string",cmd);
		TCLAP::ValueArg<string> mappabilityFilenameArg("m","map","Mappability BedGraph Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> regionsFilenameArg("r","regions","Excluded Regions Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		numSamples = numSamplesArg.getValue();
		faiFilename = faiFilenameArg.getValue();
		mappabilityFilename = mappabilityFilenameArg.getValue();
		regionsFilename = regionsFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	const int neighbourhood = 100;
	const double minMappability = 0.01;
	
	cerr << "Reading fasta index" << endl;
	
	vector<string> chromosomeNames;
	vector<long> chromosomeLengths;
	ReadFAI(faiFilename, chromosomeNames, chromosomeLengths);
	
	cerr << "Reading mappability" << endl;
	
	unordered_map<string,vector<uint8_t> > mappability;
	ReadMappabilityWig(mappabilityFilename, mappability);
	
	cerr << "Reading excluded regions" << endl;
	
	RegionDB excludedRegions;
	excludedRegions.Add(regionsFilename);
	
	cerr << "Sampling Positions" << endl;
	
	RandomGenomicPositionGenerator randomPosition(chromosomeLengths);
	
	int sampleCount = 0;
	while (sampleCount < numSamples)
	{
		int chrIdx;
		long position;
		randomPosition.Next(chrIdx, position);
		
		double positionMappability = CalculateMappability(mappability, chromosomeNames[chrIdx], position, neighbourhood);
		
		bool excluded = excludedRegions.Overlapped(chromosomeNames[chrIdx], position - neighbourhood, position + neighbourhood);
		
		if (excluded || positionMappability < minMappability)
		{
			continue;
		}
		
		cout << chromosomeNames[chrIdx] << "\t" << position << endl;
		
		sampleCount++;
	}
}

