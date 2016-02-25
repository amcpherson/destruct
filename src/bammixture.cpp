/*
 *  bammixture.cpp
 *
 *  Created by Andrew McPherson on 10/07/13.
 *
 */

#include "Common.h"
#include "DebugCheck.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <cstdlib>
#include <tclap/CmdLine.h>
#include <bamtools/api/BamReader.h>
#include <bamtools/api/BamWriter.h>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;

using namespace BamTools;


class PairedBamReader
{
public:
	explicit PairedBamReader(const string& bamFilename) : mRowCount(0)
	{
		if (!mBamReader.Open(bamFilename))
		{
			cerr << "Error: Unable to open bam file " << bamFilename << endl;
			exit(1);
		}
	}
	
	bool Next(BamAlignment& alignment1, BamAlignment& alignment2)
	{
		BamAlignment alignment;
		while (mBamReader.GetNextAlignment(alignment))
		{
			if (++mRowCount % 2000000 == 0)
			{
				cerr << ".";
				cerr.flush();
			}
			
			int readEnd = alignment.IsFirstMate() ? 0 : 1;
			int otherReadEnd = OtherReadEnd(readEnd);
			
			unordered_map<string,BamAlignment>::iterator otherEndIter = mReadBuffer[otherReadEnd].find(alignment.Name);
			
			if (otherEndIter != mReadBuffer[otherReadEnd].end())
			{
				if (alignment.IsFirstMate())
				{
					alignment1 = alignment;
					alignment2 = otherEndIter->second;
				}
				else
				{
					alignment1 = otherEndIter->second;
					alignment2 = alignment;
				}
				
				mReadBuffer[otherReadEnd].erase(otherEndIter);
				
				return true;
			}
			else
			{
				mReadBuffer[readEnd].insert(make_pair(alignment.Name, alignment));
			}
		}
		
		cerr << endl;
		
		return false;
	}
	
private:
	BamReader mBamReader;
	int mRowCount;
	unordered_map<string,BamAlignment> mReadBuffer[2];
};

int CountPairedReads(const string& bamFilename)
{
	BamReader bamReader;
	if (!bamReader.Open(bamFilename))
	{
		cerr << "Error: Unable to open bam file " << bamFilename << endl;
		exit(1);
	}
	
	BamAlignment alignment;
	int readCount = 0;
	while (bamReader.GetNextAlignmentCore(alignment))
	{
		readCount++;
	}
	
	return readCount / 2;
}

vector<int> SelectRandomSubset(int size, int subset)
{
	assert(size >= subset);
	
	vector<int> indices;
	for (int i = 0; i < size; i++)
	{
		indices.push_back(i);
	}
	
	random_shuffle(indices.begin(), indices.end());
	
	indices.resize(subset);
	
	sort(indices.rbegin(), indices.rend());
	
	return indices;
}

void SaveSelectedReads(BamWriter& bamWriter, const string& bamFilename, vector<int>& selected, const string& appendToName)
{
	PairedBamReader pairedBamReader(bamFilename);
	BamAlignment alignment1;
	BamAlignment alignment2;
	int readPairIndex = 0;
	while (pairedBamReader.Next(alignment1, alignment2))
	{
		if (selected.empty())
		{
			break;
		}
		
		assert(selected.back() >= readPairIndex);
		
		if (selected.back() == readPairIndex)
		{
			alignment1.Name += appendToName;
			alignment2.Name += appendToName;
			
			bamWriter.SaveAlignment(alignment1);
			bamWriter.SaveAlignment(alignment2);
			
			selected.pop_back();
		}
		
		readPairIndex++;
	}
}

int main(int argc, char* argv[])
{
	string bamFilenameA;
	string bamFilenameB;
	string bamOutFilename;
	int numReadsA;
	int numReadsB;
	
	try
	{
		TCLAP::CmdLine cmd("Bam Mixing");
		TCLAP::ValueArg<string> bamFilenameAArg("","bama","Bam A Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> bamFilenameBArg("","bamb","Bam B Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> numReadsAArg("","numa","Number of A Reads",true,0,"integer",cmd);
		TCLAP::ValueArg<int> numReadsBArg("","numb","Number of B Reads",true,0,"integer",cmd);
		TCLAP::ValueArg<string> bamOutFilenameArg("o","out","Bam Out Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		bamFilenameA = bamFilenameAArg.getValue();
		bamFilenameB = bamFilenameBArg.getValue();
		numReadsA = numReadsAArg.getValue();
		numReadsB = numReadsBArg.getValue();
		bamOutFilename = bamOutFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	srand(2013);
	
	cerr << "Counting reads from " << bamFilenameA << endl;
	
	int readCountA = CountPairedReads(bamFilenameA);
	
	cerr << "Bam file " << bamFilenameA << " has " << readCountA << " paired reads" << endl;
	
	cerr << "Counting reads from " << bamFilenameB << endl;
	
	int readCountB = CountPairedReads(bamFilenameB);
	
	cerr << "Bam file " << bamFilenameB << " has " << readCountB << " paired reads" << endl;
	
	if (numReadsA > readCountA || numReadsB > readCountB)
	{
		if (numReadsA > readCountA)
		{
			cerr << "Error: bam file " << bamFilenameA << " has only " << readCountA << " reads when " << numReadsA << " were requested" << endl;
		}
		
		if (numReadsB > readCountB)
		{
			cerr << "Error: bam file " << bamFilenameB << " has only " << readCountB << " reads when " << numReadsB << " were requested" << endl;
		}
		
		exit(1);
	}
	
	cerr << "Selecting reads for " << bamFilenameA << endl;
	
	vector<int> selectedA = SelectRandomSubset(readCountA, numReadsA);
	
	cerr << "Selecting reads for " << bamFilenameB << endl;
	
	vector<int> selectedB = SelectRandomSubset(readCountB, numReadsB);
	
	cerr << "Creating mixed bam" << endl;
	
	BamReader bamReaderA;
	if (!bamReaderA.Open(bamFilenameA))
	{
		cerr << "Error: Unable to open bam file " << bamFilenameA << endl;
		exit(1);
	}
	
	BamWriter bamOutWriter;
	if (!bamOutWriter.Open(bamOutFilename, bamReaderA.GetHeader(), bamReaderA.GetReferenceData()))
	{
		cerr << "Error: Unable to write to bam file " << bamOutFilename << endl;
		exit(1);
	}
	
	bamReaderA.Close();
	
	cerr << "Adding reads from " << bamFilenameA << endl;
	
	SaveSelectedReads(bamOutWriter, bamFilenameA, selectedA, "A");
	
	cerr << "Adding reads from " << bamFilenameB << endl;
	
	SaveSelectedReads(bamOutWriter, bamFilenameB, selectedB, "B");
	
	bamOutWriter.Close();
}

