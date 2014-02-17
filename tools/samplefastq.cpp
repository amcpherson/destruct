/*
 *  samplefastq.cpp
 *
 */

#include "Common.h"
#include "DebugCheck.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>

using namespace boost;
using namespace std;


bool IsInSample(RandomNumberGenerator& rng, int previousReadCount, int numSamples, int& sampleIndex)
{
	if (previousReadCount < numSamples)
	{
		sampleIndex = previousReadCount;
		
		return true;
	}
	
	sampleIndex = rng.Next(0, previousReadCount - 1);
	
	if (sampleIndex < numSamples)
	{
		return true;
	}
	
	return false;
}

struct RawRead
{
	const string& readname()
	{
		return lines[0];
	}
	
	const string& sequence()
	{
		return lines[1];
	}
	
	const string& comment()
	{
		return lines[2];
	}
	
	const string& qualities()
	{
		return lines[3];
	}
	
	void trim(int length)
	{
		lines[1] = lines[1].substr(0, length);
		lines[3] = lines[3].substr(0, length);
	}
	
	string lines[4];
};

std::ostream& operator<<(std::ostream& stream, const RawRead& rawRead)
{
	stream << rawRead.lines[0] << endl;
	stream << rawRead.lines[1] << endl;
	stream << rawRead.lines[2] << endl;
	stream << rawRead.lines[3] << endl;
	return stream;
}

class FastqReadStream
{
public:
	FastqReadStream(istream& in) : mStream(in)
	{
	}
	
	bool Good()
	{
		return mStream.good();
	}
	
	bool Next(RawRead& read)
	{
		int lineIndex = 0;
		while (lineIndex < 4 && getline(mStream, read.lines[lineIndex]))
		{
			lineIndex++;
		}
		
		if (lineIndex < 4)
		{
			return false;
		}
		
		return true;
	}
	
private:
	istream& mStream;
};

int main(int argc, char* argv[])
{
	string inSeqs1Filename;
	string inSeqs2Filename;
	string outSeqs1Filename;
	string outSeqs2Filename;
	int numSamples;
	
	try
	{
		TCLAP::CmdLine cmd("Sample fastq tool");
		TCLAP::UnlabeledValueArg<string> inSeqs1FilenameArg("in1","Input End 1 Sequences",true,"","filename",cmd);
		TCLAP::UnlabeledValueArg<string> inSeqs2FilenameArg("in2","Input End 2 Sequences",true,"","filename",cmd);
		TCLAP::UnlabeledValueArg<string> outSeqs1FilenameArg("out1","Output End 1 Sequences",true,"","filename",cmd);
		TCLAP::UnlabeledValueArg<string> outSeqs2FilenameArg("out2","Output End 2 Sequences",true,"","filename",cmd);
		TCLAP::ValueArg<int> numSamplesArg("n","num","Number of Samples",true,0,"integer",cmd);
		cmd.parse(argc,argv);
		
		inSeqs1Filename = inSeqs1FilenameArg.getValue();
		inSeqs2Filename = inSeqs2FilenameArg.getValue();
		outSeqs1Filename = outSeqs1FilenameArg.getValue();
		outSeqs2Filename = outSeqs2FilenameArg.getValue();
		numSamples = numSamplesArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	RandomNumberGenerator rng;
	
	ifstream inSeqs1File(inSeqs1Filename.c_str());
	ifstream inSeqs2File(inSeqs2Filename.c_str());
	
	CheckFile(inSeqs1File, inSeqs1Filename);
	CheckFile(inSeqs2File, inSeqs2Filename);
	
	FastqReadStream inSeqs1Stream(inSeqs1File);
	FastqReadStream inSeqs2Stream(inSeqs2File);
	
	int readCount = 0;
	vector<RawRead> samples1(numSamples);
	vector<RawRead> samples2(numSamples);
	
	RawRead rawRead1;
	RawRead rawRead2;
	while (inSeqs1Stream.Next(rawRead1) && inSeqs2Stream.Next(rawRead2))
	{
		int sampleIndex;
		if (IsInSample(rng, readCount, numSamples, sampleIndex))
		{
			samples1[sampleIndex] = rawRead1;
			samples2[sampleIndex] = rawRead2;
		}
		
		readCount++;
	}
	
	if (readCount < numSamples)
	{
		cerr << "Error: more samples requested than reads" << endl;
		exit(1);
	}
	
	ofstream outSeqs1File(outSeqs1Filename.c_str());
	ofstream outSeqs2File(outSeqs2Filename.c_str());
	
	CheckFile(outSeqs1File, outSeqs1Filename);
	CheckFile(outSeqs2File, outSeqs2Filename);
	
	for (int sampleIndex = 0; sampleIndex < samples1.size(); sampleIndex++)
	{
		outSeqs1File << samples1[sampleIndex];
		outSeqs2File << samples2[sampleIndex];
	}
}

