/*
 *  qualtrimfastq
 *
 *  Created by Andrew McPherson on 14/10/12.
 *
 */

#include "Common.h"
#include "DebugCheck.h"

#include <fstream>
#include <iostream>
#include <string>
#include <tclap/CmdLine.h>

using namespace boost;
using namespace std;


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

int GetTrimmedLength(const string& qualities, int offset, int threshold)
{
	int trimmedLength = (int)qualities.length();
	while (trimmedLength > 0 && (int)qualities[trimmedLength-1] < offset + threshold)
	{
		trimmedLength--;
	}
	return trimmedLength;
}

int main(int argc, char* argv[])
{
	string inSeqs1Filename;
	string inSeqs2Filename;
	string outSeqs1Filename;
	string outSeqs2Filename;
	int qualityOffset;
	int qualityThreshold;
	int minLength;
	
	try
	{
		TCLAP::CmdLine cmd("Trim low quality ends of read sequences");
		TCLAP::UnlabeledValueArg<string> inSeqs1FilenameArg("in1","Input End 1 Sequences",true,"","filename",cmd);
		TCLAP::UnlabeledValueArg<string> inSeqs2FilenameArg("in2","Input End 2 Sequences",true,"","filename",cmd);
		TCLAP::UnlabeledValueArg<string> outSeqs1FilenameArg("out1","Output End 1 Sequences",true,"","filename",cmd);
		TCLAP::UnlabeledValueArg<string> outSeqs2FilenameArg("out2","Output End 2 Sequences",true,"","filename",cmd);
		TCLAP::ValueArg<int> qualityOffsetArg("o","offset","Quality Offset, Sanger:33, Illumina:64, Illumina1.8+:33",true,-1,"integer",cmd);
		TCLAP::ValueArg<int> qualityThresholdArg("q","qualmin","Quality Threshold",true,-1,"integer",cmd);
		TCLAP::ValueArg<int> minLengthArg("l","lenmin","Minimum trimmed length",true,-1,"integer",cmd);
		cmd.parse(argc,argv);
		
		inSeqs1Filename = inSeqs1FilenameArg.getValue();
		inSeqs2Filename = inSeqs2FilenameArg.getValue();
		outSeqs1Filename = outSeqs1FilenameArg.getValue();
		outSeqs2Filename = outSeqs2FilenameArg.getValue();
		qualityOffset = qualityOffsetArg.getValue();
		qualityThreshold = qualityThresholdArg.getValue();
		minLength = minLengthArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	ifstream inSeqs1File(inSeqs1Filename.c_str());
	ifstream inSeqs2File(inSeqs2Filename.c_str());
	
	CheckFile(inSeqs1File, inSeqs1Filename);
	CheckFile(inSeqs2File, inSeqs2Filename);
	
	ofstream outSeqs1File(outSeqs1Filename.c_str());
	ofstream outSeqs2File(outSeqs2Filename.c_str());
	
	CheckFile(outSeqs1File, outSeqs1Filename);
	CheckFile(outSeqs2File, outSeqs2Filename);
	
	FastqReadStream inSeqs1Stream(inSeqs1File);
	FastqReadStream inSeqs2Stream(inSeqs2File);
	
	RawRead rawRead1;
	RawRead rawRead2;
	while (inSeqs1Stream.Next(rawRead1) && inSeqs2Stream.Next(rawRead2))
	{
		int trimmedLength1 = GetTrimmedLength(rawRead1.qualities(), qualityOffset, qualityThreshold);
		int trimmedLength2 = GetTrimmedLength(rawRead2.qualities(), qualityOffset, qualityThreshold);
		
		if (trimmedLength1 < minLength || trimmedLength2 < minLength)
		{
			continue;
		}
		
		rawRead1.trim(trimmedLength1);
		rawRead2.trim(trimmedLength2);
		
		outSeqs1File << rawRead1;
		outSeqs2File << rawRead2;
	}
}



