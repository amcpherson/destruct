/*
 *  SplitAlignStream.cpp
 *
 */

#include "SplitAlignStream.h"

#include <fstream>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


SplitAlignmentStream::SplitAlignmentStream(const string& alignFilename) : mStream(0), mFileStream(0), mLineNumber(0)
{
	if (alignFilename == "-")
	{
		mStream = &cin;
	}
	else
	{
		mFileStream = new ifstream(alignFilename.c_str());
		mStream = mFileStream;
		
		if (!mStream->good())
		{
			cerr << "Error: Unable to open alignment file " << alignFilename << endl;
			exit(1);
		}
	}
}

SplitAlignmentStream::~SplitAlignmentStream()
{
	delete mFileStream;
}

bool SplitAlignmentStream::GetNextAlignment(RawSplitAlignment& alignment)
{
	string line;
	while (getline(*mStream, line))
	{
		mLineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty alignment line " << mLineNumber << endl;
			exit(1);
		}
		
		vector<string> alignmentFields;
		split(alignmentFields, line, is_any_of("\t"));
		
		if (alignmentFields.size() < 6)
		{
			cerr << "Error: Format error for alignment line " << mLineNumber << endl;
			exit(1);
		}
		
		alignment.fragment = alignmentFields[0];
		alignment.readEnd = (alignmentFields[1] == "0") ? 0 : 1;
		alignment.positions[0].reference = alignmentFields[2];
		alignment.positions[0].strand = (alignmentFields[3] == "-") ? MinusStrand : PlusStrand;
		alignment.positions[0].position = SAFEPARSE(int, alignmentFields[4]);
		alignment.positions[1].reference = alignmentFields[5];
		alignment.positions[1].strand = (alignmentFields[6] == "-") ? MinusStrand : PlusStrand;
		alignment.positions[1].position = SAFEPARSE(int, alignmentFields[7]);
		alignment.inserted = alignmentFields[8];
		alignment.line = line;
		
		return true;
	}
	
	return false;	
}


