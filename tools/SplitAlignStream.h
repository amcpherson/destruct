/*
 *  SplitAlignStream.h
 *
 */

#ifndef SPLITALIGNSTREAM_H_
#define SPLITALIGNSTREAM_H_

#include "Common.h"

#include <iostream>
#include <fstream>

using namespace std;


struct RawPosition
{
	string reference;
	int strand;
	int position;
};

struct RawSplitAlignment
{
	string fragment;
	int readEnd;
	RawPosition positions[2];
	string inserted;
	string line;
};


class SplitAlignmentStream
{
public:
	SplitAlignmentStream(const string& alignFilename);
	~SplitAlignmentStream();
	
	bool GetNextAlignment(RawSplitAlignment& alignment);
	
protected:
	istream* mStream;
	istream* mFileStream;
	int mLineNumber;
};


#endif
