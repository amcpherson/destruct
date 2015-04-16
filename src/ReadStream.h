/*
 *  ReadStream.h
 *  findbreaks
 *
 *  Created by Andrew McPherson on 12/10/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef READSTREAM_H_
#define READSTREAM_H_

#include "Common.h"

#include <iostream>
#include <fstream>

using namespace std;

struct RawRead
{
	string fragment;
	int readEnd;
	string sequence;
	string quality;
};

class IReadStream
{
public:	
	virtual ~IReadStream() {};
	
	virtual bool Good() = 0;
	virtual bool GetNextRead(RawRead& read) = 0;
};

class FastqReadStream : public IReadStream
{
public:
	explicit FastqReadStream(istream& in);	
	bool Good();
	bool GetNextRead(RawRead& read);
private:
	istream& mStream;
};

#endif
