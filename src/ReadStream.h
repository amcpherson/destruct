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
#include <boost/shared_ptr.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace std;
using namespace boost;

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
	explicit FastqReadStream(const string& filename);
	bool Good();
	bool GetNextRead(RawRead& read);
private:
	boost::shared_ptr<ifstream> mFile;
	boost::shared_ptr<iostreams::filtering_streambuf<iostreams::input> > mStreamBuf;
	boost::shared_ptr<istream> mStream;
};

#endif
