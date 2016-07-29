/*
 *  ReadStream.cpp
 *
 *  Created by Andrew McPherson on 12/10/09.
 *
 */

#include "ReadStream.h"

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace std;
using namespace boost;

FastqReadStream::FastqReadStream(const string& filename)
{
	if (boost::algorithm::ends_with(filename, ".gz"))
	{
		mFile = boost::shared_ptr<ifstream>(new ifstream(filename.c_str(), ios::in|ios::binary));
		CheckFile(*mFile, filename);
		mStreamBuf = boost::shared_ptr<iostreams::filtering_streambuf<iostreams::input> >(new iostreams::filtering_streambuf<iostreams::input>());
		mStreamBuf->push(iostreams::gzip_decompressor());
		mStreamBuf->push(*mFile);
		mStream = boost::shared_ptr<istream>(new istream(&*mStreamBuf));
	}
	else
	{
		mStream = boost::shared_ptr<istream>(new ifstream(filename.c_str()));
	}
}

bool FastqReadStream::Good()
{
	return mFile->good();
}

bool FastqReadStream::GetNextRead(RawRead& read)
{
	string line[4];
	int lineIndex = 0;

	while (lineIndex < 4 && getline(*mStream, line[lineIndex]))
	{
		lineIndex++;
	}
	
	if (lineIndex < 4)
	{
		return false;
	}
	
	if (line[0][0] != '@')
	{
		cerr << "Error: Unable to interpret read name " << line[0] << endl;
		exit(1);
	}
	
	int readEnd = 0;
	string::size_type readEndStart = line[0].find_first_of('/');
	if (readEndStart != std::string::npos && readEndStart + 1 < line[0].length())
	{
		if (line[0][readEndStart + 1] != '1' && line[0][readEndStart + 1] != '2')
		{
			cerr << "Error: Unable to interpret read end " << line[0] << endl;
			exit(1);
		}
		
		readEnd = (line[0][readEndStart + 1] == '1') ? 0 : 1;
	}
	
	read.fragment = line[0].substr(1, readEndStart - 1);
	read.readEnd = readEnd;
	read.sequence = line[1];
	read.quality = line[3];
	
	return true;
}

