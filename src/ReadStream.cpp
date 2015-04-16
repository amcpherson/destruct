/*
 *  ReadStream.cpp
 *
 *  Created by Andrew McPherson on 12/10/09.
 *
 */

#include "ReadStream.h"

#include <fstream>
#include <boost/algorithm/string.hpp>

using namespace boost;

FastqReadStream::FastqReadStream(istream& in) : mStream(in)
{
}

bool FastqReadStream::Good()
{
	return mStream.good();
}

bool FastqReadStream::GetNextRead(RawRead& read)
{
	string line[4];
	int lineIndex = 0;
	
	while (lineIndex < 4 && getline(mStream, line[lineIndex]))
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

