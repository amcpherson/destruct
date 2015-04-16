/*
 *  RegionDB.cpp
 *
 *  Created by Andrew McPherson on 12-01-12.
 *
 */

#include "RegionDB.h"
#include "Common.h"
#include "DebugCheck.h"

#include <fstream>
#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;

	
void RegionDB::Add(const CompactSimpleRegion& region)
{
	for (int bin = region.start / mBinSpacing; bin <= region.end / mBinSpacing; bin++)
	{
		mBinned[RefBin(region.ref,bin)].push_back(StartEnd(region.start,region.end));
	}
}

void RegionDB::Add(const CompactSimpleRegionVec& regions)
{
	for (CompactSimpleRegionVec::const_iterator regionIter = regions.begin(); regionIter != regions.end(); regionIter++)
	{
		Add(*regionIter);
	}
}

void RegionDB::Add(const string& regionsFilename)
{
	// Open regions file
	ifstream regionsFile(regionsFilename.c_str());
	if (!regionsFile)
	{
		cerr << "Error: unable to read from regions file " << regionsFilename << endl;		
		exit(1);
	}
	
	// Parse file contents
	string line;
	int lineNumber = 0;
	while (getline(regionsFile, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty break regions line " << lineNumber << endl;
			exit(1);
		}
		
		vector<string> regionFields;
		split(regionFields, line, is_any_of("\t"));
		
		if (regionFields.size() < 3)
		{
			cerr << "Error: Format error for regions line " << lineNumber << endl;
			exit(1);
		}
		
		CompactSimpleRegion region;
		region.ref = regionFields[0];
		region.start = SAFEPARSE(int, regionFields[1]);
		region.end = SAFEPARSE(int, regionFields[2]);
		
		Add(region);
	}
	
	regionsFile.close();
}	

bool RegionDB::Contained(string ref, int start, int end) const
{
	StartEnd queryStartEnd(start,end);
	
	for (int bin = start / mBinSpacing; bin <= end / mBinSpacing; bin++)
	{
		unordered_map<RefBin,StartEndVec>::const_iterator findIter = mBinned.find(RefBin(ref,bin));
		
		if (findIter != mBinned.end())
		{
			for (StartEndVec::const_iterator regionIter = findIter->second.begin(); regionIter != findIter->second.end(); regionIter++)
			{
				if (Contains(*regionIter,queryStartEnd))
				{
					return true;
				}
			}
		}			
	}
	
	return false;
}

bool RegionDB::Overlapped(string ref, int start, int end) const
{
	StartEnd queryStartEnd(start,end);
	
	for (int bin = start / mBinSpacing; bin <= end / mBinSpacing; bin++)
	{
		unordered_map<RefBin,StartEndVec>::const_iterator findIter = mBinned.find(RefBin(ref,bin));
		
		if (findIter != mBinned.end())
		{
			for (StartEndVec::const_iterator regionIter = findIter->second.begin(); regionIter != findIter->second.end(); regionIter++)
			{
				if (Overlaps(*regionIter,queryStartEnd))
				{
					return true;
				}
			}
		}			
	}
	
	return false;
}

void RegionDB::Clear()
{
	mBinned.clear();
}

