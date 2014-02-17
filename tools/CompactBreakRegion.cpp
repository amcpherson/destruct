/*
 *  CompactBreakRegion.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-08-11.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "CompactBreakRegion.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;

bool BreakRegionStartLessThan(const CompactBreakRegion& br1, const CompactBreakRegion& br2)
{
	if (br1.refStrand.referenceIndex < br2.refStrand.referenceIndex)
	{
		return true;
	}
	else if (br1.refStrand.referenceIndex > br2.refStrand.referenceIndex)
	{
		return false;
	}
	
	if (br1.start < br2.start)
	{
		return true;
	}
	
	return false;
}

void ReadBreakRegions(const string& breakRegionsFilename, NameIndex& referenceNames, BrRegVec& breakRegions)
{
	// Open breakpoint regions file
	ifstream breakRegionsFile(breakRegionsFilename.c_str());
	if (!breakRegionsFile)
	{
		cerr << "Error: unable to read from breakpoint regions file" << endl;		
		exit(1);
	}
	
	// Parse file contents
	string line;
	int lineNumber = 0;
	while (getline(breakRegionsFile, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty break regions line " << lineNumber << endl;
			exit(1);
		}
		
		vector<string> breakRegionFields;
		split(breakRegionFields, line, is_any_of("\t"));
		
		if (breakRegionFields.size() < 5)
		{
			cerr << "Error: Format error for break regions line " << lineNumber << endl;
			exit(1);
		}
		
		int clusterID = SAFEPARSE(int, breakRegionFields[0]);
		int clusterEnd = SAFEPARSE(int, breakRegionFields[1]);
		string referenceName = breakRegionFields[2];
		string strandName = breakRegionFields[3];
		int start = SAFEPARSE(int, breakRegionFields[4]);
		int end = SAFEPARSE(int, breakRegionFields[5]);
		
		// Interpret strand
		int strand;
		if (strandName == "+")
		{
			strand = PlusStrand;
		}
		else if (strandName == "-")
		{
			strand = MinusStrand;
		}
		else
		{
			cerr << "Error: Unable to interpret strand for break regions line " << lineNumber << endl;
			exit(1);
		}
		
		CompactBreakRegion breakRegion;
		breakRegion.clustEnd.clusterID = clusterID;
		breakRegion.clustEnd.clusterEnd = clusterEnd;
		breakRegion.refStrand.referenceIndex = referenceNames.Index(referenceName);
		breakRegion.refStrand.strand = strand;
		breakRegion.start = start;
		breakRegion.end = end;
		
		breakRegions.push_back(breakRegion);
	}
	
	breakRegionsFile.close();
}
