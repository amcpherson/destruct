/*
 *  filterreads.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "DebugCheck.h"
#include "AlignmentRecord.h"
#include "Indexer.h"
#include "RegionDB.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


bool IsFiltered(const vector<SpanningAlignmentRecord>& alignments, const RegionDB& excludedRegions, int numEnds)
{
	int excluded[] = {0,0};
	for (vector<SpanningAlignmentRecord>::const_iterator alignmentIter = alignments.begin(); alignmentIter != alignments.end(); alignmentIter++)
	{
		if (excludedRegions.Overlapped(alignmentIter->chromosome, alignmentIter->position, alignmentIter->position))
		{
			excluded[alignmentIter->readEnd] = 1;
		}
	}
	
	if (excluded[0] + excluded[1] >= numEnds)
	{
		return true;
	}
	
	return false;
}

int main(int argc, char* argv[])
{
	string alignmentsFilename;
	string regionsFilename;
	int numEnds;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Pair Filtering Tool");
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Read Sorted Compact Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> regionsFilenameArg("r","regions","Excluded Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> numEndsArg("n","numend","Number of Ends for Exclusion",true,0,"integer",cmd);
		cmd.parse(argc,argv);
		
		alignmentsFilename = alignmentsFilenameArg.getValue();
		regionsFilename = regionsFilenameArg.getValue();
		numEnds = numEndsArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	RegionDB excludedRegions;
	excludedRegions.Add(regionsFilename);
	
	ifstream alignmentsFile(alignmentsFilename.c_str());
	CheckFile(alignmentsFile, alignmentsFilename);

	GroupedRecordsStream<SpanningAlignmentRecord> alignmentsStream(alignmentsFile);

	vector<SpanningAlignmentRecord> alignments;
	while (alignmentsStream.Next(alignments, ReadEqual<SpanningAlignmentRecord>))
	{
		if (IsFiltered(alignments, excludedRegions, numEnds))
		{
			continue;
		}
		
		cout << alignments;
	}
}

