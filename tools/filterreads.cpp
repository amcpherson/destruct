/*
 *  filterreads.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "DebugCheck.h"
#include "Indexer.h"
#include "AlignmentStream.h"
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


bool IsFiltered(const RawAlignmentVec& alignments, const RegionDB& excludedRegions, int numEnds)
{
	int excluded[] = {0,0};
	for (RawAlignmentVec::const_iterator alignmentIter = alignments.begin(); alignmentIter != alignments.end(); alignmentIter++)
	{
		if (excludedRegions.Overlapped(alignmentIter->reference, alignmentIter->region.start, alignmentIter->region.end))
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

void Write(ostream& out, const RawAlignmentVec& alignments)
{
	for (RawAlignmentVecConstIter alignmentIter = alignments.begin(); alignmentIter != alignments.end(); alignmentIter++)
	{
		out << alignmentIter->line << endl;
	}
}

int main(int argc, char* argv[])
{
	string samFilename;
	string alignFilename;
	string regionsFilename;
	int numEnds;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Pair Filtering Tool");
		TCLAP::ValueArg<string> samFilenameArg("s","sam","Read Sorted Sam Filename",false,"","string");
		TCLAP::ValueArg<string> alignFilenameArg("a","align","Read Sorted Compact Alignments Filename",true,"","string");
		TCLAP::ValueArg<string> regionsFilenameArg("r","regions","Excluded Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> numEndsArg("n","numend","Number of Ends for Exclusion",true,0,"integer",cmd);
		
		vector<TCLAP::Arg*> alignmentArgs;
		alignmentArgs.push_back(&samFilenameArg);
		alignmentArgs.push_back(&alignFilenameArg);
		cmd.xorAdd(alignmentArgs);
		
		cmd.parse(argc,argv);
		
		samFilename = samFilenameArg.getValue();
		alignFilename = alignFilenameArg.getValue();
		regionsFilename = regionsFilenameArg.getValue();
		numEnds = numEndsArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	AlignmentStream* alignmentStream = 0;
	if (!samFilename.empty())
	{
		alignmentStream = new SamAlignmentStream(samFilename);
	}
	else if (!alignFilename.empty())
	{
		alignmentStream = new CompactAlignmentStream(alignFilename);
	}
	
	FragmentAlignmentStream fragmentAlignmentStream(alignmentStream);
	
	RegionDB excludedRegions;
	excludedRegions.Add(regionsFilename);
	
	RawAlignmentVec alignments;
	while (fragmentAlignmentStream.GetNextAlignments(alignments))
	{
		if (IsFiltered(alignments, excludedRegions, numEnds))
		{
			continue;
		}
		
		Write(cout, alignments);
	}
}

