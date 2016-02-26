/*
 *  bamreheader.cpp
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/SamSequenceDictionary.h"

#include <fstream>
#include <iostream>
#include <string>
#include <tclap/CmdLine.h>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;

using namespace BamTools;


int main(int argc, char* argv[])
{
	string fastaIndexFilename;
	string bamInputFilename;
	string bamOutputFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Bam Reheader Tool");
		TCLAP::ValueArg<string> fastaIndexFilenameArg("f","faidx","Fasta index filename",true,"","string",cmd);
		TCLAP::ValueArg<string> bamInputFilenameArg("i","inbam","Input bam filename",true,"","string",cmd);
		TCLAP::ValueArg<string> bamOutputFilenameArg("o","outbam","Output bam filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		fastaIndexFilename = fastaIndexFilenameArg.getValue();
		bamInputFilename = bamInputFilenameArg.getValue();
		bamOutputFilename = bamOutputFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}

	vector<string> referenceNames;
	vector<long> referenceLengths;
	ReadFAI(fastaIndexFilename, referenceNames, referenceLengths);

	SamHeader header;
	RefVector referenceData;

	unordered_map<string,int> refNameLookup;

	for (int idx = 0; idx < (int)referenceNames.size(); idx++)
	{
		header.Sequences.Add(referenceNames[idx], referenceLengths[idx]);
		referenceData.push_back(RefData(referenceNames[idx], referenceLengths[idx]));
		refNameLookup[referenceNames[idx]] = idx;
	}

	BamReader bamInput;
	if (!bamInput.Open(bamInputFilename))
	{
		cerr << "Error: Unable to read bam file " << bamInputFilename << endl;
		exit(1);
	}

	vector<int> refIDRemap(bamInput.GetReferenceCount(), -1);

	for (int idx = 0; idx < bamInput.GetReferenceCount(); idx++)
	{
		const string& refName = bamInput.GetReferenceData()[idx].RefName;

		unordered_map<string,int>::const_iterator refNameIter = refNameLookup.find(refName);

		if (refNameIter != refNameLookup.end())
		{
			refIDRemap[idx] = refNameIter->second;
		}
	}
	
	BamWriter bamOutput;
	if (!bamOutput.Open(bamOutputFilename, header, referenceData))
	{
		cerr << "Error: Unable to write to bam file " << bamOutputFilename << endl;
		exit(1);
	}

	BamAlignment alignment;
	while (bamInput.GetNextAlignmentCore(alignment))
	{
		if (alignment.RefID < 0 || alignment.MateRefID < 0)
		{
			continue;
		}

		alignment.RefID = refIDRemap[alignment.RefID];
		alignment.MateRefID = refIDRemap[alignment.MateRefID];

		if (alignment.RefID < 0 || alignment.MateRefID < 0)
		{
			continue;
		}
		
		bamOutput.SaveAlignment(alignment);
	}
}

