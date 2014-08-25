/*
 *  dosplitalign.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "SplitAlignmentGen.h"
#include "Common.h"
#include "DebugCheck.h"
#include "ReadStream.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>

using namespace boost;
using namespace std;


int main(int argc, char* argv[])
{
	string alignFilename;
	string referenceFasta;
	string alignmentsFilename;
	double fragmentLengthMean;
	double fragmentLengthStdDev;
	string splitAlignmentsFilename;
	string readSeqs1Filename;
	string readSeqs2Filename;

	try
	{
		TCLAP::CmdLine cmd("Fusion sequence prediction by split reads");
		TCLAP::ValueArg<string> alignFilenameArg("i","input","Input Alignment Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("f","fasta","Reference Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> fragmentLengthMeanArg("u","ufrag","Fragment Length Mean",true,0.0,"float",cmd);
		TCLAP::ValueArg<double> fragmentLengthStdDevArg("s","sfrag","Fragment Length Standard Deviation",true,0.0,"float",cmd);
		TCLAP::ValueArg<string> splitAlignmentsFilenameArg("o","out","Output Split Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqs1FilenameArg("1","seq1","End 1 Sequences",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqs2FilenameArg("2","seq2","End 2 Sequences",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		alignFilename = alignFilenameArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		alignmentsFilename = alignmentsFilenameArg.getValue();
		fragmentLengthMean = fragmentLengthMeanArg.getValue();
		fragmentLengthStdDev = fragmentLengthStdDevArg.getValue();
		splitAlignmentsFilename = splitAlignmentsFilenameArg.getValue();
		readSeqs1Filename = readSeqs1FilenameArg.getValue();
		readSeqs2Filename = readSeqs2FilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	cerr << "Reading input alignment regions" << endl;
	
	LocationVecMap alignRegionPairs;
	ReadAlignRegionPairs(alignFilename, alignRegionPairs);
	
	cerr << "Reading input reference fasta" << endl;
	
	Sequences reference;
	reference.Read(referenceFasta);
	
	cerr << "Initializing split alignment" << endl;
	
	SplitAlignment::SplitAlignmentMap splitAlignments;
	
	const int minReadLength = 50;
	const int maxReadLength = 100;
	
	for (LocationVecMapConstIter pairIter = alignRegionPairs.begin(); pairIter != alignRegionPairs.end(); pairIter++)
	{
		int id = pairIter->first;
		const LocationVec& alignRegionPair = pairIter->second;
		
		splitAlignments[id].Initialize(alignRegionPair, fragmentLengthMean, fragmentLengthStdDev, minReadLength, maxReadLength);
	}
	
	cerr << "Finding candidate split reads" << endl;
	
	AlignmentStream* alignmentStream = new CompactAlignmentStream(alignmentsFilename);
	
	SplitAlignment::FindCandidates(alignmentStream, splitAlignments);
	
	cerr << "Reading candidate read sequences" << endl;
	
	ifstream readSeqs1File(readSeqs1Filename.c_str());
	ifstream readSeqs2File(readSeqs2Filename.c_str());
	
	CheckFile(readSeqs1File, readSeqs1Filename);
	CheckFile(readSeqs2File, readSeqs2Filename);
	
	FastqReadStream readSeqs1Stream(readSeqs1File);
	FastqReadStream readSeqs2Stream(readSeqs2File);
	
	StringVec readSeqs;
	SplitAlignment::ReadCandidateSequences(&readSeqs1Stream, splitAlignments, readSeqs);
	SplitAlignment::ReadCandidateSequences(&readSeqs2Stream, splitAlignments, readSeqs);
	
	cerr << "Split alignment" << endl;
	
	for (SplitAlignment::SplitAlignmentMapIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		SplitAlignment& splitAlignment = splitAlignIter->second;		
		splitAlignment.Align(reference, readSeqs);
	}
	
	ofstream splitAlignmentsFile(splitAlignmentsFilename.c_str());
	CheckFile(splitAlignmentsFile, splitAlignmentsFilename);
	
	SplitAlignment::WriteAlignments(splitAlignmentsFile, splitAlignments);
}



