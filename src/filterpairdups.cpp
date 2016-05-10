/*
 *  filterpairdups.cpp
 *
 */

#include "DebugCheck.h"
#include "ReadStream.h"

#include <fstream>
#include <tclap/CmdLine.h>

using namespace boost;
using namespace std;


int main(int argc, char* argv[])
{
	string inReads1Filename;
	string inReads2Filename;
	string outReads1Filename;
	string outReads2Filename;
	
	try
	{
		TCLAP::CmdLine cmd("Filter read pairs dups");
		TCLAP::ValueArg<string> inReads1FilenameArg("","in1","Input End 1 Fastq Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> inReads2FilenameArg("","in2","Input End 2 Fastq Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> outReads1FilenameArg("","out1","Output End 1 Fastq Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> outReads2FilenameArg("","out2","Output End 2 Fastq Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		inReads1Filename = inReads1FilenameArg.getValue();
		inReads2Filename = inReads2FilenameArg.getValue();
		outReads1Filename = outReads1FilenameArg.getValue();
		outReads2Filename = outReads2FilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	ifstream inReads1File(inReads1Filename.c_str());
	CheckFile(inReads1File, inReads1Filename);
	FastqReadStream inReads1Stream(inReads1File);
	
	ifstream inReads2File(inReads2Filename.c_str());
	CheckFile(inReads2File, inReads2Filename);
	FastqReadStream inReads2Stream(inReads2File);

	ofstream outReads1File(outReads1Filename.c_str());
	CheckFile(outReads1File, outReads1Filename);
	
	ofstream outReads2File(outReads2Filename.c_str());
	CheckFile(outReads2File, outReads2Filename);
	
	RawRead read1;
	RawRead read2;
	while (inReads1Stream.GetNextRead(read1) && inReads2Stream.GetNextRead(read2))
	{
		if (read1.sequence == read2.sequence)
		{
			outReads1File << read1.fragment << "/" << read1.readEnd << endl;
			outReads1File << "+" << endl;
			outReads1File << read1.sequence << endl;
			outReads1File << read1.quality << endl;
			
			outReads2File << read2.fragment << "/" << read2.readEnd << endl;
			outReads2File << "+" << endl;
			outReads2File << read2.sequence << endl;
			outReads2File << read2.quality << endl;
		} 
	}
}

