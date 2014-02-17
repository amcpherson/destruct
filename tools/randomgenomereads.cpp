/*
 *  randomgenomereads.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Sequences.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>

using namespace boost;
using namespace std;


void GetRandomSequence(int length, string& sequence)
{
	char ntchars[] = {'A','C','T','G'};
	
	while (sequence.length() < length)
	{
		sequence = sequence + ntchars[rand() % 4];
	}
}

int main(int argc, char* argv[])
{
	string genomeFasta;
	int readLength;
	int numReads;
	
	try
	{
		TCLAP::CmdLine cmd("Generate random reads in fastq format");
		TCLAP::ValueArg<string> genomeFastaArg("g","genome","Genome Fasta",true,"","string",cmd);
		TCLAP::ValueArg<int> readLengthArg("l","length","Read Length",true,0,"int",cmd);
		TCLAP::ValueArg<int> numReadsArg("n","num","Number of Reads",true,0,"int",cmd);
		cmd.parse(argc,argv);
		
		genomeFasta = genomeFastaArg.getValue();
		readLength = readLengthArg.getValue();
		numReads = numReadsArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	cerr << "Read reference fasta" << endl;
	
	Sequences genomeSequences;
	genomeSequences.Read(genomeFasta);
	
	for (int i = 0; i < numReads; i++)
	{
		string sequence;
		string refName;
		int strand;
		int start;
		int end;
		genomeSequences.GetRandomSequence(readLength, sequence, refName, strand, start, end);
		
		cout << "@" << i << "_" << refName << "_" << strand << "_" << start << "_" << end << endl;
		cout << sequence << endl;
		cout << "+" << endl;
		cout << string(readLength, '2') << endl;
	}
}


