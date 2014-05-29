/*
 *  bamextractsimreads.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "DebugCheck.h"
#include "Indexer.h"
#include "AlignmentStream.h"
#include "RegionDB.h"
#include "Parsers.h"
#include "Sequences.h"
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_utilities.h"
#include "utils/bamtools_fasta.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace boost;
using namespace std;

using namespace BamTools;


// Get read end as 0 or 1
inline int GetReadEnd(const BamAlignment& alignment)
{
	return alignment.IsFirstMate() ? 0 : 1;
}

// Get sequence from an alignment accounting for strand
inline string GetSequence(const BamAlignment& alignment)
{
	string sequence = alignment.QueryBases;
	if (alignment.IsReverseStrand())
	{
		ReverseComplement(sequence);
	}
	return sequence;
}

// Get qualities from an alignment accounting for strand
inline string GetQualities(const BamAlignment& alignment)
{
	string qualities = alignment.Qualities;
	if (alignment.IsReverseStrand())
	{
		reverse(qualities.begin(), qualities.end());
	}
	return qualities;
}

struct BamSimReader : PileupVisitor
{
	BamSimReader(const string& bamFilename, const string& fastaFilename)
	{
		if (!mBamReader.Open(bamFilename))
		{
			cerr << "Error: Unable to open bam file " << bamFilename << endl;
			exit(1);
		}
		
		if (!mBamReader.LocateIndex())
		{
			cerr << "Error: Unable to find index for bam file " << bamFilename << endl;
			exit(1);
		}

		string fastaIndexFilename = "";
		if (Utilities::FileExists(fastaFilename + ".fai"))
		{
			fastaIndexFilename = fastaFilename + ".fai";
		}
		
		if (!mFasta.Open(fastaFilename, fastaIndexFilename))
		{
			cerr << "Error: Unable to open reference fasta file " << fastaFilename << endl;
			exit(1);
		}
	}
	
	void Read(const string& chromosome, int position, const string& simSequence)
	{
		mStartPos = position;
		mEndPos = position + simSequence.size();
		mSimSequence = simSequence;

		// Set the region for the bam index
		int bamRefID = mBamReader.GetReferenceID(chromosome);
		mBamReader.SetRegion(BamRegion(bamRefID, position, bamRefID, mEndPos));

		// Get the fasta chromosome index
		vector<string> fastaReferenceNames = mFasta.GetReferenceNames();
		vector<string>::const_iterator faChrIter = find(fastaReferenceNames.begin(), fastaReferenceNames.end(), chromosome);
		if (faChrIter == fastaReferenceNames.end())
		{
			cerr << "Error: Unable to find chromosome " << chromosome << " in fasta" << endl;
			exit(1);
		}
		mFastaRefID = faChrIter - fastaReferenceNames.begin();

		PileupEngine pileupEngine;
		
		pileupEngine.AddVisitor(dynamic_cast<PileupVisitor*>(this));

		BamAlignment alignment;
		while (mBamReader.GetNextAlignment(alignment))
		{
			if (alignment.Position < mStartPos || alignment.GetEndPosition() >= mEndPos)
			{
				continue;
			}

			int readEnd = GetReadEnd(alignment);

			mAlignments[readEnd][alignment.Name] = alignment;

			pileupEngine.AddAlignment(alignment);
		}

		pileupEngine.Flush();

		for (unordered_map<string,BamAlignment>::const_iterator alignmentIter = mAlignments[0].begin(); alignmentIter != mAlignments[0].end(); alignmentIter++)
		{
			if (mAlignments[1].find(alignmentIter->first) != mAlignments[1].end())
			{
				mReadNames.push_back(alignmentIter->first);
			}
		}
	}

	void Visit(const PileupPosition& pileupData)
	{		
		int localPosition = pileupData.Position - mStartPos;

		if (localPosition < 0 || localPosition >= mSimSequence.size())
		{
			assert(false);
			return;
		}

		char refBase;
		mFasta.GetBase(mFastaRefID, mStartPos, refBase);

		char simBase = mSimSequence[localPosition];

		refBase = toupper(refBase);
		simBase = toupper(simBase);

		for (vector<PileupAlignment>::const_iterator pileupIter = pileupData.PileupAlignments.begin(); pileupIter != pileupData.PileupAlignments.end(); ++pileupIter)
		{
			const PileupAlignment& pileupAlignment = (*pileupIter);
			const BamAlignment& alignment = pileupAlignment.Alignment;
			
			if (pileupAlignment.IsCurrentDeletion)
			{
				continue;
			}

			char base = toupper(alignment.QueryBases.at(pileupAlignment.PositionInAlignment));
			
			// By default set the base to that of the simulation base at this position
			char simReadBase = simBase;

			// If there is a mismatch in the read at this position, produce a mismatch in the new read
			if (refBase != base)
			{
				// Try the original base of the read first
				simReadBase = base;

				// If the original base equals the simulated sequence base, move to the next base
				if (simReadBase == simBase)
				{
					switch (simReadBase)
					{
						case 'A': simReadBase = 'C'; break;
						case 'C': simReadBase = 'T'; break;
						case 'T': simReadBase = 'G'; break;
						case 'G': simReadBase = 'A'; break;
						default: simReadBase = 'N';
					}
				}
			}

			int readEnd = GetReadEnd(alignment);

			mModifications[readEnd][alignment.Name].push_back(make_pair(pileupAlignment.PositionInAlignment, simReadBase));
		}
	}

	void WriteFastq(const string& namePrefix, int readEnd, ostream& fastq)
	{
		for (int readID = 0; readID < mReadNames.size(); readID++)
		{
			const string& readName = mReadNames[readID];
			BamAlignment alignment = mAlignments[readEnd][readName];

			const vector<pair<int,char> >& modifications = mModifications[readEnd][readName];

			for (vector<pair<int,char> >::const_iterator modIter = modifications.begin(); modIter != modifications.end(); modIter++)
			{
				alignment.QueryBases[modIter->first] = modIter->second;
			}

			fastq << "@" << namePrefix << "_" << readID << "/" << readID+1 << endl;
			fastq << GetSequence(alignment) << endl;
			fastq << "+" << readName << endl;
			fastq << GetQualities(alignment) << endl;
		}
	}

	BamReader mBamReader;
	Fasta mFasta;
	
	int mFastaRefID;
	int mStartPos;
	int mEndPos;

	string mSimSequence;

	vector<string> mReadNames;
	unordered_map<string,BamAlignment> mAlignments[2];
	unordered_map<string,vector<pair<int,char> > > mModifications[2];
};


int main(int argc, char* argv[])
{
	string bamFilename;
	string refFasta;
	string seqFasta;
	string fastq1Filename;
	string fastq2Filename;
	
	try
	{
		TCLAP::CmdLine cmd("Bam Read Simulator");
		TCLAP::ValueArg<string> bamFilenameArg("b","bam","Bam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> refFastaArg("r","ref","Reference Fasta Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> seqFastaArg("s","seq","Sequence to match reads to",true,"","string",cmd);
		TCLAP::ValueArg<string> fastq1FilenameArg("1","fastq1","Fastq 1 filename",true,"","string",cmd);
		TCLAP::ValueArg<string> fastq2FilenameArg("2","fastq2","Fastq 2 filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		bamFilename = bamFilenameArg.getValue();
		refFasta = refFastaArg.getValue();
		seqFasta = seqFastaArg.getValue();
		fastq1Filename = fastq1FilenameArg.getValue();
		fastq2Filename = fastq2FilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	BamSimReader bamSimReader(bamFilename, refFasta);

	ofstream fastq1File(fastq1Filename.c_str());
	ofstream fastq2File(fastq2Filename.c_str());

	CheckFile(fastq1File, fastq1Filename);
	CheckFile(fastq2File, fastq2Filename);

	Sequences simSequences;
	simSequences.Read(seqFasta);

	vector<long> chromosomeLengths;
	for (int chrIdx = 0; chrIdx < bamSimReader.mBamReader.GetReferenceCount(); chrIdx++)
	{
		chromosomeLengths.push_back(bamSimReader.mBamReader.GetReferenceData()[chrIdx].RefLength);
	}
	
	RandomGenomicPositionGenerator randomPosition(chromosomeLengths);
	
	const vector<string>& seqNames = simSequences.GetNames();
	for (int nameIdx = 0; nameIdx < seqNames.size(); nameIdx++)
	{
		const string& seqName = seqNames[nameIdx];
		const string& simSeq = simSequences.Get(seqName);

		int chrIdx;
		long position;
		randomPosition.Next(chrIdx, position);
		
		const string& chromosome = bamSimReader.mBamReader.GetReferenceData()[chrIdx].RefName;

		bamSimReader.Read(chromosome, position, simSeq);

		bamSimReader.WriteFastq(seqName, 0, fastq1File);
		bamSimReader.WriteFastq(seqName, 1, fastq2File);
	}
}


