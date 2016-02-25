/*
 *  bamdiscordantfastq.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "DebugCheck.h"
#include "Indexer.h"
#include "AlignmentStream.h"
#include "RegionDB.h"
#include "DiskPriorityQueue.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <bamtools/api/BamReader.h>
#include <boost/algorithm/string.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

using namespace boost;
using namespace std;

using namespace BamTools;

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

// Calculate number of soft clipped reads for an alignment
inline int GetNumSoftClipped(const BamAlignment& alignment)
{
	int numSoftClipped = 0;
	for (vector<CigarOp>::const_iterator cigarOpIter = alignment.CigarData.begin(); cigarOpIter != alignment.CigarData.end(); cigarOpIter++)
	{
		if (cigarOpIter->Type == 'S')
		{
			numSoftClipped += cigarOpIter->Length;
		}
	}
	return numSoftClipped;
}

// Check if an alignment is concordant accoring to proper pair flag, insert size, and soft clipping
inline bool IsConcordant(const BamAlignment& alignment, int maxFragmentLength, int maxSoftClipped)
{
	bool properPair = alignment.IsProperPair() && abs(alignment.InsertSize) <= maxFragmentLength;
	bool concordant = properPair && GetNumSoftClipped(alignment) <= maxSoftClipped;
	return concordant;
}

// Store minimal read info
struct ReadInfo
{
	ReadInfo() {}

	ReadInfo(const BamAlignment& alignment)
		: Name(alignment.Name), 
		  Sequence(GetSequence(alignment)),
		  Qualities(GetQualities(alignment)),
		  IsFailedQC(alignment.IsFailedQC())
	{
	}

	string Name;
	string Sequence;
	string Qualities;
	bool IsFailedQC;

	// Used for sorting by read name
	bool operator<(const ReadInfo& other) const
	{
		return Name < other.Name;
	}

	// Serialization for on disk sort
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & Name;
        ar & Sequence;
        ar & Qualities;
        ar & IsFailedQC;
    }
};

class PairedBamReader
{
public:
	PairedBamReader(BamReader& bamReader, const string& tempsPrefix)
		: mBamReader(bamReader),
		  mBamReadFinished(false),
		  mDiscordantReadQueue1(tempsPrefix + "_1_", 1024*1024),
		  mDiscordantReadQueue2(tempsPrefix + "_2_", 1024*1024)
	{
	}
	
	bool NextConcordant(BamAlignment& alignment1, BamAlignment& alignment2)
	{
		if (!mBamReadFinished)
		{
			BamAlignment alignment;
			while (mBamReader.GetNextAlignment(alignment))
			{
				if (alignment.IsProperPair())
				{
					// Proper pairs should be close to each other in the bam file
					// store proper pairs in a buffer and match them on the fly

					int readEnd = alignment.IsFirstMate() ? 0 : 1;
					int otherReadEnd = 1 - readEnd;

					// Search for other end already in the buffer
					unordered_map<string,BamAlignment>::iterator otherEndIter = mConcordantReadBuffer[otherReadEnd].find(alignment.Name);
					
					if (otherEndIter != mConcordantReadBuffer[otherReadEnd].end())
					{
						// Return this alignment and other end alignment
						// erase other end from the buffer

						if (alignment.IsFirstMate())
						{
							alignment1 = alignment;
							alignment2 = otherEndIter->second;
						}
						else
						{
							alignment1 = otherEndIter->second;
							alignment2 = alignment;
						}
						
						mConcordantReadBuffer[otherReadEnd].erase(otherEndIter);
						
						return true;
					}
					else
					{
						// Insert alignment into the buffer
						mConcordantReadBuffer[readEnd].insert(make_pair(alignment.Name, alignment));
					}
				}
				else
				{
					// Improper pairs are not guaranteed to be close to each other
					// store improper pairs in an on disk priority queue

					if (alignment.IsFirstMate())
					{
						mDiscordantReadQueue1.Push(ReadInfo(alignment));
					}
					else
					{
						mDiscordantReadQueue2.Push(ReadInfo(alignment));
					}
				}
			}

			// If we just read the last alignment, finalize the on disk priority queue
			mDiscordantReadQueue1.Finalize();
			mDiscordantReadQueue2.Finalize();

			mBamReadFinished = true;
		}

		return false;
	}

	bool NextDiscordant(ReadInfo& read1, ReadInfo& read2)
	{
		assert(mBamReadFinished);

		// Tandem iteration to match pairs of reads
		while (!mDiscordantReadQueue1.Empty() && !mDiscordantReadQueue2.Empty())
		{
			if (mDiscordantReadQueue1.Top().Name < mDiscordantReadQueue2.Top().Name)
			{
				mDiscordantReadQueue1.Pop();
			}
			else if (mDiscordantReadQueue1.Top().Name > mDiscordantReadQueue2.Top().Name)
			{
				mDiscordantReadQueue2.Pop();
			}
			else
			{
				// Return matched read pair and pop from queue

				read1 = mDiscordantReadQueue1.Top();
				read2 = mDiscordantReadQueue2.Top();

				mDiscordantReadQueue1.Pop();
				mDiscordantReadQueue2.Pop();

				return true;
			}
		}

		return false;
	}
	
private:
	BamReader& mBamReader;
	bool mBamReadFinished;

	unordered_map<string,BamAlignment> mConcordantReadBuffer[2];

	DiskPriorityQueue<ReadInfo> mDiscordantReadQueue1;
	DiskPriorityQueue<ReadInfo> mDiscordantReadQueue2;
};

int main(int argc, char* argv[])
{
	string bamFilename;
	string fastq1Filename;
	string fastq2Filename;
	int maxSoftClipped;
	int maxFragmentLength;
	string statsFilename;
	string tempsPrefix;
	bool renameReads;
	
	try
	{
		TCLAP::CmdLine cmd("Bam to Fastq Tool");
		TCLAP::ValueArg<string> bamFilenameArg("b","bam","Bam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> fastq1FilenameArg("1","fastq1","Fastq End 1 Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> fastq2FilenameArg("2","fastq2","Fastq End 2 Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> maxSoftClippedArg("c","clipmax","Maximum Allowable Soft Clipped",true,0,"integer",cmd);
		TCLAP::ValueArg<int> maxFragmentLengthArg("f","flen","Maximum Fragment Length",true,0,"integer",cmd);
		TCLAP::ValueArg<string> statsFilenameArg("s","stats","Concordant Stats Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> tempsPrefixArg("t","temp","Filename Prefix for Temporary Files",true,"","string",cmd);
		TCLAP::SwitchArg renameReadsArg("r","rename","Rename With Integer IDs",cmd);
		cmd.parse(argc,argv);
		
		bamFilename = bamFilenameArg.getValue();
		fastq1Filename = fastq1FilenameArg.getValue();
		fastq2Filename = fastq2FilenameArg.getValue();
		maxSoftClipped = maxSoftClippedArg.getValue();
		maxFragmentLength = maxFragmentLengthArg.getValue();
		statsFilename = statsFilenameArg.getValue();
		tempsPrefix = tempsPrefixArg.getValue();
		renameReads = renameReadsArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	int concordantReadCount = 0;
	int discordantReadCount = 0;
	unordered_map<int,int> readLengthHist;
	unordered_map<int,int> fragmentLengthHist;
	
	BamReader bamReader;
	if (!bamReader.Open(bamFilename))
	{
		cerr << "Error: Unable to open bam file " << bamFilename << endl;
		exit(1);
	}

	PairedBamReader pairedReader(bamReader, tempsPrefix);

	ofstream fastq1File(fastq1Filename.c_str());
	ofstream fastq2File(fastq2Filename.c_str());
	
	CheckFile(fastq1File, fastq1Filename);
	CheckFile(fastq2File, fastq2Filename);
	
	int fragmentIndex = 0;

	BamAlignment alignment1;
	BamAlignment alignment2;
	while (pairedReader.NextConcordant(alignment1, alignment2))
	{
		// Ignore all failed reads
		if (alignment1.IsFailedQC() || alignment2.IsFailedQC())
		{
			continue;
		}

		concordantReadCount++;

		// Update read length histogram for all reads
		readLengthHist.insert(make_pair(alignment1.Length, 0)).first->second++;
		readLengthHist.insert(make_pair(alignment2.Length, 0)).first->second++;

		if (IsConcordant(alignment1, maxFragmentLength, maxSoftClipped) && IsConcordant(alignment2, maxFragmentLength, maxSoftClipped))
		{
			// Update fragment length histogram for concordant reads
			fragmentLengthHist.insert(make_pair(abs(alignment1.InsertSize), 0)).first->second++;
		}
		else
		{
			// Optionally change the fragment name
			string fragment = alignment1.Name;
			if (renameReads)
			{
				stringstream fragmentStream;
				fragmentStream << fragmentIndex;
				fragment = fragmentStream.str();
			}

			// Write fastq
			
			fastq1File << "@" << fragment << "/1" << endl;
			fastq1File << GetSequence(alignment1) << endl;
			fastq1File << "+" << alignment1.Name << endl;
			fastq1File << GetQualities(alignment1) << endl;
			
			fastq2File << "@" << fragment << "/2" << endl;
			fastq2File << GetSequence(alignment2) << endl;
			fastq2File << "+" << alignment2.Name << endl;
			fastq2File << GetQualities(alignment2) << endl;
			
			fragmentIndex++;
		}
	}

	ReadInfo discordantRead1;
	ReadInfo discordantRead2;
	while (pairedReader.NextDiscordant(discordantRead1, discordantRead2))
	{
		// Ignore all failed reads
		if (discordantRead1.IsFailedQC || discordantRead2.IsFailedQC)
		{
			continue;
		}

		discordantReadCount++;

		// Update read length histogram for all reads
		readLengthHist.insert(make_pair(alignment1.Length, 0)).first->second++;
		readLengthHist.insert(make_pair(alignment2.Length, 0)).first->second++;

		// Optionally change the fragment name
		string fragment = discordantRead1.Name;
		if (renameReads)
		{
			stringstream fragmentStream;
			fragmentStream << fragmentIndex;
			fragment = fragmentStream.str();
		}

		// Write fastq
		
		fastq1File << "@" << fragment << "/1" << endl;
		fastq1File << discordantRead1.Sequence << endl;
		fastq1File << "+" << discordantRead1.Name << endl;
		fastq1File << discordantRead1.Qualities << endl;
		
		fastq2File << "@" << fragment << "/2" << endl;
		fastq2File << discordantRead2.Sequence << endl;
		fastq2File << "+" << discordantRead2.Name << endl;
		fastq2File << discordantRead2.Qualities << endl;
		
		fragmentIndex++;
	}

	// Check for an empty bam file (fail, somethings wrong)
	if (concordantReadCount + discordantReadCount == 0)
	{
		cerr << "Error: No reads" << endl;
		exit(1);
	}
	
	// Check for a bam file with no concordant reads (fail, somethings wrong)
	if (concordantReadCount == 0)
	{
		cerr << "Error: No concordant reads" << endl;
		exit(1);
	}
	
	// Check for a bam file with no discordant reads (usually somethings wrong)
	if (discordantReadCount == 0)
	{
		cerr << "Error: No discordant reads" << endl;
		exit(1);
	}
	
	// Output stats
	ofstream statsFile(statsFilename.c_str());
	CheckFile(statsFile, statsFilename);
	statsFile << "type\tkey\tvalue\n";

	// Output read counts
	statsFile << "read_count\ttotal\t" << concordantReadCount + discordantReadCount << endl;
	statsFile << "read_count\tconcordant\t" << concordantReadCount << endl;
	statsFile << "read_count\tdiscordant\t" << discordantReadCount << endl;

	// Output read lengths
	for (unordered_map<int,int>::const_iterator iter = readLengthHist.begin(); iter != readLengthHist.end(); iter++)
	{
		statsFile << "read_length\t" << iter->first << "\t" << iter->second << endl;
	}
	
	// Output fragment lengths
	for (unordered_map<int,int>::const_iterator iter = fragmentLengthHist.begin(); iter != fragmentLengthHist.end(); iter++)
	{
		statsFile << "fragment_length\t" << iter->first << "\t" << iter->second << endl;
	}
}

