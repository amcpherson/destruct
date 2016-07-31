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
#include "api/BamReader.h"
#include "DiskPriorityQueue.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

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

template<typename T>
struct ReservoirSampler
{
	ReservoirSampler(int numSamples) : mNumSamples(numSamples), mNumValues(0) {}

	void AddSample(const T& value)
	{
		mNumValues++;

		if (mSamples.size() < mNumSamples)
		{
			mSamples.push_back(value);
		}
		else
		{
			int sampleIndex = mRNG.Next(0, mNumValues - 1);
			if (sampleIndex < mNumSamples)
			{
				mSamples[sampleIndex] = value;
			}
		}
	}

	int mNumSamples;
	int mNumValues;
	vector<T> mSamples;
	RandomNumberGenerator mRNG;	
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
	string sample1Filename;
	string sample2Filename;
	int numSamples;
	bool renameReads;
	
	try
	{
		TCLAP::CmdLine cmd("Bam to Fastq Tool");
		TCLAP::ValueArg<string> bamFilenameArg("b","bam","Bam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> fastq1FilenameArg("","fastq1","Fastq End 1 Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> fastq2FilenameArg("","fastq2","Fastq End 2 Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> maxSoftClippedArg("c","clipmax","Maximum Allowable Soft Clipped",true,0,"integer",cmd);
		TCLAP::ValueArg<int> maxFragmentLengthArg("f","flen","Maximum Fragment Length",true,0,"integer",cmd);
		TCLAP::ValueArg<string> statsFilenameArg("s","stats","Concordant Stats Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> tempsPrefixArg("t","temp","Filename Prefix for Temporary Files",true,"","string",cmd);
		TCLAP::ValueArg<string> sample1FilenameArg("","sample1","Sample Fastq End 1 Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> sample2FilenameArg("","sample2","Sample Fastq End 2 Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> numSamplesArg("n","num","Number of Samples",true,0,"integer",cmd);
		TCLAP::SwitchArg renameReadsArg("r","rename","Rename With Integer IDs",cmd);
		cmd.parse(argc,argv);
		
		bamFilename = bamFilenameArg.getValue();
		fastq1Filename = fastq1FilenameArg.getValue();
		fastq2Filename = fastq2FilenameArg.getValue();
		maxSoftClipped = maxSoftClippedArg.getValue();
		maxFragmentLength = maxFragmentLengthArg.getValue();
		statsFilename = statsFilenameArg.getValue();
		tempsPrefix = tempsPrefixArg.getValue();
		sample1Filename = sample1FilenameArg.getValue();
		sample2Filename = sample2FilenameArg.getValue();
		numSamples = numSamplesArg.getValue();
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

	ofstream fastq1File(fastq1Filename.c_str(), std::ios_base::out | std::ios_base::binary);
	ofstream fastq2File(fastq2Filename.c_str(), std::ios_base::out | std::ios_base::binary);

	CheckFile(fastq1File, fastq1Filename);
	CheckFile(fastq2File, fastq2Filename);
	
	iostreams::filtering_ostream fastq1Stream;
	iostreams::filtering_ostream fastq2Stream;

	fastq1Stream.push(iostreams::gzip_compressor());
	fastq2Stream.push(iostreams::gzip_compressor());

	fastq1Stream.push(fastq1File);
	fastq2Stream.push(fastq2File);

	ReservoirSampler<pair<ReadInfo,ReadInfo> > sampledReads(numSamples);

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

		sampledReads.AddSample(make_pair(ReadInfo(alignment1), ReadInfo(alignment2)));

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
			
			fastq1Stream << "@" << fragment << "/1" << endl;
			fastq1Stream << GetSequence(alignment1) << endl;
			fastq1Stream << "+" << alignment1.Name << endl;
			fastq1Stream << GetQualities(alignment1) << endl;
			
			fastq2Stream << "@" << fragment << "/2" << endl;
			fastq2Stream << GetSequence(alignment2) << endl;
			fastq2Stream << "+" << alignment2.Name << endl;
			fastq2Stream << GetQualities(alignment2) << endl;
			
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

		sampledReads.AddSample(make_pair(discordantRead1, discordantRead2));

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
		
		fastq1Stream << "@" << fragment << "/1" << endl;
		fastq1Stream << discordantRead1.Sequence << endl;
		fastq1Stream << "+" << discordantRead1.Name << endl;
		fastq1Stream << discordantRead1.Qualities << endl;
		
		fastq2Stream << "@" << fragment << "/2" << endl;
		fastq2Stream << discordantRead2.Sequence << endl;
		fastq2Stream << "+" << discordantRead2.Name << endl;
		fastq2Stream << discordantRead2.Qualities << endl;
		
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

	// Write out sampled reads
	ofstream sample1File(sample1Filename.c_str(), std::ios_base::out | std::ios_base::binary);
	ofstream sample2File(sample2Filename.c_str(), std::ios_base::out | std::ios_base::binary);

	CheckFile(sample1File, sample1Filename);
	CheckFile(sample2File, sample2Filename);
	
	iostreams::filtering_ostream sample1Stream;
	iostreams::filtering_ostream sample2Stream;

	sample1Stream.push(iostreams::gzip_compressor());
	sample2Stream.push(iostreams::gzip_compressor());

	sample1Stream.push(sample1File);
	sample2Stream.push(sample2File);

	for (int sampleIndex = 0; sampleIndex < sampledReads.mSamples.size(); sampleIndex++)
	{
		const ReadInfo& read1 = sampledReads.mSamples[sampleIndex].first;
		const ReadInfo& read2 = sampledReads.mSamples[sampleIndex].second;

		sample1Stream << "@" << sampleIndex << "/1" << endl;
		sample1Stream << read1.Sequence << endl;
		sample1Stream << "+" << read1.Name << endl;
		sample1Stream << read1.Qualities << endl;
		
		sample2Stream << "@" << sampleIndex << "/2" << endl;
		sample2Stream << read2.Sequence << endl;
		sample2Stream << "+" << read2.Name << endl;
		sample2Stream << read2.Qualities << endl;
	}
}

