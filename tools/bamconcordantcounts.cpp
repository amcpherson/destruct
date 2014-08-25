/*
 *  bamconcordantcounts.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "DebugCheck.h"
#include "Indexer.h"
#include "AlignmentStream.h"
#include "RegionDB.h"
#include "api/BamReader.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;

using namespace BamTools;


class ConcordantReadCounts
{
public:
	explicit ConcordantReadCounts(const vector<string>& refNames, const vector<int>& refLengths) : mRefNames(refNames), mRefLengths(refLengths)
	{
		for (int refID = 0; refID < refNames.size(); refID++)
		{
			mRefIDMap.insert(make_pair(refNames[refID], refID));
		}
	}
	
	void AddVariants(const string& variantFilename)
	{
		ifstream variantFile(variantFilename.c_str());
		CheckFile(variantFile, variantFilename);
		
		StringVec fields;
		int line = 0;
		while (ReadTSV(variantFile, fields))
		{
			line++;
			
			if (fields.size() < 7)
			{
				cerr << "Error: unable to parser " << variantFilename << ":" << line << endl;
				exit(1);
			}
			
			for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
			{
				const string& chromosome = fields[1 + 3 * clusterEnd];
				int breakend = SAFEPARSE(int, fields[2 + 3 * clusterEnd]);
				const string& strand = fields[3 + 3 * clusterEnd];
				
				int refID = mRefIDMap.find(chromosome)->second;
				mPositions[refID].insert(breakend);
			}
		}
	}
	
	void AddTelomeres(const string& telomereFilename)
	{
		ifstream telomereFile(telomereFilename.c_str());
		CheckFile(telomereFile, telomereFilename);
		
		StringVec fields;
		int line = 0;
		while (ReadTSV(telomereFile, fields))
		{
			line++;
			
			if (fields.size() < 4)
			{
				cerr << "Error: unable to parser " << telomereFilename << ":" << line << endl;
				exit(1);
			}
			
			const string& chromosome = fields[1];
			int position = SAFEPARSE(int, fields[2]);
			
			int refID = mRefIDMap.find(chromosome)->second;
			mPositions[refID].insert(position);
		}
	}
	
	void InitializeCounts(int maxFragmentLength)
	{
		for (unordered_map<int,set<int> >::iterator refPositionsIter = mPositions.begin(); refPositionsIter != mPositions.end(); refPositionsIter++)
		{
			int refID = refPositionsIter->first;
			set<int>& positions = refPositionsIter->second;
			
			set<int>::const_iterator positionIter = positions.begin();
			set<int>::const_iterator nextPositionIter = positions.begin();
			nextPositionIter++;
			while (nextPositionIter != positions.end())
			{
				mIntervalCounts[refID].insert(make_pair(make_pair(*positionIter, *nextPositionIter), 0));
				
				positionIter++;
				nextPositionIter++;
			}
			
			set<int>::const_iterator firstPositionIter = positions.begin();
			set<int>::const_iterator lastPositionIter = positions.end();
			lastPositionIter--;
			
			positionIter = positions.begin();
			while (positionIter != positions.end())
			{
				if (positionIter != firstPositionIter && positionIter != lastPositionIter)
				{
					mReferenceCounts[refID].insert(make_pair(make_pair(*positionIter, *positionIter), 0));
				}
				
				positionIter++;
			}
		}
	}
	
	void AddConcordantRead(int refID, int start, int end)
	{
		unordered_map<int,set<int> >::const_iterator refPositionsIter = mPositions.find(refID);
		
		if (refPositionsIter == mPositions.end())
		{
			return;
		}
		
		if ((start < *(refPositionsIter->second.begin()) + 1) || (end > *(refPositionsIter->second.rbegin()) - 1) || end < start)
		{
			cerr << "unexpected read " << mRefNames[refID] << " " << start << " " << end << endl;
			return;
		}
		
		set<int>::const_iterator positionIter = refPositionsIter->second.lower_bound(start);
		
		assert(positionIter != refPositionsIter->second.begin());
		assert(positionIter != refPositionsIter->second.end());
		
		positionIter--;
		
		set<int> localPositions;
		while (*positionIter <= end)
		{
			localPositions.insert(*positionIter);
			positionIter++;
		}
		localPositions.insert(*positionIter);
		
		assert(localPositions.size() >= 2);
		
		if (localPositions.size() == 2)
		{
			set<int>::const_iterator leftIter = localPositions.begin();
			set<int>::const_iterator rightIter = leftIter;
			rightIter++;
			
			unordered_map<pair<int,int>,int>::iterator intervalIter = mIntervalCounts[refID].find(make_pair(*leftIter, *rightIter));
			assert(intervalIter != mIntervalCounts[refID].end());
			
			intervalIter->second++;
		}
		else
		{
			set<int>::const_iterator lastIter = localPositions.end();
			lastIter--;
			
			positionIter = localPositions.begin();
			positionIter++;
			
			while (positionIter != lastIter)
			{
				unordered_map<pair<int,int>,int>::iterator referenceIter = mReferenceCounts[refID].find(make_pair(*positionIter, *positionIter));
				assert(referenceIter != mReferenceCounts[refID].end());
				
				referenceIter->second++;
				
				positionIter++;
			}
		}
	}
	
	void WriteIntervalCounts(const string& filename) const
	{
		WriteCounts(mIntervalCounts, filename, '-', '+');
	}
	
	void WriteReferenceCounts(const string& filename) const
	{
		WriteCounts(mReferenceCounts, filename, '+', '-');
	}
	
private:
	void WriteCounts(const unordered_map<int,unordered_map<pair<int,int>,int> >& counts, const string& filename, char leftStrand, char rightStrand) const
	{
		ofstream file(filename.c_str());
		CheckFile(file, filename);
		
		int id = 0;
		
		for (int refID = 0; refID < mRefNames.size(); refID++)
		{
			if (counts.find(refID) == counts.end())
			{
				continue;
			}
			
			map<pair<int,int>,int> regions(counts.find(refID)->second.begin(), counts.find(refID)->second.end());
			
			for (map<pair<int,int>,int>::const_iterator regionIter = regions.begin(); regionIter != regions.end(); regionIter++)
			{
				file << id << "\t";
				file << mRefNames[refID] << "\t";
				file << regionIter->first.first << "\t";
				file << leftStrand << "\t";
				file << mRefNames[refID] << "\t";
				file << regionIter->first.second << "\t";
				file << rightStrand << "\t";
				file << regionIter->second << endl;
				
				id++;
			}
		}
		
		file.close();
	}
	
	const vector<string>& mRefNames;
	const vector<int>& mRefLengths;
	unordered_map<string,int> mRefIDMap;
	
	unordered_map<int,set<int> > mPositions;
	unordered_map<int,unordered_map<pair<int,int>,int> > mIntervalCounts;
	unordered_map<int,unordered_map<pair<int,int>,int> > mReferenceCounts;
};

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

class PairedBamReader
{
public:
	explicit PairedBamReader(BamReader& bamReader) : mBamReader(bamReader), mRowCount(0) {}
	
	bool Next(BamAlignment& alignment1, BamAlignment& alignment2)
	{
		BamAlignment alignment;
		while (mBamReader.GetNextAlignment(alignment))
		{
			if (++mRowCount % 2000000 == 0)
			{
				cerr << ".";
				cerr.flush();
			}
			
			int readEnd = alignment.IsFirstMate() ? 0 : 1;
			int otherReadEnd = OtherReadEnd(readEnd);
			
			unordered_map<string,BamAlignment>::iterator otherEndIter = mReadBuffer[otherReadEnd].find(alignment.Name);
			
			if (otherEndIter != mReadBuffer[otherReadEnd].end())
			{
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
				
				mReadBuffer[otherReadEnd].erase(otherEndIter);
				
				return true;
			}
			else
			{
				mReadBuffer[readEnd].insert(make_pair(alignment.Name, alignment));
			}
		}
		
		cerr << endl;
		
		return false;
	}
	
private:
	BamReader& mBamReader;
	int mRowCount;
	unordered_map<string,BamAlignment> mReadBuffer[2];
};


int main(int argc, char* argv[])
{
	string variantFilename;
	string telomereFilename;
	string bamFilename;
	int maxSoftClipped;
	int maxFragmentLength;
	string intervalFilename;
	string referenceFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Bam Concordant Read Counter");
		TCLAP::ValueArg<string> variantFilenameArg("v","variant","Variant Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> telomereFilenameArg("t","telomere","Telomere Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> bamFilenameArg("b","bam","Bam Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> maxSoftClippedArg("","clipmax","Maximum Allowable Soft Clipped",true,0,"integer",cmd);
		TCLAP::ValueArg<int> maxFragmentLengthArg("","flen","Maximum Fragment Length",true,0,"integer",cmd);
		TCLAP::ValueArg<string> intervalFilenameArg("i","interval","Interval Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> referenceFilenameArg("r","reference","Reference Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		variantFilename = variantFilenameArg.getValue();
		telomereFilename = telomereFilenameArg.getValue();
		bamFilename = bamFilenameArg.getValue();
		maxSoftClipped = maxSoftClippedArg.getValue();
		maxFragmentLength = maxFragmentLengthArg.getValue();
		intervalFilename = intervalFilenameArg.getValue();
		referenceFilename = referenceFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	BamReader bamReader;
	if (!bamReader.Open(bamFilename))
	{
		cerr << "Error: Unable to open bam file " << bamFilename << endl;
		exit(1);
	}
	
	PairedBamReader pairedBamReader(bamReader);
	
	vector<string> refNames;
	vector<int> refLengths;
	for (int refID = 0; refID < bamReader.GetReferenceCount(); refID++)
	{
		refNames.push_back(bamReader.GetReferenceData()[refID].RefName);
		refLengths.push_back(bamReader.GetReferenceData()[refID].RefLength);
	}
	
	ConcordantReadCounts readCounts(refNames, refLengths);
	
	cerr << "Adding breakpoints" << endl;
	
	readCounts.AddVariants(variantFilename);
	
	cerr << "Adding telomeres" << endl;
	
	readCounts.AddTelomeres(telomereFilename);
	
	cerr << "Initializing counts" << endl;
	
	readCounts.InitializeCounts(maxFragmentLength);
	
	cerr << "Counting reads from bam" << endl;
	
	int rowCount = 0;
	
	BamAlignment alignment1;
	BamAlignment alignment2;
	while (pairedBamReader.Next(alignment1, alignment2))
	{
		bool properPair = alignment1.IsProperPair() && alignment2.IsProperPair() && abs(alignment1.InsertSize) <= maxFragmentLength;
		bool concordant = properPair && GetNumSoftClipped(alignment1) <= maxSoftClipped && GetNumSoftClipped(alignment2) <= maxSoftClipped;
		bool failedQC = alignment1.IsFailedQC() || alignment2.IsFailedQC();
		bool failedMapQual = alignment1.MapQuality <= 0 || alignment2.MapQuality <= 0;
		
		if (concordant && !failedQC && !failedMapQual)
		{
			int fragmentStart = min(alignment1.Position, alignment2.Position);
			int fragmentEnd = fragmentStart + abs(alignment1.InsertSize) - 1;
			
			readCounts.AddConcordantRead(alignment1.RefID, fragmentStart, fragmentEnd);
		}
	}
	
	readCounts.WriteIntervalCounts(intervalFilename);
	readCounts.WriteReferenceCounts(referenceFilename);
}
