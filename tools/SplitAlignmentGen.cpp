/*
 *  SplitAlignment.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "SplitAlignmentGen.h"
#include "DebugCheck.h"
#include "SplitReadAligner.h"

#include <map>
#include <list>
#include <fstream>
#include <sstream>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;

const int numBreakPadding = 10;
const int matchScore = 2;
const int mismatchScore = -1;
const int gapScore = -2;
const int minAnchor = 4;

bool SplitAlignment::Initialize(const LocationVec& alignPair, double fragmentLengthMean, double fragmentLengthStdDev, int minReadLength, int maxReadLength)
{
	int minFragmentLength = (int)(fragmentLengthMean - 3 * fragmentLengthStdDev);
	int maxFragmentLength = (int)(fragmentLengthMean + 3 * fragmentLengthStdDev);

	if (alignPair.size() != 2)
	{
		cerr << "Error: Incorrect input for SplitAlignment::Calculate()" << endl;
		return false;
	}
	
	for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
	{
		int alignStart = alignPair[clusterEnd].start;
		int alignEnd = alignPair[clusterEnd].end;
		
		mAlignRefName[clusterEnd] = alignPair[clusterEnd].refName;
		mAlignStrand[clusterEnd] = alignPair[clusterEnd].strand;
		mSplitSeqStrand[clusterEnd] = (clusterEnd == 0) ? mAlignStrand[clusterEnd] : OtherStrand(mAlignStrand[clusterEnd]);
		
		bool revCompReads = (clusterEnd == 0) ? 1 : 0;
		
		int breakRegionStart;
		int breakRegionLength;
		CalculateBreakRegion(minReadLength, maxReadLength, maxFragmentLength, alignStart, alignEnd, mAlignStrand[clusterEnd], breakRegionStart, breakRegionLength);
		
		// Assumption: the break region we get from the paired end analysis will not involve a read that pushes more than 50%
		// into the breakpoint.  This is aligner dependent.
		if (mAlignStrand[clusterEnd] == PlusStrand)
		{
			mSplitAlignSeqStart[clusterEnd] = breakRegionStart - maxReadLength;
			mSplitAlignSeqLength[clusterEnd] = breakRegionLength + maxReadLength;
		}
		else
		{
			mSplitAlignSeqStart[clusterEnd] = breakRegionStart - breakRegionLength + 1;
			mSplitAlignSeqLength[clusterEnd] = breakRegionLength + maxReadLength;
		}
		
		const string& chromosome = mAlignRefName[clusterEnd];
		int genomeAlignStrand = mAlignStrand[clusterEnd];
		int genomeBreakRegionStart = breakRegionStart;
		
		// Find min and max range of mates
		int mateMin = minFragmentLength - breakRegionLength - maxReadLength + 1;
		int mateMax = maxFragmentLength - minReadLength;
		
		// Find start and end of mate region in the genome
		Region genomeMateRegion;
		if (genomeAlignStrand == PlusStrand)
		{
			genomeMateRegion.start = genomeBreakRegionStart - mateMax;
			genomeMateRegion.end = genomeBreakRegionStart - mateMin;
		}
		else
		{
			genomeMateRegion.start = genomeBreakRegionStart + mateMin;
			genomeMateRegion.end = genomeBreakRegionStart + mateMax;
		}
		
		mMateRegions[clusterEnd].refName = chromosome;
		mMateRegions[clusterEnd].strand = genomeAlignStrand;
		mMateRegions[clusterEnd].start = genomeMateRegion.start;
		mMateRegions[clusterEnd].end = genomeMateRegion.end;
	}
	
	return true;
}

class BinnedLocations
{
public:
	BinnedLocations(int binSpacing) : mBinSpacing(binSpacing) {}
	
	void Add(int id, const Location& location)
	{
		int startBin = location.start / mBinSpacing;
		int endBin = location.end / mBinSpacing;
		
		for (int bin = startBin; bin <= endBin; bin++)
		{
			mBinned[location.strand][location.refName][bin].push_back(id);
		}
		
		mRegions[id].start = location.start;
		mRegions[id].end = location.end;		
	}
	
	void Overlapping(const RawAlignment& alignment, unordered_set<int>& indices) const
	{
		unordered_map<string,unordered_map<int,IntegerVec> >::const_iterator findRefIter = mBinned[alignment.strand].find(alignment.reference);
		if (findRefIter != mBinned[alignment.strand].end())
		{
			int startBin = alignment.region.start / mBinSpacing;
			int endBin = alignment.region.end / mBinSpacing;
			
			for (int bin = startBin; bin <= endBin; bin++)
			{
				unordered_map<int,IntegerVec>::const_iterator findBinIter = findRefIter->second.find(bin);
				if (findBinIter != findRefIter->second.end())
				{
					for (IntegerVecConstIter iter = findBinIter->second.begin(); iter != findBinIter->second.end(); iter++)
					{
						const Region& region = mRegions.find(*iter)->second;
						
						if (region.start <= alignment.region.end && region.end >= alignment.region.start)
						{
							indices.insert(*iter);
						}
					}
				}
			}
		}
	}
	
private:	
	int mBinSpacing;
	unordered_map<string,unordered_map<int,IntegerVec> > mBinned[2];
	unordered_map<int,Region> mRegions;
};

bool SplitAlignment::FindCandidates(AlignmentStream* alignments, SplitAlignmentMap& splitAlignments)
{
	BinnedLocations binnedMateRegions(2000);
	
	IntegerVec ids;
	IntegerVec revComps;
	for (SplitAlignmentMapConstIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		const SplitAlignment& splitAlignment = splitAlignIter->second;
		
		for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
		{
			int revCompReads = (clusterEnd == 0) ? 1 : 0;
			
			int mateRegionID = ids.size();
			
			ids.push_back(id);
			revComps.push_back(revCompReads);
			
			binnedMateRegions.Add(mateRegionID, splitAlignment.mMateRegions[clusterEnd]);
		}
	}
	
	unordered_map<int,unordered_set<IntegerPair> > candidateUnique;
	
	RawAlignment alignment;
	while (alignments->GetNextAlignment(alignment))
	{
		unordered_set<int> overlapping;
		binnedMateRegions.Overlapping(alignment, overlapping);
		
		for (unordered_set<int>::const_iterator olapIter = overlapping.begin(); olapIter != overlapping.end(); olapIter++)
		{
			int id = ids[*olapIter];
			int revComp = revComps[*olapIter];
			
			ReadID candidateReadID;
			candidateReadID.fragmentIndex = SAFEPARSE(int, alignment.fragment);
			candidateReadID.readEnd = (alignment.readEnd == 0) ? 1 : 0;
			
			if (candidateUnique[id].insert(IntegerPair(candidateReadID.id,revComp)).second)
			{
				splitAlignments[id].mCandidateReadID.push_back(candidateReadID.id);
				splitAlignments[id].mCandidateRevComp.push_back(revComp);
			}
		}
	}
}

void SplitAlignment::ReadCandidateSequences(IReadStream* readStream, SplitAlignmentMap& splitAlignments, StringVec& sequences)
{
	unordered_map<int,IntegerPairVec> candidateReadMap;
	
	for (SplitAlignmentMapIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		SplitAlignment& splitAlignment = splitAlignIter->second;
		
		splitAlignment.mCandidateSeqIndex.resize(splitAlignment.mCandidateReadID.size());
		
		for (int candidateIndex = 0; candidateIndex < splitAlignment.mCandidateReadID.size(); candidateIndex++)
		{
			int readID = splitAlignment.mCandidateReadID[candidateIndex];
			candidateReadMap[readID].push_back(IntegerPair(id,candidateIndex));
		}
	}
	
	RawRead rawRead;
	while (readStream->GetNextRead(rawRead))
	{
		ReadID readID;
		readID.fragmentIndex = SAFEPARSE(int, rawRead.fragment);
		readID.readEnd = rawRead.readEnd;
		
		if (candidateReadMap.find(readID.id) == candidateReadMap.end())
		{
			continue;
		}
		
		for (IntegerPairVecConstIter candidateIter = candidateReadMap[readID.id].begin(); candidateIter != candidateReadMap[readID.id].end(); candidateIter++)
		{
			int id = candidateIter->first;
			int candidateIndex = candidateIter->second;
			
			splitAlignments.find(id)->second.mCandidateSeqIndex[candidateIndex] = sequences.size();
		}
		
		sequences.push_back(rawRead.sequence);
	}
}

void SplitAlignment::ReadCandidateSequences(const ReadIndex& readIndex, SplitAlignmentMap& splitAlignments, StringVec& sequences)
{
	unordered_set<int> readIDs;
	for (SplitAlignmentMapIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		SplitAlignment& splitAlignment = splitAlignIter->second;
		
		readIDs.insert(splitAlignment.mCandidateReadID.begin(), splitAlignment.mCandidateReadID.end());
	}
	
	unordered_map<int,int> sequenceIndices;
	for (unordered_set<int>::const_iterator readIDIter = readIDs.begin(); readIDIter != readIDs.end(); readIDIter++)
	{
		ReadID readID;
		readID.id = *readIDIter;
		
		sequenceIndices[readID.id] = sequences.size();
		sequences.push_back(string());
		readIndex.Find(readID.fragmentIndex, readID.readEnd, sequences.back());
	}
	
	for (SplitAlignmentMapIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		SplitAlignment& splitAlignment = splitAlignIter->second;
		
		splitAlignment.mCandidateSeqIndex.resize(splitAlignment.mCandidateReadID.size());
		
		for (int candidateIndex = 0; candidateIndex < splitAlignment.mCandidateReadID.size(); candidateIndex++)
		{
			ReadID readID;
			readID.id = splitAlignment.mCandidateReadID[candidateIndex];
			
			splitAlignment.mCandidateSeqIndex[candidateIndex] = sequenceIndices[readID.id];
		}
	}
}

bool SplitAlignment::Align(const Sequences& reference, const StringVec& sequences)
{
	if (mCandidateReadID.size() == 0)
	{
		return false;
	}
	
	string splitAlignSeq[2];
	for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
	{
		reference.Get(mAlignRefName[clusterEnd], mSplitSeqStrand[clusterEnd], mSplitAlignSeqStart[clusterEnd], mSplitAlignSeqLength[clusterEnd], splitAlignSeq[clusterEnd]);
	}
	
	SplitReadAligner splitReadAligner(matchScore, mismatchScore, gapScore, false, minAnchor * matchScore, splitAlignSeq[0], splitAlignSeq[1]);
	
	for (int candidateIndex = 0; candidateIndex < mCandidateReadID.size(); candidateIndex++)
	{
		ReadID readID;
		readID.id = mCandidateReadID[candidateIndex];
		int revComp = mCandidateRevComp[candidateIndex];
		int seqIndex = mCandidateSeqIndex[candidateIndex];
		string readSeq = sequences[seqIndex];
		
		if (revComp)
		{
			ReverseComplement(readSeq);
		}
		
		splitReadAligner.Align(readSeq);
		
		SplitReadAlignVec splitAlignments;
		splitReadAligner.GetAlignments(splitAlignments, (int)((float)readSeq.length() * (float)matchScore * 0.90), true, false, false);
		
		unordered_set<IntegerPair> readSplits;
		
		for (int splitAlignIndex = 0; splitAlignIndex < splitAlignments.size(); splitAlignIndex++)
		{
			const SplitReadAlignment& splitAlignment = splitAlignments[splitAlignIndex];
			
			if (readSplits.find(splitAlignment.refSplit) != readSplits.end())
			{
				continue;
			}
			
			readSplits.insert(splitAlignment.refSplit);
			
			IntegerPair breakPos;
			if (mSplitSeqStrand[0] == PlusStrand)
			{
				breakPos.first = mSplitAlignSeqStart[0] + splitAlignment.refSplit.first - 1;
			}
			else
			{
				breakPos.first = mSplitAlignSeqStart[0] + mSplitAlignSeqLength[0] - splitAlignment.refSplit.first;				
			}
			
			if (mSplitSeqStrand[1] == PlusStrand)
			{
				breakPos.second = mSplitAlignSeqStart[1] + splitAlignment.refSplit.second + 1;
			}
			else
			{
				breakPos.second = mSplitAlignSeqStart[1] + mSplitAlignSeqLength[1] - splitAlignment.refSplit.second - 2;
			}
			
			mAlignmentReadID.push_back(readID.id);
			mAlignmentBreakPos.push_back(breakPos);
			mAlignmentReadSplit.push_back(splitAlignment.readSplit);
			mAlignmentScore.push_back(splitAlignment.score1 + splitAlignment.score2);
		}
	}
	
	return true;
}

void SplitAlignment::WriteAlignments(ostream& out, SplitAlignmentMap& splitAlignments)
{
	for (SplitAlignmentMapConstIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		const SplitAlignment& splitAlignment = splitAlignIter->second;
		
		for (int alignmentIndex = 0; alignmentIndex < splitAlignment.mAlignmentReadID.size(); alignmentIndex++)
		{
			ReadID readID;
			readID.id = splitAlignment.mAlignmentReadID[alignmentIndex];
			
			out << id << "\t";
			out << readID.fragmentIndex << "\t" << readID.readEnd << "\t";
			out << splitAlignment.mAlignmentBreakPos[alignmentIndex].first << "\t" << splitAlignment.mAlignmentBreakPos[alignmentIndex].second << "\t";
			out << splitAlignment.mAlignmentReadSplit[alignmentIndex].first << "\t" << splitAlignment.mAlignmentReadSplit[alignmentIndex].second << "\t";
			out << splitAlignment.mAlignmentScore[alignmentIndex] << "\t";
			out << endl;
		}
	}
}

void SplitAlignment::CalculateBreakRegion(int minReadLength, int maxReadLength, int maxFragmentLength, int alignStart, int alignEnd, int strand, int& breakStart, int& breakLength)
{
	int alignRegionLength = alignEnd - alignStart + 1;
	
	// Push back the break region maximum read length nucleotides from 3' end of reads
	// or half way into the reads, whichever is minimum
	int pushBreakRegion = min(maxReadLength, (int)(0.5 * alignRegionLength));
	
	breakLength = maxFragmentLength - alignRegionLength - minReadLength + 2 * pushBreakRegion;
	
	if (strand == PlusStrand)
	{
		breakStart = alignEnd - pushBreakRegion + 1;
	}
	else
	{
		breakStart = alignStart + pushBreakRegion - 1;
	}
}

void SplitAlignment::CalculateSplitMateRegion(int minReadLength, int maxReadLength, int minFragmentLength, int maxFragmentLength, int breakStart, int breakLength, int strand, int& mateRegionStart, int& mateRegionEnd)
{
	if (strand == PlusStrand)
	{
		int breakEnd = breakStart + breakLength - 1;

		mateRegionStart = breakStart - maxFragmentLength + minReadLength;
		mateRegionEnd = breakEnd - minFragmentLength + maxReadLength;
	}
	else
	{
		int breakEnd = breakStart - breakLength + 1;

		mateRegionStart = breakEnd + minFragmentLength - maxReadLength;
		mateRegionEnd = breakStart + maxFragmentLength - minReadLength;
	}
}


