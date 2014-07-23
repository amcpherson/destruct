/*
 *  DiscordantAlignments.h
 */

#ifndef DISCORDANTALIGNMENTS_H_
#define DISCORDANTALIGNMENTS_H_

#include "DebugCheck.h"
#include "Indexer.h"
#include "AlignmentRecord.h"

#include <string>
#include <map>
#include <iostream>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;


bool IsConcordant(const vector<SpanningAlignmentRecord>& alignments, int maxFragmentLength)
{
	for (vector<SpanningAlignmentRecord>::const_iterator alignmentIter1 = alignments.begin(); alignmentIter1 != alignments.end(); alignmentIter1++)
	{
		for (vector<SpanningAlignmentRecord>::const_iterator alignmentIter2 = alignmentIter1 + 1; alignmentIter2 != alignments.end(); alignmentIter2++)
		{
			if (alignmentIter1->readEnd == alignmentIter2->readEnd)
			{
				continue;
			}
			
			if (alignmentIter1->chromosome != alignmentIter2->chromosome)
			{
				continue;
			}
			
			if (alignmentIter1->strand == "+" && alignmentIter2->strand == "-")
			{
				int inferredLength = alignmentIter2->end - alignmentIter1->start + 1;
				
				if (inferredLength >= 0 && inferredLength < maxFragmentLength)
				{
					return true;
				}
			}
			
			if (alignmentIter2->strand == "+" && alignmentIter1->strand == "-")
			{
				int inferredLength = alignmentIter1->end - alignmentIter2->start + 1;
				
				if (inferredLength >= 0 && inferredLength < maxFragmentLength)
				{
					return true;
				}
			}
		}
	}
	
	return false;
}


class ChromosomeStrandIndex
{
public:
	uint32_t Index(const string& chromosome, const string& strand)
	{
		RefStrand refStrand;
		refStrand.referenceIndex = mRefNameIndex.Index(chromosome);
		refStrand.strand = InterpretStrand(strand);

		return refStrand.id;
	}

	string GetChromosome(uint32_t idx) const
	{
		RefStrand refStrand;
		refStrand.id = idx;

		return mRefNameIndex.Get(refStrand.referenceIndex);
	}

	string GetStrand(uint32_t idx) const
	{
		RefStrand refStrand;
		refStrand.id = idx;

		return ((refStrand.strand == PlusStrand) ? "+" : "-");
	}

private:
	NameIndex mRefNameIndex;
};


class DiscordantAlignments
{
public:
	DiscordantAlignments(const vector<double>& fragmentMeans, const vector<double>& fragmentStdDevs, int maxFragmentLength)
	: mFragmentMeans(fragmentMeans), mFragmentStdDevs(fragmentStdDevs), mMaxFragmentLength(maxFragmentLength), mFragmentCount(0)
	{}

	void SetIncludedChromosomePair(const string& chromosome1, const string& chromosome2)
	{
		mIncludedChromosomePair = pair<string,string>(chromosome1, chromosome2);
	}

	void SetExcludedChromosomePairs(const vector<string>& chromosomes)
	{
		mExcludedChromosomePairs = unordered_set<string>(chromosomes.begin(), chromosomes.end());
	}

	bool IsExcluded(const string& chromosome1, const string& chromosome2) const
	{
		if (!mIncludedChromosomePair.first.empty() && !mIncludedChromosomePair.second.empty())
		{
			if ((chromosome1 == mIncludedChromosomePair.first && chromosome2 == mIncludedChromosomePair.second) || 
				(chromosome2 == mIncludedChromosomePair.first && chromosome1 == mIncludedChromosomePair.second))
			{
				return false;
			}
			else
			{
				return true;
			}
		}
		else if (mExcludedChromosomePairs.find(chromosome1) != mExcludedChromosomePairs.end() &&
			     mExcludedChromosomePairs.find(chromosome2) != mExcludedChromosomePairs.end())
		{
			return true;
		}

		return false;
	}

	vector<pair<size_t,size_t> > CalculateAlignmentPairs(const vector<SpanningAlignmentRecord>& alignments) const
	{
		vector<pair<size_t,size_t> > alignmentPairs;

		vector<size_t> readEndIdxs[2];
		for (size_t idx = 0; idx < alignments.size(); idx++)
		{
			readEndIdxs[alignments[idx].readEnd].push_back(idx);
		}

		for (vector<size_t>::const_iterator idx1Iter = readEndIdxs[0].begin(); idx1Iter != readEndIdxs[0].end(); idx1Iter++)
		{
			size_t idx1 = *idx1Iter;

			for (vector<size_t>::const_iterator idx2Iter = readEndIdxs[1].begin(); idx2Iter != readEndIdxs[1].end(); idx2Iter++)
			{
				size_t idx2 = *idx2Iter;

				if (IsExcluded(alignments[idx1].chromosome, alignments[idx2].chromosome))
				{
					continue;
				}

				alignmentPairs.push_back(pair<size_t,size_t>(idx1, idx2));
			}
		}

		return alignmentPairs;
	}

	vector<SpanningAlignmentRecord> CalculateFilteredAlignments(const vector<SpanningAlignmentRecord>& alignments) const
	{
		vector<pair<size_t,size_t> > alignmentPairs = CalculateAlignmentPairs(alignments);

		vector<bool> isIncluded(alignments.size(), false);

		for (vector<pair<size_t,size_t> >::const_iterator pairIter = alignmentPairs.begin(); pairIter != alignmentPairs.end(); pairIter++)
		{
			size_t idx1 = pairIter->first;
			size_t idx2 = pairIter->second;

			isIncluded[idx1] = true;
			isIncluded[idx2] = true;
		}

		vector<SpanningAlignmentRecord> filteredAlignments;

		for (int idx = 0; idx < alignments.size(); idx++)
		{
			if (isIncluded[idx])
			{
				filteredAlignments.push_back(alignments[idx]);
			}
		}

		return filteredAlignments;
	}

	void AddFragmentAlignments(const vector<SpanningAlignmentRecord>& fullAlignments)
	{
		if (fullAlignments.empty())
		{
			return;
		}
		
		if (IsConcordant(fullAlignments, mMaxFragmentLength))
		{
			return;
		}

		vector<SpanningAlignmentRecord> alignments = CalculateFilteredAlignments(fullAlignments);

		if (alignments.empty())
		{
			return;
		}
		
		size_t idxOffset = mLibIDs.size();
		DebugCheck(mLibIDs.size() == mReadIDs.size());
		DebugCheck(mLibIDs.size() == mReadEnds.size());
		DebugCheck(mLibIDs.size() == mAlignIDs.size());
		DebugCheck(mLibIDs.size() == mPositions.size());

		vector<pair<size_t,size_t> > alignmentPairs = CalculateAlignmentPairs(alignments);

		for (vector<pair<size_t,size_t> >::const_iterator pairIter = alignmentPairs.begin(); pairIter != alignmentPairs.end(); pairIter++)
		{
			size_t idx1 = pairIter->first;
			size_t idx2 = pairIter->second;

			const SpanningAlignmentRecord& alignment1 = alignments[idx1];
			const SpanningAlignmentRecord& alignment2 = alignments[idx2];

			uint32_t chrStrIdx1 = mChrStrIndex.Index(alignment1.chromosome, alignment1.strand);
			uint32_t chrStrIdx2 = mChrStrIndex.Index(alignment2.chromosome, alignment2.strand);

			pair<uint32_t,uint32_t> chrStrIdxPair;
			pair<size_t,size_t> alignmentIdxPair;
			if ((chrStrIdx1 < chrStrIdx2) || 
				((chrStrIdx1 == chrStrIdx2) && alignment1.GetOuterPosition() < alignment2.GetOuterPosition()))
			{
				chrStrIdxPair = pair<uint32_t,uint32_t>(chrStrIdx1, chrStrIdx2);
				alignmentIdxPair = pair<size_t,size_t>(idx1 + idxOffset, idx2 + idxOffset);
			}
			else
			{
				chrStrIdxPair = pair<uint32_t,uint32_t>(chrStrIdx2, chrStrIdx1);
				alignmentIdxPair = pair<size_t,size_t>(idx2 + idxOffset, idx1 + idxOffset);
			}

			mPaired[chrStrIdxPair].push_back(alignmentIdxPair);
		}

		for (vector<SpanningAlignmentRecord>::const_iterator alignmentIter = alignments.begin(); alignmentIter != alignments.end(); alignmentIter++)
		{
			mLibIDs.push_back(alignmentIter->libID);
			mReadIDs.push_back(alignmentIter->readID);
			mReadEnds.push_back(alignmentIter->readEnd);
			mAlignIDs.push_back(alignmentIter->alignID);
			mPositions.push_back(alignmentIter->GetOuterPosition());
		}

		mFragmentCount++;
	}

	vector<pair<uint32_t,uint32_t> > GetChrStrIdxPairs()
	{
		vector<pair<uint32_t,uint32_t> > chrStrIdxPairs;

		for (unordered_map<pair<uint32_t,uint32_t>,vector<pair<size_t,size_t> > >::const_iterator pairIter = mPaired.begin(); pairIter != mPaired.end(); pairIter++)
		{
			chrStrIdxPairs.push_back(pairIter->first);
		}

		return chrStrIdxPairs;
	}

	pair<string,string> GetChromosomePair(const pair<uint32_t,uint32_t>& chrStrIdxPair)
	{
		return pair<string,string>(mChrStrIndex.GetChromosome(chrStrIdxPair.first),
		                           mChrStrIndex.GetChromosome(chrStrIdxPair.second));
	}

	pair<string,string> GetStrandPair(const pair<uint32_t,uint32_t>& chrStrIdxPair)
	{
		return pair<string,string>(mChrStrIndex.GetStrand(chrStrIdxPair.first),
		                           mChrStrIndex.GetStrand(chrStrIdxPair.second));
	}

	vector<MatePair> CreateMatePairs(const pair<uint32_t,uint32_t>& chrStrIdxPair) const
	{
		vector<MatePair> matePairs;

		const vector<pair<size_t,size_t> >& alignmentIdxPairs = mPaired.find(chrStrIdxPair)->second;

		for (vector<pair<size_t,size_t> >::const_iterator alignmentIdxPairIter = alignmentIdxPairs.begin(); alignmentIdxPairIter != alignmentIdxPairs.end(); alignmentIdxPairIter++)
		{
			size_t alignmentIdx1 = alignmentIdxPairIter->first;
			size_t alignmentIdx2 = alignmentIdxPairIter->second;

			string strand1 = mChrStrIndex.GetStrand(chrStrIdxPair.first);
			string strand2 = mChrStrIndex.GetStrand(chrStrIdxPair.second);

			DebugCheck(mLibIDs[alignmentIdx1] == mLibIDs[alignmentIdx2]);

			MatePair matePair;
			matePair.x = (strand1 == "+") ? mPositions[alignmentIdx1] : -mPositions[alignmentIdx1];
			matePair.y = (strand2 == "+") ? mPositions[alignmentIdx2] : -mPositions[alignmentIdx2];
			matePair.u = mFragmentMeans[mLibIDs[alignmentIdx1]];
			matePair.s = mFragmentStdDevs[mLibIDs[alignmentIdx1]];

			matePairs.push_back(matePair);
		}

		return matePairs;
	}

	vector<pair<ReadInfo,ReadInfo> > CreateReadInfos(const pair<uint32_t,uint32_t>& chrStrIdxPair) const
	{
		vector<pair<ReadInfo,ReadInfo> > readInfos;

		const vector<pair<size_t,size_t> >& alignmentIdxPairs = mPaired.find(chrStrIdxPair)->second;

		for (vector<pair<size_t,size_t> >::const_iterator alignmentIdxPairIter = alignmentIdxPairs.begin(); alignmentIdxPairIter != alignmentIdxPairs.end(); alignmentIdxPairIter++)
		{
			size_t alignmentIdx1 = alignmentIdxPairIter->first;
			size_t alignmentIdx2 = alignmentIdxPairIter->second;

			DebugCheck(mLibIDs[alignmentIdx1] == mLibIDs[alignmentIdx2]);
			DebugCheck(mReadIDs[alignmentIdx1] == mReadIDs[alignmentIdx2]);
			DebugCheck(mReadEnds[alignmentIdx1] != mReadEnds[alignmentIdx2]);

			ReadInfo readInfo1;
			readInfo1.libID = mLibIDs[alignmentIdx1];
			readInfo1.readID = mReadIDs[alignmentIdx1];
			readInfo1.readEnd = mReadEnds[alignmentIdx1];
			readInfo1.alignID = mAlignIDs[alignmentIdx1];

			ReadInfo readInfo2;
			readInfo2.libID = mLibIDs[alignmentIdx2];
			readInfo2.readID = mReadIDs[alignmentIdx2];
			readInfo2.readEnd = mReadEnds[alignmentIdx2];
			readInfo2.alignID = mAlignIDs[alignmentIdx2];

			readInfos.push_back(make_pair(readInfo1, readInfo2));
		}

		return readInfos;
	}

	int GetFragmentCount() const
	{
		return mFragmentCount;
	}

	int GetAlignmentCount() const
	{
		return (int)mLibIDs.size();
	}
	
private:
	int mMaxFragmentLength;
	const vector<double> mFragmentMeans;
	const vector<double> mFragmentStdDevs;
	int mFragmentCount;

	ChromosomeStrandIndex mChrStrIndex;

	vector<uint16_t> mLibIDs;
	vector<uint32_t> mReadIDs;
	vector<uint8_t> mReadEnds;
	vector<uint16_t> mAlignIDs;
	vector<uint32_t> mPositions;

	unordered_map<pair<uint32_t,uint32_t>,vector<pair<size_t,size_t> > > mPaired;

	pair<string,string> mIncludedChromosomePair;
	unordered_set<string> mExcludedChromosomePairs;
};


#endif

