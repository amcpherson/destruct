/*
 *  DiscordantAlignments.h
 */

#ifndef DISCORDANTALIGNMENTS_H_
#define DISCORDANTALIGNMENTS_H_

#include "DebugCheck.h"

#include <string>
#include <map>
#include <iostream>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;


inline bool operator==(const RefStrand& rs1, const RefStrand& rs2)
{
	return rs1.id == rs2.id;
}

inline size_t hash_value(const RefStrand& rs)
{
	size_t seed = 0;
	hash_combine(seed, rs.id);
	return seed;
}

struct AlignID
{
	union
	{
		struct
		{
			unsigned alignIndex : 31;
			unsigned readEnd : 1;
		};
		
		int id;
	};
};


class DiscordantAlignments
{
public:
	DiscordantAlignments()
	{}
	
	void AddFragmentAlignments(const CompAlignVec& alignments, float fragmentMean, float fragmentStdDev)
	{
		if (alignments.empty())
		{
			return;
		}
		
		int fragmentRemap = mFragmentIndices.size();
		
		mFragmentIndices.push_back(alignments.front().readID.fragmentIndex);
		mFragmentMeans.push_back(fragmentMean);
		mFragmentStdDevs.push_back(fragmentStdDev);
		
		int alignIndexBegin[2];
		for (int readEnd = 0; readEnd < 2; readEnd++)
		{
			alignIndexBegin[readEnd] = mRefStrands[readEnd].size();
		}
		
		for (CompAlignVecConstIter alignIter = alignments.begin(); alignIter != alignments.end(); alignIter++)
		{
			const CompactAlignment& alignment = *alignIter;
			
			DebugCheck(mFragmentIndices.back() == alignment.readID.fragmentIndex);
			
			mFragmentRemaps[alignment.readID.readEnd].push_back(fragmentRemap);
			
			mRefStrands[alignment.readID.readEnd].push_back(alignment.refStrand);
			mStarts[alignment.readID.readEnd].push_back(alignment.region.start);
			mLengths[alignment.readID.readEnd].push_back(alignment.region.end - alignment.region.start + 1);
			mAlignProb[alignment.readID.readEnd].push_back((uint16_t)((float)numeric_limits<uint16_t>::max() * alignment.alignProb));
			mChimericProb[alignment.readID.readEnd].push_back((uint16_t)((float)numeric_limits<uint16_t>::max() * alignment.chimericProb));
			mValidProb[alignment.readID.readEnd].push_back((uint16_t)((float)numeric_limits<uint16_t>::max() * alignment.validProb));
		}
		
		int alignIndexEnd[2];
		for (int readEnd = 0; readEnd < 2; readEnd++)
		{
			alignIndexEnd[readEnd] = mRefStrands[readEnd].size();
		}
		
		for (int readEnd = 0; readEnd < 2; readEnd++)
		{
			mMateBegin[readEnd].push_back(alignIndexBegin[OtherReadEnd(readEnd)]);
			mMateEnd[readEnd].push_back(alignIndexEnd[OtherReadEnd(readEnd)]);
			
			for (int alignIndex = alignIndexBegin[readEnd]; alignIndex < alignIndexEnd[readEnd]; alignIndex++)
			{
				AlignID alignID;
				
				alignID.alignIndex = alignIndex;
				alignID.readEnd = readEnd;
				
				const RefStrand& refStrand = mRefStrands[readEnd][alignIndex];
				
				mBinnedReads[refStrand].push_back(alignID);
			}
		}
	}
	
	void StartRefStrandIteration() const
	{
		mFinishedSectorIteration = DoSectorIteration(false);
	}
	
	void NextRefStrandIteration() const
	{
		if (!mFinishedSectorIteration)
		{
			mFinishedSectorIteration = DoSectorIteration(true);
		}
	}
	
	void RetrieveRefStrandAlignments(CompAlignVec& alignments1, CompAlignVec& alignments2) const
	{
		const vector<pair<AlignID,AlignID> >& binnedMates = mRefStrandIter2->second;
		
		for (vector<pair<AlignID,AlignID> >::const_iterator mateIter = binnedMates.begin(); mateIter != binnedMates.end(); mateIter++)
		{
			AlignID alignID1 = mateIter->first;
			AlignID alignID2 = mateIter->second;
			
			CompactAlignment alignment1;
			CompactAlignment alignment2;
			
			DebugCheck(alignID1.readEnd != alignID2.readEnd);
			
			int fragmentRemap = mFragmentRemaps[alignID1.readEnd][alignID1.alignIndex];
			int fragmentIndex = mFragmentIndices[fragmentRemap];
			
			alignment1.readID.fragmentIndex = fragmentIndex;
			alignment1.readID.readEnd = alignID1.readEnd;
			alignment1.refStrand = mRefStrands[alignID1.readEnd][alignID1.alignIndex];
			alignment1.region.start = mStarts[alignID1.readEnd][alignID1.alignIndex];
			alignment1.region.end = alignment1.region.start + mLengths[alignID1.readEnd][alignID1.alignIndex] - 1;
			alignment1.alignProb = (float)mAlignProb[alignID1.readEnd][alignID1.alignIndex] / (float)numeric_limits<uint16_t>::max();
			alignment1.chimericProb = (float)mChimericProb[alignID1.readEnd][alignID1.alignIndex] / (float)numeric_limits<uint16_t>::max();
			alignment1.validProb = (float)mValidProb[alignID1.readEnd][alignID1.alignIndex] / (float)numeric_limits<uint16_t>::max();
			
			alignment2.readID.fragmentIndex = fragmentIndex;
			alignment2.readID.readEnd = alignID2.readEnd;
			alignment2.refStrand = mRefStrands[alignID2.readEnd][alignID2.alignIndex];
			alignment2.region.start = mStarts[alignID2.readEnd][alignID2.alignIndex];
			alignment2.region.end = alignment2.region.start + mLengths[alignID2.readEnd][alignID2.alignIndex] - 1;
			alignment2.alignProb = (float)mAlignProb[alignID2.readEnd][alignID2.alignIndex] / (float)numeric_limits<uint16_t>::max();
			alignment2.chimericProb = (float)mChimericProb[alignID2.readEnd][alignID2.alignIndex] / (float)numeric_limits<uint16_t>::max();
			alignment2.validProb = (float)mValidProb[alignID2.readEnd][alignID2.alignIndex] / (float)numeric_limits<uint16_t>::max();
			
			alignments1.push_back(alignment1);
			alignments2.push_back(alignment2);
		}
	}
	
	bool FinishedRefStrandIteration() const
	{
		return mFinishedSectorIteration;	
	}
	
private:
	void GenerateBinnedMates() const
	{
		mBinnedMates.clear();
		
		const vector<AlignID>& alignIDs = mRefStrandIter1->second;
		
		for (vector<AlignID>::const_iterator alignIDIter = alignIDs.begin(); alignIDIter != alignIDs.end(); alignIDIter++)
		{
			AlignID alignID = *alignIDIter;
			
			int fragmentRemap = mFragmentRemaps[alignID.readEnd][alignID.alignIndex];
			RefStrand refStrand = mRefStrands[alignID.readEnd][alignID.alignIndex];
			int start = mStarts[alignID.readEnd][alignID.alignIndex];
			
			for (int mateAlignIndex = mMateBegin[alignID.readEnd][fragmentRemap]; mateAlignIndex < mMateEnd[alignID.readEnd][fragmentRemap]; mateAlignIndex++)
			{
				AlignID mateAlignID;
				mateAlignID.alignIndex = mateAlignIndex;
				mateAlignID.readEnd = OtherReadEnd(alignID.readEnd);
				
				RefStrand mateRefStrand = mRefStrands[mateAlignID.readEnd][mateAlignID.alignIndex];
				
				int mateStart = mStarts[mateAlignID.readEnd][mateAlignID.alignIndex];
				
				if (refStrand.id != mateRefStrand.id || start < mateStart)
				{
					mBinnedMates[mateRefStrand].push_back(make_pair(alignID,mateAlignID));
				}
			}
		}
	}
	
	bool DoSectorIteration(bool resume) const
	{
		if (mBinnedReads.empty())
		{
			return true;
		}
		
		if (resume)
		{
			goto resumeiteration;
		}
		
		for (mRefStrandIter1 = mBinnedReads.begin(); mRefStrandIter1 != mBinnedReads.end(); mRefStrandIter1++)
		{
			GenerateBinnedMates();
			
			for (mRefStrandIter2 = mBinnedMates.begin(); mRefStrandIter2 != mBinnedMates.end(); mRefStrandIter2++)
			{
				if (mRefStrandIter2->first.id > mRefStrandIter1->first.id)
				{
					continue;
				}
				
				return false;
				
				resumeiteration: (void)(0);
			}
		}
		
		return true;
	}	
	
	vector<int> mFragmentRemaps[2];
	vector<RefStrand> mRefStrands[2];
	vector<int> mStarts[2];
	vector<unsigned char> mLengths[2];
	vector<uint16_t> mAlignProb[2];
	vector<uint16_t> mChimericProb[2];
	vector<uint16_t> mValidProb[2];
	
	vector<int> mFragmentIndices;
	vector<float> mFragmentMeans;
	vector<float> mFragmentStdDevs;
	vector<int> mMateBegin[2];
	vector<int> mMateEnd[2];
	
	unordered_map<RefStrand,vector<AlignID> > mBinnedReads;
	
	mutable bool mFinishedSectorIteration;
	mutable unordered_map<RefStrand,vector<AlignID> >::const_iterator mRefStrandIter1;
	mutable unordered_map<RefStrand,vector<pair<AlignID,AlignID> > > mBinnedMates;
	mutable unordered_map<RefStrand,vector<pair<AlignID,AlignID> > >::const_iterator mRefStrandIter2;
};


#endif

