/*
 *  SplitAlignment.h
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#ifndef SPLITALIGNMENT_H_
#define SPLITALIGNMENT_H_

#include "AlignmentStream.h"
#include "ExonRegions.h"
#include "Sequences.h"
#include "ReadStream.h"
#include "ReadIndex.h"
#include <boost/unordered_set.hpp>

#include <vector>

using namespace std;
using namespace boost;
 
class SplitAlignment
{
public:
	typedef unordered_map<int,SplitAlignment> SplitAlignmentMap;
	typedef unordered_map<int,SplitAlignment>::iterator SplitAlignmentMapIter;
	typedef unordered_map<int,SplitAlignment>::const_iterator SplitAlignmentMapConstIter;

	bool Initialize(const LocationVec& alignPair, double fragmentLengthMean, 
					double fragmentLengthStdDev, int maxReadLength, int minReadLength);
	
	static bool FindCandidates(AlignmentStream* alignments, SplitAlignmentMap& splitAlignments);
	
	static void ReadCandidateSequences(IReadStream* readStream, SplitAlignmentMap& splitAlignments, StringVec& sequences);
	static void ReadCandidateSequences(const ReadIndex& readIndex, SplitAlignmentMap& splitAlignments, StringVec& sequences);

	bool Align(const Sequences& reference, const StringVec& sequences);

	static void WriteAlignments(ostream& out, SplitAlignmentMap& splitAlignments);
	
private:	
	inline void CalculateBreakRegion(int minReadLength, int maxReadLength, int maxFragmentLength, int alignStart, 
									 int alignEnd, int strand, int& breakStart, int& breakLength);
	inline void CalculateSplitMateRegion(int minReadLength, int maxReadLength, int minFragmentLength, 
										 int maxFragmentLength, int breakRegionStart, int breakRegionLength, 
										 int strand, int& mateRegionStart, int& mateRegionEnd);

	string mAlignRefName[2];
	int mAlignStrand[2];
	int mSplitAlignSeqStart[2];
	int mSplitAlignSeqLength[2];
	int mSplitSeqStrand[2];
	
	Location mMateRegions[2];
	
	IntegerVec mCandidateReadID;
	IntegerVec mCandidateRevComp;
	IntegerVec mCandidateSeqIndex;
	
	IntegerVec mAlignmentReadID;
	IntegerPairVec mAlignmentBreakPos;
	IntegerPairVec mAlignmentReadSplit;
	IntegerVec mAlignmentScore;
};

#endif

