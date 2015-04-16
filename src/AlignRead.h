/*
 *  AlignRead.h
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Sequences.h"
#include "AlignmentStream.h"
#include "SimpleAligner.h"
#include "ReadStream.h"

#include <iostream>
#include <string>

using namespace boost;
using namespace std;


class PreppedReads
{
public:
	void Prep(FastqReadStream& readSeqsStream)
	{
		RawRead rawRead;
		while (readSeqsStream.GetNextRead(rawRead))
		{
			ReadID readID;
			readID.fragmentIndex = SAFEPARSE(int, rawRead.fragment);
			readID.readEnd = rawRead.readEnd;
			
			mReadSequences += string(16,'X');
			
			string readSeqPlus = rawRead.sequence;
			reverse(readSeqPlus.begin(), readSeqPlus.end());
			
			mReadSeqInfo[readID].start[PlusStrand] = mReadSequences.size();
			mReadSequences += readSeqPlus;
			mReadSeqInfo[readID].end[PlusStrand] = mReadSequences.size();
			
			mReadSequences += string(16,'X');
			
			string readSeqMinus = rawRead.sequence;
			ReverseComplement(readSeqMinus);
			reverse(readSeqMinus.begin(), readSeqMinus.end());
			
			mReadSeqInfo[readID].start[MinusStrand] = mReadSequences.size();
			mReadSequences += readSeqMinus;
			mReadSeqInfo[readID].end[MinusStrand] = mReadSequences.size();
		}
		
		mReadSequences += string(16,'X');
	}
	
	void SetCurrentRead(int fragmentIndex)
	{
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			ReadID readID;
			readID.fragmentIndex = fragmentIndex;
			readID.readEnd = readEnd;
			
			unordered_map<ReadID,ReadSeqInfo>::const_iterator infoIter = mReadSeqInfo.find(readID);
			
			if (infoIter == mReadSeqInfo.end())
			{
				cerr << "Error: Could not find sequence for read " << readID.fragmentIndex << " end " << readID.readEnd << endl;
				exit(1);
			}	
			
			for (int strand = 0; strand <= 1; strand++)
			{
				mCurrentSeqStartPtr[strand][readEnd] = &mReadSequences[infoIter->second.start[strand]];
				mCurrentSeqEndPtr[strand][readEnd] = &mReadSequences[infoIter->second.end[strand]];
			}
			
			mCurrentSeq5PrimeSeed16Ptr[PlusStrand][readEnd] = &mReadSequences[infoIter->second.end[PlusStrand] - 16];
			mCurrentSeq5PrimeSeed16Ptr[MinusStrand][readEnd] = &mReadSequences[infoIter->second.start[MinusStrand]];
			
			mCurrentSeq3PrimeSeed16Ptr[PlusStrand][readEnd] = &mReadSequences[infoIter->second.start[PlusStrand]];
			mCurrentSeq3PrimeSeed16Ptr[MinusStrand][readEnd] = &mReadSequences[infoIter->second.end[MinusStrand] - 16];
		}
	}
	
	const char* StartPtr(int readEnd, int strand) const
	{
		return mCurrentSeqStartPtr[strand][readEnd];
	}
	
	const char* EndPtr(int readEnd, int strand) const
	{
		return mCurrentSeqEndPtr[strand][readEnd];
	}
	
	const char* StartPtr5PrimeSeed16(int readEnd, int strand) const
	{
		return mCurrentSeq5PrimeSeed16Ptr[strand][readEnd];
	}
	
	const char* StartPtr3PrimeSeed16(int readEnd, int strand) const
	{
		return mCurrentSeq3PrimeSeed16Ptr[strand][readEnd];
	}
	
	int ReadLength(int readEnd) const
	{
		return mCurrentSeqEndPtr[0][readEnd] - mCurrentSeqStartPtr[0][readEnd];
	}
	
	string Sequence(int readEnd) const
	{
		string sequence(StartPtr(readEnd, PlusStrand), EndPtr(readEnd, PlusStrand));
		reverse(sequence.begin(), sequence.end());
		return sequence;
	}
	
private:
	struct ReadSeqInfo
	{
		ReadSeqInfo()
		{
			start[0] = 0;
			start[1] = 0;
			end[0] = 0;
			end[1] = 0;
		}
		
		size_t start[2];
		size_t end[2];
	};
	
	string mReadSequences;
	unordered_map<ReadID,ReadSeqInfo> mReadSeqInfo;
	const char* mCurrentSeqStartPtr[2][2];
	const char* mCurrentSeqEndPtr[2][2];
	const char* mCurrentSeq5PrimeSeed16Ptr[2][2];
	const char* mCurrentSeq3PrimeSeed16Ptr[2][2];
};

int AlignSelfScoreSSE(SimpleAligner& aligner, const RawAlignment& alignment, const Sequences& references, const PreppedReads& reads)
{
	int score;
	
	if (alignment.strand == PlusStrand)
	{
		const char* refPtr = references.Get(alignment.reference, alignment.region.start);
		score = aligner.AlignBandedSSE2BW7ScoreFwd(refPtr, reads.StartPtr(alignment.readEnd, PlusStrand), reads.EndPtr(alignment.readEnd, PlusStrand));
	}
	else
	{
		const char* refPtr = references.Get(alignment.reference, alignment.region.end);
		score = aligner.AlignBandedSSE2BW7ScoreRev(refPtr, reads.StartPtr(alignment.readEnd, MinusStrand), reads.EndPtr(alignment.readEnd, MinusStrand));
	}
	
	return score;
}

class AlignInfo
{
public:
	AlignInfo() : refStart(-1), strand(0) {}
	AlignInfo(int refStart, int strand, int readLength) : refStart(refStart), strand(strand), seqScores(readLength + 1 + 32), refLengths(readLength + 1 + 32) {}
	
	int BestPartialSeqLength() const
	{
		return max_element(seqScores.begin(), seqScores.end() - 32) - seqScores.begin();
	}

	int OuterPosition() const
	{
		return refStart;
	}
	
	int BreakPosition(int seqLength) const
	{
		if (strand == PlusStrand)
		{
			return refStart + refLengths[seqLength] - 1;
		}
		else
		{
			return refStart - refLengths[seqLength] + 1;
		}
	}
	
	int AlignmentStart(int seqLength) const
	{
		if (strand == PlusStrand)
		{
			return refStart;
		}
		else
		{
			return refStart - refLengths[seqLength] + 1;
		}
	}
	
	int AlignmentEnd(int seqLength) const
	{
		if (strand == PlusStrand)
		{
			return refStart + refLengths[seqLength] - 1;
		}
		else
		{
			return refStart;
		}
	}
	
	int AlignmentPosition(int seqLength) const
	{
		if (strand == PlusStrand)
		{
			return refStart + refLengths[seqLength] - 1;
		}
		else
		{
			return refStart - refLengths[seqLength] + 1;
		}
	}
	
	const short int* SeqScores() const
	{
		return &seqScores.front();
	}
	
	int SeqScoresLength() const
	{
		return seqScores.size() - 32;
	}
	
	const short int* RefLengths() const
	{
		return &refLengths.front();
	}
	
	int RefLengthsLength() const
	{
		return refLengths.size() - 32;
	}
	
private:
	friend AlignInfo AlignSelfFullSSE(SimpleAligner& aligner, const RawAlignment& alignment, const Sequences& references, const PreppedReads& reads);
	friend AlignInfo AlignFwdMateFullSSE(SimpleAligner& aligner, const RawAlignment& alignment, const Sequences& references, const PreppedReads& reads, int refPosition);
	friend AlignInfo AlignRevMateFullSSE(SimpleAligner& aligner, const RawAlignment& alignment, const Sequences& references, const PreppedReads& reads, int refPosition);
	
	int refStart;
	int strand;
	vector<short int> seqScores;
	vector<short int> refLengths;
};

AlignInfo AlignSelfFullSSE(SimpleAligner& aligner, const RawAlignment& alignment, const Sequences& references, const PreppedReads& reads)
{
	AlignInfo alignInfo;
	
	if (alignment.strand == PlusStrand)
	{
		alignInfo = AlignInfo(alignment.region.start, alignment.strand, reads.ReadLength(alignment.readEnd));
		
		const char* refPtr = references.Get(alignment.reference, alignment.region.start);
		aligner.AlignBandedSSE2BW7ScoreFwd(refPtr, reads.StartPtr(alignment.readEnd, PlusStrand), reads.EndPtr(alignment.readEnd, PlusStrand), &alignInfo.seqScores.front(), &alignInfo.refLengths.front());
	}
	else
	{
		alignInfo = AlignInfo(alignment.region.end, alignment.strand, reads.ReadLength(alignment.readEnd));
		
		const char* refPtr = references.Get(alignment.reference, alignment.region.end);
		aligner.AlignBandedSSE2BW7ScoreRev(refPtr, reads.StartPtr(alignment.readEnd, MinusStrand), reads.EndPtr(alignment.readEnd, MinusStrand), &alignInfo.seqScores.front(), &alignInfo.refLengths.front());
	}
	
	return alignInfo;
}

void AlignMate3PrimeSeed16SSE(SimpleAligner& aligner, const RawAlignment& alignment, const Sequences& references, const PreppedReads& reads, int searchLength, int& score, int& refPosition)
{
	int mateEnd = OtherReadEnd(alignment.readEnd);
	
	if (alignment.strand == PlusStrand)
	{
		const char* refPtr = references.Get(alignment.reference, alignment.region.start);
		
		vector<short int> refLengthScores(searchLength + 1 + 32);
		aligner.Align16baseSSE2Rev(refPtr, refPtr + searchLength, reads.StartPtr3PrimeSeed16(mateEnd, MinusStrand), &refLengthScores.front());
		
		int refLength = max_element(refLengthScores.begin(), refLengthScores.end() - 32) - refLengthScores.begin();
		
		score = refLengthScores[refLength];
		refPosition = alignment.region.start + searchLength - 1 - refLength + 1;
	}
	else
	{
		const char* refPtr = references.Get(alignment.reference, alignment.region.end);
		
		vector<short int> refLengthScores(searchLength + 1 + 32);
		aligner.Align16baseSSE2Fwd(refPtr - searchLength + 1, refPtr + 1, reads.StartPtr3PrimeSeed16(mateEnd, PlusStrand), &refLengthScores.front());
		
		int refLength = max_element(refLengthScores.begin(), refLengthScores.end() - 32) - refLengthScores.begin();
		
		score = refLengthScores[refLength];
		refPosition = alignment.region.end - searchLength + 1 + refLength - 1;
	}
}

void AlignMate3PrimeSeed16SSE(SimpleAligner& aligner, const RawAlignment& alignment, const Sequences& references, const PreppedReads& reads, int minFragmentLength, int maxFragmentLength, int& score, int& refPosition)
{
	int mateEnd = OtherReadEnd(alignment.readEnd);
	
	int searchLength = maxFragmentLength - minFragmentLength + 16;
	
	if (alignment.strand == PlusStrand)
	{
		int searchStart = alignment.region.start + minFragmentLength - reads.ReadLength(mateEnd);
		
		const char* refPtr = references.Get(alignment.reference, searchStart);
		
		vector<short int> refLengthScores(searchLength + 1 + 32);
		aligner.Align16baseSSE2Rev(refPtr, refPtr + searchLength, reads.StartPtr3PrimeSeed16(mateEnd, MinusStrand), &refLengthScores.front());
		
		int refLength = (int)(max_element(refLengthScores.begin(), refLengthScores.end() - 32) - refLengthScores.begin());
		
		score = refLengthScores[refLength];
		refPosition = searchStart + searchLength - 1 - refLength + 1;
	}
	else
	{
		int searchEnd = alignment.region.end - minFragmentLength + reads.ReadLength(mateEnd);
		
		const char* refPtr = references.Get(alignment.reference, searchEnd);
		
		vector<short int> refLengthScores(searchLength + 1 + 32);
		aligner.Align16baseSSE2Fwd(refPtr - searchLength + 1, refPtr + 1, reads.StartPtr3PrimeSeed16(mateEnd, PlusStrand), &refLengthScores.front());
		
		int refLength = (int)(max_element(refLengthScores.begin(), refLengthScores.end() - 32) - refLengthScores.begin());
		
		score = refLengthScores[refLength];
		refPosition = searchEnd - searchLength + 1 + refLength - 1;
	}
}

AlignInfo AlignFwdMateFullSSE(SimpleAligner& aligner, const RawAlignment& alignment, const Sequences& references, const PreppedReads& reads, int refPosition)
{
	int mateEnd = OtherReadEnd(alignment.readEnd);
	
	AlignInfo alignInfo(refPosition, alignment.strand, reads.ReadLength(mateEnd));
	
	const char* refPtr = references.Get(alignment.reference, refPosition);
	
	if (alignment.strand == PlusStrand)
	{
		aligner.AlignBandedSSE2BW7ScoreFwd(refPtr, reads.StartPtr(mateEnd, MinusStrand), reads.EndPtr(mateEnd, MinusStrand), &alignInfo.seqScores.front(), &alignInfo.refLengths.front());
	}
	else
	{
		aligner.AlignBandedSSE2BW7ScoreRev(refPtr, reads.StartPtr(mateEnd, PlusStrand), reads.EndPtr(mateEnd, PlusStrand), &alignInfo.seqScores.front(), &alignInfo.refLengths.front());
	}
	
	return alignInfo;
}

AlignInfo AlignRevMateFullSSE(SimpleAligner& aligner, const RawAlignment& alignment, const Sequences& references, const PreppedReads& reads, int refPosition)
{
	int mateEnd = OtherReadEnd(alignment.readEnd);
	
	AlignInfo alignInfo(refPosition, OtherStrand(alignment.strand), reads.ReadLength(mateEnd));
	
	const char* refPtr = references.Get(alignment.reference, refPosition);
	
	if (alignment.strand == PlusStrand)
	{
		aligner.AlignBandedSSE2BW7ScoreRev(refPtr, reads.StartPtr(mateEnd, MinusStrand), reads.EndPtr(mateEnd, MinusStrand), &alignInfo.seqScores.front(), &alignInfo.refLengths.front());
	}
	else
	{
		aligner.AlignBandedSSE2BW7ScoreFwd(refPtr, reads.StartPtr(mateEnd, PlusStrand), reads.EndPtr(mateEnd, PlusStrand), &alignInfo.seqScores.front(), &alignInfo.refLengths.front());
	}
	
	return alignInfo;
}



