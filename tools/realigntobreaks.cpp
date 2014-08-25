/*
 *  realign.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Indexer.h"
#include "ReadStream.h"
#include "AlignmentStream.h"
#include "Sequences.h"
#include "SimpleAligner.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


void Execute(int matchScore, int misMatchScore, int gapScore, int maxPaired, const string& referenceFasta, const string& readSeqsFilename, const string& alignmentsFilename, const string& splitsFilename, const string& breakpointsFilename);

void SystemTests()
{
	ofstream testbreaks("testbreaks.txt");
	testbreaks << "0\t0\t1\t+\t140\t180" << endl;
	testbreaks << "0\t1\t1\t-\t5821\t5880" << endl;
	testbreaks.close();
	
	ofstream testalign("testalignments.txt");
	testalign << "0\t0\t1\t+\t61\t120\t0\t0\t0" << endl;
	testalign << "1\t1\t1\t-\t5881\t5940\t0\t0\t0" << endl;
	testalign.close();
	
	ofstream testreads("testreads.fq");
	testreads << "@0/1" << endl << "TTTCATTAGTGCTTGAAGTTCGCAACACTGAATTCCTGTTTCCTAAGCTTATCTAAATGC" << endl << "+" << endl << string(60,'N') << endl;
	testreads << "@0/2" << endl << "TATGGATTAGTGAAACCCTACATAGAGTAAAATCCAGTAGGAACTGTAACTGTTGTTCAA" << endl << "+" << endl << string(60,'N') << endl;
	testreads << "@1/1" << endl << "TTGAACAACAGTTACAGTTCCTACTGGATTTTACTCTATGTAGGGTTTCACTAATCCATA" << endl << "+" << endl << string(60,'N') << endl;
	testreads << "@1/2" << endl << "GTAACCTTGGGCAAGTTATTCAACCTCTCTAAATCATGATTTATTCCTAGAAAAGGGAAG" << endl << "+" << endl << string(60,'N') << endl;
	testreads.close();
	
	ofstream testgenome("testgenome.fa");
	testgenome << ">1" << endl;
	testgenome << "GTATTCTCCTGGATTGAAACATTAAGACTTAGAAAACTCCGTAGTATCTAATTTGGTATA" << endl;
	testgenome << "TTTCATTAGTGCTTGAAGTTCGCAACACTGAATTCCTGTTTCCTAAGCTTATCTAAATGC" << endl;
	testgenome << "TAAATGAAATGTATGAAATATCTTTCTGTTGAACAACAGTTACAGTTCCTACTGGATTTT" << endl;
	testgenome << "CTTTTGGTTGGAAGACTGAAAGAGTGCTCAGAGAGAAGTCAGAGAGATCCAGTTCACATG" << endl;
	testgenome << "CTAATACATGATTTATTATTCTTGATTTATTATTCTGTGAAATGTGAAAAAAAATCCTTC" << endl;
	testgenome << "CCTAGATGGTGTCACTGATTGTTTTGTTTCCAAAGTATTAGTGATAAGAAAAATATTATG" << endl;
	for (int i = 0; i < 100 - 12; i++)
	{
		testgenome << string(60,'N') << endl;
	}
	testgenome << "GAGGTGTCTCTCTTCAGTCCTCCCGTAGTCTGAGGGAAACTCAACAGGGTACTGTCACTG" << endl;
	testgenome << "GCTGTTAGTCTCACCAGCTGAGAGAATGAGTACAGTGTCGAAGAGGGACCTGGCCAGTGC" << endl;
	testgenome << "TCCACAGCATCTACTACACCATTATTGAATGACTACATGATGTTGGCTATTGCATAAGTC" << endl;
	testgenome << "ACTCTATGTAGGGTTTCACTAATCCATAGAGTAAACACACACAAAAAAAGGTGCTGATGT" << endl;
	testgenome << "CTTCCCTTTTCTAGGAATAAATCATGATTTAGAGAGGTTGAATAACTTGCCCAAGGTTAC" << endl;
	testgenome << "ATAAAAATGGCAAGTCTGTCTGGCTAAAAAGCTCATTTTTTTTCATAAAGCTGCTCTTTA" << endl;
	testgenome.close();
	
	Execute(2, -3, -4, 100, "testgenome.fa", "testreads.fq", "testalignments.txt", "testsplits.txt", "testbreaks.txt");
	
	ifstream testsplits("testsplits.txt");
	string line;
	getline(testsplits, line); DebugCheck(line == "0\t1\t1\t+\t180\t1\t-\t5821\t\t32\t28\t64\t56\t120");
	getline(testsplits, line); DebugCheck(line == "1\t0\t1\t-\t5821\t1\t+\t180\t\t28\t32\t56\t64\t120");
}

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
		
		int start[2];
		int end[2];
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

struct AlignInfo
{
	int refStart;
	int strand;
	IntegerVec seqScores;
	IntegerVec refLengths;
	
	int BestPartialSeqLength() const
	{
		return max_element(seqScores.begin(), seqScores.end()) - seqScores.begin();
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
};

AlignInfo AlignSelfFullSSE(SimpleAligner& aligner, const RawAlignment& alignment, const Sequences& references, const PreppedReads& reads)
{
	static vector<short int> scoresBuffer(1024);
	static vector<short int> lengthsBuffer(1024);
	
	short int* scoresPtr = &scoresBuffer.front();
	short int* lengthsPtr = &lengthsBuffer.front();
	
	AlignInfo alignInfo;
	
	alignInfo.strand = alignment.strand;
	
	if (alignment.strand == PlusStrand)
	{
		alignInfo.refStart = alignment.region.start;
		const char* refPtr = references.Get(alignment.reference, alignment.region.start);
		aligner.AlignBandedSSE2BW7ScoreFwd(refPtr, reads.StartPtr(alignment.readEnd, PlusStrand), reads.EndPtr(alignment.readEnd, PlusStrand), scoresPtr, lengthsPtr);
	}
	else
	{
		alignInfo.refStart = alignment.region.end;
		const char* refPtr = references.Get(alignment.reference, alignment.region.end);
		aligner.AlignBandedSSE2BW7ScoreRev(refPtr, reads.StartPtr(alignment.readEnd, MinusStrand), reads.EndPtr(alignment.readEnd, MinusStrand), scoresPtr, lengthsPtr);
	}
	
	alignInfo.seqScores = IntegerVec(scoresPtr, scoresPtr + reads.ReadLength(alignment.readEnd) + 1);
	alignInfo.refLengths = IntegerVec(lengthsPtr, lengthsPtr + reads.ReadLength(alignment.readEnd) + 1);
	
	return alignInfo;
}

void Align3PrimeSeed16SSE(SimpleAligner& aligner, int readEnd, const string& refName, int strand, int start, const Sequences& references, const PreppedReads& reads, int searchLength, int& score, int& refPosition)
{
	if (strand == PlusStrand)
	{
		const char* refPtr = references.Get(refName, start);
		
		IntegerVec refLengthScores;
		aligner.AlignEndToEndRev(refPtr + searchLength - 1, searchLength, reads.StartPtr3PrimeSeed16(readEnd, MinusStrand), reads.StartPtr3PrimeSeed16(readEnd, MinusStrand) + 16, refLengthScores);
		int refLength = max_element(refLengthScores.begin(), refLengthScores.end()) - refLengthScores.begin();
		
		score = refLengthScores[refLength];
		refPosition = start + searchLength - 1 - refLength + 1;
	}
	else
	{
		const char* refPtr = references.Get(refName, start);
		
		IntegerVec refLengthScores;
		aligner.AlignEndToEndFwd(refPtr - searchLength + 1, searchLength, reads.StartPtr3PrimeSeed16(readEnd, PlusStrand), reads.StartPtr3PrimeSeed16(readEnd, PlusStrand) + 16, refLengthScores);
		int refLength = max_element(refLengthScores.begin(), refLengthScores.end()) - refLengthScores.begin();
		
		score = refLengthScores[refLength];
		refPosition = start - searchLength + 1 + refLength - 1;
	}
}

AlignInfo Align3PrimeFullSSE(SimpleAligner& aligner, int readEnd, const string& refName, int strand, const Sequences& references, const PreppedReads& reads, int refPosition)
{
	AlignInfo alignInfo;
	
	const char* refPtr = references.Get(refName, refPosition);
	
	static vector<short int> scoresBuffer(1024);
	static vector<short int> lengthsBuffer(1024);
	
	short int* scoresPtr = &scoresBuffer.front();
	short int* lengthsPtr = &lengthsBuffer.front();
	
	if (strand == PlusStrand)
	{
		aligner.AlignBandedSSE2BW7ScoreFwd(refPtr, reads.StartPtr(readEnd, MinusStrand), reads.EndPtr(readEnd, MinusStrand), scoresPtr, lengthsPtr);
	}
	else
	{
		aligner.AlignBandedSSE2BW7ScoreRev(refPtr, reads.StartPtr(readEnd, PlusStrand), reads.EndPtr(readEnd, PlusStrand), scoresPtr, lengthsPtr);
	}
	
	alignInfo.strand = strand;
	alignInfo.refStart = refPosition;
	alignInfo.seqScores = IntegerVec(scoresPtr, scoresPtr + reads.ReadLength(readEnd) + 1);
	alignInfo.refLengths = IntegerVec(lengthsPtr, lengthsPtr + reads.ReadLength(readEnd) + 1);
	
	return alignInfo;
}

void Align5PrimeSeed16SSE(SimpleAligner& aligner, int readEnd, const string& refName, int strand, int start, const Sequences& references, const PreppedReads& reads, int searchLength, int& score, int& refPosition)
{
	if (strand == PlusStrand)
	{
		const char* refPtr = references.Get(refName, start);
		
		IntegerVec refLengthScores;
		aligner.AlignEndToEndRev(refPtr + searchLength - 1, searchLength, reads.StartPtr5PrimeSeed16(readEnd, PlusStrand), reads.StartPtr5PrimeSeed16(readEnd, PlusStrand) + 16, refLengthScores);
		int refLength = max_element(refLengthScores.begin(), refLengthScores.end()) - refLengthScores.begin();
		
		score = refLengthScores[refLength];
		refPosition = start + searchLength - 1 - refLength + 1;
	}
	else
	{
		const char* refPtr = references.Get(refName, start);
		
		IntegerVec refLengthScores;
		aligner.AlignEndToEndFwd(refPtr - searchLength + 1, searchLength, reads.StartPtr5PrimeSeed16(readEnd, MinusStrand), reads.StartPtr5PrimeSeed16(readEnd, MinusStrand) + 16, refLengthScores);
		int refLength = max_element(refLengthScores.begin(), refLengthScores.end()) - refLengthScores.begin();
		
		score = refLengthScores[refLength];
		refPosition = start - searchLength + 1 + refLength - 1;
	}
}

AlignInfo Align5PrimeFullSSE(SimpleAligner& aligner, int readEnd, const string& refName, int strand, const Sequences& references, const PreppedReads& reads, int refPosition)
{
	AlignInfo alignInfo;
	
	const char* refPtr = references.Get(refName, refPosition);
	
	static vector<short int> scoresBuffer(1024);
	static vector<short int> lengthsBuffer(1024);
	
	short int* scoresPtr = &scoresBuffer.front();
	short int* lengthsPtr = &lengthsBuffer.front();
	
	if (strand == PlusStrand)
	{
		aligner.AlignBandedSSE2BW7ScoreFwd(refPtr, reads.StartPtr(readEnd, PlusStrand), reads.EndPtr(readEnd, PlusStrand), scoresPtr, lengthsPtr);
	}
	else
	{
		aligner.AlignBandedSSE2BW7ScoreRev(refPtr, reads.StartPtr(readEnd, MinusStrand), reads.EndPtr(readEnd, MinusStrand), scoresPtr, lengthsPtr);
	}
	
	alignInfo.strand = strand;
	alignInfo.refStart = refPosition;
	alignInfo.seqScores = IntegerVec(scoresPtr, scoresPtr + reads.ReadLength(readEnd) + 1);
	alignInfo.refLengths = IntegerVec(lengthsPtr, lengthsPtr + reads.ReadLength(readEnd) + 1);
	
	return alignInfo;
}

int BestSplitAlignment(const IntegerVec& scores1Fwd, const IntegerVec& scores2Rev, int breakInsertScore, IntegerVec& seq1Lengths, IntegerVec& seq2Lengths)
{
	IntegerVec scores1FwdBreakInsert(scores1Fwd.size());
	IntegerVec scores1FwdPrevMax(scores1Fwd.size());
	
	scores1FwdBreakInsert[0] = scores1Fwd[0];
	scores1FwdPrevMax[0] = 0;
	for (int idx = 1; idx < scores1Fwd.size(); idx++)
	{
		if (scores1FwdBreakInsert[idx-1] + breakInsertScore >= scores1Fwd[idx])
		{
			scores1FwdBreakInsert[idx] = scores1FwdBreakInsert[idx-1] + breakInsertScore;
			scores1FwdPrevMax[idx] = scores1FwdPrevMax[idx-1];
		}
		else
		{
			scores1FwdBreakInsert[idx] = scores1Fwd[idx];
			scores1FwdPrevMax[idx] = idx;
		}
	}
	
	IntegerVec seq1LengthMax;
	
	int score = numeric_limits<int>::min();
	for (int seq1Length = 0; seq1Length < scores1Fwd.size(); seq1Length++)
	{
		int seq2Length = scores2Rev.size() - seq1Length - 1;
		
		if (scores1FwdBreakInsert[seq1Length] + scores2Rev[seq2Length] > score)
		{
			score = scores1FwdBreakInsert[seq1Length] + scores2Rev[seq2Length];
			seq1LengthMax.clear();
		}
		
		if (scores1FwdBreakInsert[seq1Length] + scores2Rev[seq2Length] >= score)
		{
			seq1LengthMax.push_back(seq1Length);
		}
	}
	
	for (int idx = 0; idx < seq1LengthMax.size(); idx++)
	{
		int seq1Length = seq1LengthMax[idx];
		int seq2Length = scores2Rev.size() - seq1Length - 1;
		
		if (scores1FwdBreakInsert[seq1Length] == scores1Fwd[seq1Length])
		{
			seq1Lengths.push_back(seq1Length);
			seq2Lengths.push_back(seq2Length);
		}
		
		if (scores1FwdPrevMax[seq1Length] != seq1Length)
		{
			seq1Lengths.push_back(scores1FwdPrevMax[seq1Length]);
			seq2Lengths.push_back(seq2Length);
		}
	}
	
	return score;
}

class BinnedLocations
{
public:
	BinnedLocations(int binSpacing) : mBinSpacing(binSpacing) {}
	
	void AddConnection(const Location& location1, const Location& location2)
	{
		AddLocation(location1);
		AddLocation(location2);
	}
	
	LocationVec GetConnected(const string& reference, int strand, int start, int end)
	{
		LocationVec connected;
		
		unordered_map<string,unordered_map<int,IntegerVec> >::const_iterator findRefIter = mBinned[strand].find(reference);
		if (findRefIter != mBinned[strand].end())
		{
			int startBin = start / mBinSpacing;
			int endBin = end / mBinSpacing;
			
			for (int bin = startBin; bin <= endBin; bin++)
			{
				unordered_map<int,IntegerVec>::const_iterator findBinIter = findRefIter->second.find(bin);
				if (findBinIter != findRefIter->second.end())
				{
					for (IntegerVecConstIter iter = findBinIter->second.begin(); iter != findBinIter->second.end(); iter++)
					{
						const Location& location = mLocations[*iter];
						
						if (location.start <= end && location.end >= start)
						{
							connected.push_back(mLocations[ConnectedLocationID(*iter)]);
						}
					}
				}
			}
		}
		
		return connected;
	}
	
private:
	void AddLocation(const Location& location)
	{
		int id = mLocations.size();
		
		mLocations.push_back(location);
		
		int startBin = location.start / mBinSpacing;
		int endBin = location.end / mBinSpacing;
		
		for (int bin = startBin; bin <= endBin; bin++)
		{
			mBinned[location.strand][location.refName][bin].push_back(id);
		}
	}
	
	int ConnectedLocationID(int id)
	{
		return id ^ 1;
	}
	
	int mBinSpacing;
	vector<Location> mLocations;
	unordered_map<string,unordered_map<int,IntegerVec> > mBinned[2];
};

int main(int argc, char* argv[])
{
	int matchScore;
	int misMatchScore;
	int gapScore;
	int maxPaired;
	string referenceFasta;
	string readSeqsFilename;
	string alignmentsFilename;
	string splitsFilename;
	string breakpointsFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Realignment Tool");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> maxPairedArg("p","pairedmax","Max Paired Alignments",true,0,"int",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqsFilenameArg("s","seqs","Read Sequences Fastq",true,"","string",cmd);
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Compact Alignments",true,"","string",cmd);
		TCLAP::ValueArg<string> splitsFilenameArg("t","split","Output Split Alignments",true,"","string",cmd);
		TCLAP::ValueArg<string> breakpointsFilenameArg("b","breakpoints","Breakpoints Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		maxPaired = maxPairedArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		readSeqsFilename = readSeqsFilenameArg.getValue();
		alignmentsFilename = alignmentsFilenameArg.getValue();
		splitsFilename = splitsFilenameArg.getValue();
		breakpointsFilename = breakpointsFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	Execute(matchScore, misMatchScore, gapScore, maxPaired, referenceFasta, readSeqsFilename, alignmentsFilename, splitsFilename, breakpointsFilename);
}

void Execute(int matchScore, int misMatchScore, int gapScore, int maxPaired, const string& referenceFasta, const string& readSeqsFilename, const string& alignmentsFilename, const string& splitsFilename, const string& breakpointsFilename)
{
	int cMateSearchLength = 1000;
	int cSeedScoreThreshold = 24;
	int cMaxBreakpointHomology = 10;
	double cMatePercentScoreThreshold = 0.8;
	int cMinAnchor = 8;
	
	cerr << "Reading breakpoints" << endl;
	
	BinnedLocations breakpoints(1000);
	{
		LocationVecMap breakpointLocations;
		ReadAlignRegionPairs(breakpointsFilename, breakpointLocations);
		
		for (LocationVecMapConstIter breakIter = breakpointLocations.begin(); breakIter != breakpointLocations.end(); breakIter++)
		{
			breakpoints.AddConnection(breakIter->second[0], breakIter->second[1]);
		}
	}
	
	cerr << "Reading reference fasta" << endl;
	
	Sequences referenceSequences(1000);
	referenceSequences.Read(referenceFasta);
	
	cerr << "Reading fastq sequences" << endl;

	ifstream readSeqsFile(readSeqsFilename.c_str());
	CheckFile(readSeqsFile, readSeqsFilename);
	
	FastqReadStream readSeqsStream(readSeqsFile);
	
	PreppedReads preppedReads;
	preppedReads.Prep(readSeqsStream);
	
	cerr << "Realigning" << endl;
	
	ofstream splitsFile(splitsFilename.c_str());
	CheckFile(splitsFile, splitsFilename);
	
	SimpleAligner aligner(matchScore, misMatchScore, gapScore);
	
	CompactAlignmentStream alignmentStream(alignmentsFilename);
	FragmentAlignmentStream fragmentAlignmentStream(&alignmentStream);
	
	RawAlignmentVec alignments;
	while (fragmentAlignmentStream.GetNextAlignments(alignments))
	{
		const string& fragment = alignments.front().fragment;
		
		preppedReads.SetCurrentRead(SAFEPARSE(int, fragment));
		
		IntegerVec alignmentIndices[2];
		vector<AlignInfo> selfAlignInfos[2];
		
		LocationVec mateAlignLocations[2];
		vector<AlignInfo> mateAlignInfos[2];
		
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			int mateEnd = OtherReadEnd(alignment.readEnd);
			
			int selfAlignStart = (alignment.strand == PlusStrand) ? alignment.region.start : alignment.region.end;
			
			int seedScore;
			int seedPosition;
			Align3PrimeSeed16SSE(aligner, mateEnd, alignment.reference, alignment.strand, selfAlignStart, referenceSequences, preppedReads, cMateSearchLength, seedScore, seedPosition);
			
			if (seedScore < cSeedScoreThreshold)
			{
				continue;
			}
			
			alignmentIndices[mateEnd].push_back(alignmentIndex);
			selfAlignInfos[mateEnd].push_back(Align3PrimeFullSSE(aligner, mateEnd, alignment.reference, alignment.strand, referenceSequences, preppedReads, seedPosition));
			
			LocationVec connected = breakpoints.GetConnected(alignment.reference, alignment.strand, seedPosition - 16, seedPosition + 16);
			
			for (LocationVecConstIter connIter = connected.begin(); connIter != connected.end(); connIter++)
			{
				int mateAlignStart = (connIter->strand == PlusStrand) ? connIter->start : connIter->end;
				
				int seedScore;
				int seedPosition;
				Align5PrimeSeed16SSE(aligner, mateEnd, connIter->refName, connIter->strand, mateAlignStart, referenceSequences, preppedReads, cMateSearchLength, seedScore, seedPosition);
				
				if (seedScore >= cSeedScoreThreshold)
				{
					mateAlignLocations[mateEnd].push_back(*connIter);
					mateAlignInfos[mateEnd].push_back(Align5PrimeFullSSE(aligner, mateEnd, connIter->refName, connIter->strand, referenceSequences, preppedReads, seedPosition));
				}
			}
		}
		
		for (int mateEnd = 0; mateEnd <= 1; mateEnd++)
		{
			for (int selfAlignIndex = 0; selfAlignIndex < alignmentIndices[mateEnd].size(); selfAlignIndex++)
			{
				const RawAlignment& alignment = alignments[alignmentIndices[mateEnd][selfAlignIndex]];
				const AlignInfo& selfAlignInfo = selfAlignInfos[mateEnd][selfAlignIndex];
				
				for (int mateAlignIndex = 0; mateAlignIndex < mateAlignLocations[mateEnd].size(); mateAlignIndex++)
				{
					const Location& mateLocation = mateAlignLocations[mateEnd][mateAlignIndex];
					const AlignInfo& mateAlignInfo = mateAlignInfos[mateEnd][mateAlignIndex];
					
					IntegerVec seq1Length;
					IntegerVec seq2Length;
					int score = BestSplitAlignment(selfAlignInfo.seqScores, mateAlignInfo.seqScores, -1, seq1Length, seq2Length);
					
					if (seq1Length.size() > cMaxBreakpointHomology)
					{
						continue;
					}
					
					string readSeq = preppedReads.Sequence(mateEnd);
					
					for (int i = 0; i < seq1Length.size(); i++)
					{
						if (seq1Length[i] < cMinAnchor || seq2Length[i] < cMinAnchor)
						{
							continue;
						}
						
						splitsFile << fragment << "\t";
						splitsFile << mateEnd << "\t";
						splitsFile << alignment.reference << "\t";
						splitsFile << ((alignment.strand == PlusStrand) ? "+" : "-") << "\t";
						splitsFile << selfAlignInfo.BreakPosition(seq1Length[i]) << "\t";
						splitsFile << mateLocation.refName << "\t";
						splitsFile << ((mateLocation.strand == PlusStrand) ? "+" : "-") << "\t";
						splitsFile << mateAlignInfo.BreakPosition(seq2Length[i]) << "\t";
						splitsFile << readSeq.substr(seq1Length[i], readSeq.size() - seq2Length[i] - seq1Length[i]) << "\t";
						splitsFile << seq1Length[i] << "\t";
						splitsFile << seq2Length[i] << "\t";
						splitsFile << selfAlignInfo.seqScores[seq1Length[i]] << "\t";
						splitsFile << mateAlignInfo.seqScores[seq2Length[i]] << "\t";
						splitsFile << score << endl;
					}
				}
			}
		}
	}
}




