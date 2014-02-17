/*
 *  realign2.cpp
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
#include "AlignmentProbability.h"
#include "AlignRead.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


int BestSplitAlignment(const short int* scores1Fwd, int scores1FwdLength, const short int* scores2Rev, int scores2RevLength, int breakInsertScore, IntegerVec& seq1Lengths, IntegerVec& seq2Lengths)
{
	IntegerVec scores1FwdBreakInsert(scores1FwdLength);
	IntegerVec scores1FwdPrevMax(scores1FwdLength);
	
	scores1FwdBreakInsert[0] = scores1Fwd[0];
	scores1FwdPrevMax[0] = 0;
	for (int idx = 1; idx < scores1FwdLength; idx++)
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
	for (int seq1Length = 0; seq1Length < scores1FwdLength; seq1Length++)
	{
		int seq2Length = scores2RevLength - seq1Length - 1;
		
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
		int seq2Length = scores2RevLength - seq1Length - 1;
		
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

int main(int argc, char* argv[])
{
	int matchScore;
	int misMatchScore;
	int gapScore;
	int minFragmentLength;
	int maxFragmentLength;
	string referenceFasta;
	string readSeqsFilename;
	string alignmentsFilename;
	string statsFilename;
	double validReadPrior;
	double validReadThreshold;
	double chimericPrior;
	double chimericThreshold;
	double alignmentThreshold;
	string spanningFilename;
	string splitFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Realignment Tool");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> minFragmentLengthArg("","flmin","Minimum Fragment Length",true,0,"int",cmd);
		TCLAP::ValueArg<int> maxFragmentLengthArg("","flmax","Maximum Fragment Length",true,0,"int",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqsFilenameArg("s","seqs","Read Sequences Fastq",true,"","string",cmd);
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Sam Alignments",true,"","string",cmd);
		TCLAP::ValueArg<string> statsFilenameArg("z","stats","Stats Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> validReadPriorArg("","pvalid","Prior Probability of Valid Read",true,0.95,"float",cmd);
		TCLAP::ValueArg<double> validReadThresholdArg("","tvalid","Valid Posterior Threshold",true,0.1,"float",cmd);
		TCLAP::ValueArg<double> chimericPriorArg("","pchimer","Prior Probility of Chimeric Read",true,0.05,"float",cmd);
		TCLAP::ValueArg<double> chimericThresholdArg("","tchimer","Chimeric Posterior Threshold",true,0.1,"float",cmd);
		TCLAP::ValueArg<double> alignmentThresholdArg("","talign","Alignment Posterior Threshold",true,0.1,"float",cmd);
		TCLAP::ValueArg<string> spanningFilenameArg("","span","Spanning Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> splitFilenameArg("","split","Splits Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		minFragmentLength = minFragmentLengthArg.getValue();
		maxFragmentLength = maxFragmentLengthArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		readSeqsFilename = readSeqsFilenameArg.getValue();
		alignmentsFilename = alignmentsFilenameArg.getValue();
		statsFilename = statsFilenameArg.getValue();
		validReadPrior = validReadPriorArg.getValue();
		validReadThreshold = validReadThresholdArg.getValue();
		chimericPrior = chimericPriorArg.getValue();
		chimericThreshold = chimericThresholdArg.getValue();
		alignmentThreshold = alignmentThresholdArg.getValue();
		spanningFilename = spanningFilenameArg.getValue();
		splitFilename = splitFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	ofstream spanningFile(spanningFilename.c_str());
	CheckFile(spanningFile, spanningFilename);
	
	ofstream splitFile(splitFilename.c_str());
	CheckFile(splitFile, splitFilename);
	
	const int cSeedScoreThreshold = 8;
	const int cMaxBreakpointHomology = 25;
	const int cMinAnchor = 8;
	const int cBreakEndAdjust = 5;
	
	cerr << "Reading alignment stats" << endl;
	
	AlignmentProbability alignProbability(matchScore);
	alignProbability.ReadDistributions(statsFilename);
	
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
	
	SimpleAligner aligner(matchScore, misMatchScore, gapScore);
	
	SamAlignmentStream alignmentStream(alignmentsFilename);
	FragmentAlignmentStream fragmentAlignmentStream(&alignmentStream);
	
	RawAlignmentVec alignments;
	while (fragmentAlignmentStream.GetNextAlignments(alignments))
	{
		preppedReads.SetCurrentRead(SAFEPARSE(int, alignments.front().fragment));
		
		// Check that both ends are mapped
		bool readEndMapped[2] = {false,false};
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			readEndMapped[alignment.readEnd] = true;
		}
		if (!readEndMapped[0] || !readEndMapped[1])
		{
			continue;
		}
		
		// Realignments
		pair<int,int> bestAlignment[2] = {pair<int,int>(0,0),pair<int,int>(0,0)};
		vector<AlignInfo> selfAlignments;
		unordered_map<int,AlignInfo> mateFwdAlignments;
		unordered_map<int,AlignInfo> mateRevAlignments;
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			AlignInfo selfAlignInfo = AlignSelfFullSSE(aligner, alignment, referenceSequences, preppedReads);
			
			int selfSeqLength = selfAlignInfo.BestPartialSeqLength();
			int selfScore = selfAlignInfo.SeqScores()[selfSeqLength];
			
			bestAlignment[alignment.readEnd] = max(bestAlignment[alignment.readEnd], pair<int,int>(selfScore, selfSeqLength));
			
			selfAlignments.push_back(selfAlignInfo);
			
			int mateEnd = OtherReadEnd(alignment.readEnd);
			
			int seedScore;
			int seedPosition;
			AlignMate3PrimeSeed16SSE(aligner, alignment, referenceSequences, preppedReads, minFragmentLength, maxFragmentLength, seedScore, seedPosition);
			
			if (seedScore >= cSeedScoreThreshold)
			{
				AlignInfo mateFwdAlignInfo = AlignFwdMateFullSSE(aligner, alignment, referenceSequences, preppedReads, seedPosition);
				
				mateFwdAlignments[alignmentIndex] = mateFwdAlignInfo;
				
				AlignInfo mateRevAlignInfo = AlignRevMateFullSSE(aligner, alignment, referenceSequences, preppedReads, mateFwdAlignInfo.AlignmentPosition(preppedReads.ReadLength(mateEnd)));
				
				mateRevAlignments[alignmentIndex] = mateRevAlignInfo;
			}
		}
		
		// Initialize posterior calculation for all alignments based on length of best alignment
		int alignedLength[2] = {0, 0};
		AlignmentPosterior alignPosteriors[2];
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			alignedLength[readEnd] = bestAlignment[readEnd].second - cBreakEndAdjust;
			alignPosteriors[readEnd].Initialize(&alignProbability, alignedLength[readEnd]);
		}
		
		// Identify best alignment score
		// Add alignment scores to posterior calculation
		int bestSelfScore[2] = {0, 0};
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			int selfScore = selfAlignments[alignmentIndex].SeqScores()[alignedLength[alignment.readEnd]];
			
			alignPosteriors[alignment.readEnd].AddAlignment(selfScore);
			bestSelfScore[alignment.readEnd] = max(bestSelfScore[alignment.readEnd], selfScore);
			
			unordered_map<int,AlignInfo>::const_iterator mateAlignIter = mateRevAlignments.find(alignmentIndex);
			if (mateAlignIter != mateRevAlignments.end())
			{
				int mateEnd = OtherReadEnd(alignment.readEnd);
				
				alignPosteriors[mateEnd].AddAlignment(mateAlignIter->second.SeqScores()[alignedLength[mateEnd]]);
			}
		}
		
		// Calculate probability the read is valid
		double validReadPosterior[2] = {0.0, 0.0};
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			validReadPosterior[readEnd] = alignProbability.Classify(alignedLength[readEnd], bestSelfScore[readEnd], validReadPrior);
		}
		
		// Apply threshold on valid read posterior
		if (validReadPosterior[0] * validReadPosterior[1] < validReadThreshold)
		{
			continue;
		}
		
		// Calculate probability for best discordant score
		double probBestDiscordant = alignPosteriors[0].Posterior(bestSelfScore[0]) * alignPosteriors[1].Posterior(bestSelfScore[1]);
		
		// Calculate probability for best concordant score
		double probBestConcordant = 0.0;
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			unordered_map<int,AlignInfo>::const_iterator mateAlignIter = mateRevAlignments.find(alignmentIndex);
			if (mateAlignIter != mateRevAlignments.end())
			{
				int mateEnd = OtherReadEnd(alignment.readEnd);
				
				int selfScore = selfAlignments[alignmentIndex].SeqScores()[alignedLength[alignment.readEnd]];
				int mateScore = mateAlignIter->second.SeqScores()[alignedLength[mateEnd]];
				
				probBestConcordant = max(probBestConcordant, alignPosteriors[alignment.readEnd].Posterior(selfScore) * alignPosteriors[mateEnd].Posterior(mateScore));
			}
		}
		
		// Calculate probability of concordance
		double chimericPosterior = probBestDiscordant * chimericPrior / (probBestDiscordant * chimericPrior + probBestConcordant * (1.0 - chimericPrior));
		
		// Apply threshold on chimeric posterior
		if (chimericPosterior < chimericThreshold)
		{
			continue;
		}
		
		// Create list of alignment indices for each end
		vector<int> alignmentIndices[2];
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			AlignInfo alignInfo = selfAlignments[alignmentIndex];
			
			int score = alignInfo.SeqScores()[alignedLength[alignment.readEnd]];
			
			double alignmentPosterior = alignPosteriors[alignment.readEnd].Posterior(score);
			
			if (alignmentPosterior < alignmentThreshold)
			{
				continue;
			}
			
			alignmentIndices[alignment.readEnd].push_back(alignmentIndex);
		}
		
		if (alignmentIndices[0].empty() || alignmentIndices[1].empty())
		{
			continue;
		}
		
		// Output spanning alignments exceeding posterior threshold
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			AlignInfo alignInfo = selfAlignments[alignmentIndex];
			
			int seqLength = alignedLength[alignment.readEnd];
			int score = alignInfo.SeqScores()[alignedLength[alignment.readEnd]];
			
			double alignmentPosterior = alignPosteriors[alignment.readEnd].Posterior(score);
			
			if (alignmentPosterior < alignmentThreshold)
			{
				continue;
			}
			
			spanningFile << alignment.fragment << "\t";
			spanningFile << alignment.readEnd << "\t";
			spanningFile << alignment.reference << "\t";
			spanningFile << ((alignment.strand == PlusStrand) ? "+" : "-") << "\t";
			spanningFile << alignInfo.AlignmentStart(seqLength) << "\t";
			spanningFile << alignInfo.AlignmentEnd(seqLength) << "\t";
			spanningFile << preppedReads.ReadLength(alignment.readEnd) << "\t";
			spanningFile << seqLength << "\t";
			spanningFile << score << "\t";
			spanningFile << alignmentPosterior << "\t";
			spanningFile << chimericPosterior << "\t";
			spanningFile << validReadPosterior[alignment.readEnd] << endl;
		}
		
		// Output split alignments exceeding posterior threshold
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			int mateEnd = OtherReadEnd(readEnd);
			
			for (vector<int>::const_iterator selfAlignmentIter = alignmentIndices[readEnd].begin(); selfAlignmentIter != alignmentIndices[readEnd].end(); selfAlignmentIter++)
			{
				const AlignInfo& selfAlignInfo = selfAlignments[*selfAlignmentIter];
				const RawAlignment& selfAlignment = alignments[*selfAlignmentIter];
				
				for (vector<int>::const_iterator mateAlignmentIter = alignmentIndices[mateEnd].begin(); mateAlignmentIter != alignmentIndices[mateEnd].end(); mateAlignmentIter++)
				{
					unordered_map<int,AlignInfo>::const_iterator mateAlignInfoIter = mateFwdAlignments.find(*mateAlignmentIter);
					
					if (mateAlignInfoIter == mateFwdAlignments.end())
					{
						continue;
					}
					
					const AlignInfo& mateAlignInfo = mateAlignInfoIter->second;
					const RawAlignment& mateAlignment = alignments[*mateAlignmentIter];
					
					IntegerVec seq1Length;
					IntegerVec seq2Length;
					int score = BestSplitAlignment(selfAlignInfo.SeqScores(), selfAlignInfo.SeqScoresLength(), mateAlignInfo.SeqScores(), mateAlignInfo.SeqScoresLength(), -1, seq1Length, seq2Length);
					
					if (seq1Length.size() > cMaxBreakpointHomology)
					{
						continue;
					}
					
					string readSeq = preppedReads.Sequence(readEnd);
					
					for (int i = 0; i < seq1Length.size(); i++)
					{
						if (seq1Length[i] < cMinAnchor || seq2Length[i] < cMinAnchor)
						{
							continue;
						}
						
						splitFile << selfAlignment.fragment << "\t";
						splitFile << selfAlignment.readEnd << "\t";
						splitFile << selfAlignment.reference << "\t";
						splitFile << ((selfAlignment.strand == PlusStrand) ? "+" : "-") << "\t";
						splitFile << selfAlignInfo.BreakPosition(seq1Length[i]) << "\t";
						splitFile << mateAlignment.reference << "\t";
						splitFile << ((mateAlignment.strand == PlusStrand) ? "+" : "-") << "\t";
						splitFile << mateAlignInfo.BreakPosition(seq2Length[i]) << "\t";
						splitFile << readSeq.substr(seq1Length[i], readSeq.size() - seq2Length[i] - seq1Length[i]) << "\t";
						splitFile << seq1Length[i] << "\t";
						splitFile << seq2Length[i] << "\t";
						splitFile << selfAlignInfo.SeqScores()[seq1Length[i]] << "\t";
						splitFile << mateAlignInfo.SeqScores()[seq2Length[i]] << "\t";
						splitFile << score << endl;
					}
				}
			}
		}
	}
}





