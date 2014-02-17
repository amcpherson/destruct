/*
 *  realignp.cpp
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
#include "ExonRegions.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


struct CompactGenomicLoci
{
	union
	{
		struct
		{
			RefStrand refStrand;
			int position;
		};
		
		uint64_t id;
	};
};

uint64_t CalculateGenomicLociID(const ExonRegions& exonRegions, NameIndex& refNameIndex, const RawAlignment& alignment, const AlignInfo& alignInfo)
{
	CompactGenomicLoci genomicLoci;
	
	string gene;
	string transcript;        
	if (ParseTranscriptID(alignment.reference, gene, transcript) && exonRegions.IsTranscript(transcript))
 	{
 		string chromosome;
 		int strand;
 		int position;
		exonRegions.RemapTranscriptToGenome(transcript, alignment.strand, alignInfo.AlignmentPosition(0), chromosome, strand, position);
		
		genomicLoci.refStrand.referenceIndex = refNameIndex.Index(chromosome);
		genomicLoci.refStrand.strand = strand;
		genomicLoci.position = position;
	}
	else
	{
		genomicLoci.refStrand.referenceIndex = refNameIndex.Index(alignment.reference);
		genomicLoci.refStrand.strand = alignment.strand;
		genomicLoci.position = alignInfo.AlignmentPosition(0);
	}
	
	return genomicLoci.id;
}

int main(int argc, char* argv[])
{
	int matchScore;
	int misMatchScore;
	int gapScore;
	string referenceFasta;
	string readSeqsFilename;
	string alignmentsFilename;
	string statsFilename;
	string exonRegionsFilename;
	double chimericPrior;
	double alignmentThreshold;
	double chimericThreshold;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Realignment Tool");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqsFilenameArg("s","seqs","Read Sequences Fastq",true,"","string",cmd);
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Sam Alignments",true,"","string",cmd);
		TCLAP::ValueArg<string> statsFilenameArg("z","stats","Stats Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> exonRegionsFilenameArg("e","exons","Exon Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> chimericPriorArg("","pchimer","Prior Probility of Chimeric Read",true,0.05,"float",cmd);
		TCLAP::ValueArg<double> alignmentThresholdArg("","talign","Alignment Posterior Threshold",true,0.1,"float",cmd);
		TCLAP::ValueArg<double> chimericThresholdArg("","tchimer","Chimeric Posterior Threshold",true,0.1,"float",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		readSeqsFilename = readSeqsFilenameArg.getValue();
		alignmentsFilename = alignmentsFilenameArg.getValue();
		statsFilename = statsFilenameArg.getValue();
		exonRegionsFilename = exonRegionsFilenameArg.getValue();
		chimericPrior = chimericPriorArg.getValue();
		alignmentThreshold = alignmentThresholdArg.getValue();
		chimericThreshold = chimericThresholdArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	const int cSeedScoreThreshold = 8;
	const int cMateSearchLength = 1000;
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
	
	cerr << "Reading exon regions" << endl;
	
	ExonRegions exonRegions;
	exonRegions.Read(exonRegionsFilename);
	
	NameIndex refNameIndex;
	
	cerr << "Realigning" << endl;
	
	SimpleAligner aligner(matchScore, misMatchScore, gapScore);
	
	SamAlignmentStream alignmentStream(alignmentsFilename);
	FragmentAlignmentStream fragmentAlignmentStream(&alignmentStream);
	
	RawAlignmentVec alignments;
	while (fragmentAlignmentStream.GetNextAlignments(alignments))
	{
		preppedReads.SetCurrentRead(SAFEPARSE(int, alignments.front().fragment));
		
		// Realignments
		pair<int,int> bestAlignment[2] = {pair<int,int>(0,0),pair<int,int>(0,0)};
		vector<AlignInfo> selfAlignments;
		unordered_map<int,AlignInfo> mateAlignments;
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
			AlignMate3PrimeSeed16SSE(aligner, alignment, referenceSequences, preppedReads, cMateSearchLength, seedScore, seedPosition);
			
			if (seedScore >= cSeedScoreThreshold)
			{
				AlignInfo mateFwdAlignInfo = AlignFwdMateFullSSE(aligner, alignment, referenceSequences, preppedReads, seedPosition);
				
				AlignInfo mateRevAlignInfo = AlignRevMateFullSSE(aligner, alignment, referenceSequences, preppedReads, mateFwdAlignInfo.AlignmentPosition(preppedReads.ReadLength(mateEnd)));
				
				mateAlignments[alignmentIndex] = mateRevAlignInfo;
			}
		}
		
		// Initialize posterior calculation for all alignments based on length of best alignment
		int alignedLength[2] = {0, 0};
		AlignmentPosterior alignPosteriors[2];
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			if (bestAlignment[readEnd].second == 0)
			{
				continue;
			}
			
			alignedLength[readEnd] = bestAlignment[readEnd].second - cBreakEndAdjust;
			alignPosteriors[readEnd].Initialize(&alignProbability, alignedLength[readEnd]);
		}
		
		// Identify best alignment scores at each remapped genomic position
		int bestSelfScore[2] = {0, 0};
		unordered_map<uint64_t,int> bestGenomicScore[2];
		for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
		{
			const RawAlignment& alignment = alignments[alignmentIndex];
			
			uint64_t selfGenomicLociID = CalculateGenomicLociID(exonRegions, refNameIndex, alignment, selfAlignments[alignmentIndex]);
			int selfScore = selfAlignments[alignmentIndex].SeqScores()[alignedLength[alignment.readEnd]];
			bestSelfScore[alignment.readEnd] = max(bestSelfScore[alignment.readEnd], selfScore);
			
			unordered_map<uint64_t,int>::iterator selfScoreIter = bestGenomicScore[alignment.readEnd].insert(make_pair(selfGenomicLociID, selfScore)).first;
			selfScoreIter->second = max(selfScoreIter->second, selfScore);
			
			int mateEnd = OtherReadEnd(alignment.readEnd);
			
			if (alignedLength[mateEnd] == 0)
			{
				continue;
			}
			
			unordered_map<int,AlignInfo>::const_iterator mateAlignIter = mateAlignments.find(alignmentIndex);
			if (mateAlignIter != mateAlignments.end())
			{
				const AlignInfo& mateAlignInfo = mateAlignIter->second;
				
				uint64_t mateGenomicLociID = CalculateGenomicLociID(exonRegions, refNameIndex, alignment, mateAlignInfo);
				int mateScore = mateAlignInfo.SeqScores()[alignedLength[mateEnd]];
				
				unordered_map<uint64_t,int>::iterator mateScoreIter = bestGenomicScore[mateEnd].insert(make_pair(mateGenomicLociID, mateScore)).first;
				mateScoreIter->second = max(mateScoreIter->second, mateScore);
			}
		}
		
		// Add best alignment scores for each genomic loci to posterior calculation
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			for (unordered_map<uint64_t,int>::const_iterator scoreIter = bestGenomicScore[readEnd].begin(); scoreIter != bestGenomicScore[readEnd].end(); scoreIter++)
			{
				alignPosteriors[readEnd].AddAlignment(scoreIter->second);
			}
		}
		
		// Check for concordant alignments
		if (alignedLength[0] != 0 && alignedLength[1] != 0)
		{
			// Calculate probability for best discordant score
			double probBestDiscordant = alignPosteriors[0].Posterior(bestSelfScore[0]) * alignPosteriors[1].Posterior(bestSelfScore[1]);
			
			// Calculate probability for best concordant score
			double probBestConcordant = 0.0;
			for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
			{
				const RawAlignment& alignment = alignments[alignmentIndex];
				
				unordered_map<int,AlignInfo>::const_iterator mateAlignIter = mateAlignments.find(alignmentIndex);
				if (mateAlignIter != mateAlignments.end())
				{
					int mateEnd = OtherReadEnd(alignment.readEnd);
					
					int selfScore = selfAlignments[alignmentIndex].SeqScores()[alignedLength[alignment.readEnd]];
					int mateScore = mateAlignIter->second.SeqScores()[alignedLength[mateEnd]];
					
					probBestConcordant = max(probBestConcordant, alignPosteriors[alignment.readEnd].Posterior(selfScore) * alignPosteriors[mateEnd].Posterior(mateScore));
				}
			}
			
			// Calculate probability of concordance
			double chimericPosterior = probBestDiscordant * chimericPrior / (probBestDiscordant * chimericPrior + probBestConcordant * (1.0 - chimericPrior));
			
			if (chimericPosterior < chimericThreshold)
			{
				continue;
			}
		}
		
		// Output alignments exceeding posterior threshold
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
			
			double probability = alignmentPosterior;
			
			cout << alignment.fragment << "\t";
			cout << alignment.readEnd << "\t";
			cout << alignment.reference << "\t";
			cout << ((alignment.strand == PlusStrand) ? "+" : "-") << "\t";
			cout << alignInfo.AlignmentStart(seqLength) << "\t";
			cout << alignInfo.AlignmentEnd(seqLength) << "\t";
			cout << preppedReads.ReadLength(alignment.readEnd) << "\t";
			cout << seqLength << "\t";
			cout << score << "\t";
			cout << probability << endl;
		}
	}
}




