/*
 *  clustermatepairs.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "AlignmentStream.h"
#include "DebugCheck.h"
#include "Indexer.h"
#include "DiscordantAlignments.h"
#include "MatePairDelauny.h"
#include "MatePairGibbs.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


bool IsConcordant(const RawAlignmentVec& alignments, int maxFragmentLength)
{
	for (RawAlignmentVecConstIter alignmentIter1 = alignments.begin(); alignmentIter1 != alignments.end(); alignmentIter1++)
	{
		for (RawAlignmentVecConstIter alignmentIter2 = alignmentIter1 + 1; alignmentIter2 != alignments.end(); alignmentIter2++)
		{
			if (alignmentIter1->readEnd == alignmentIter2->readEnd)
			{
				continue;
			}
			
			if (alignmentIter1->reference != alignmentIter2->reference)
			{
				continue;
			}
			
			if (alignmentIter1->strand == PlusStrand && alignmentIter2->strand == MinusStrand)
			{
				int inferredLength = alignmentIter2->region.end - alignmentIter1->region.start + 1;
				
				if (inferredLength >= 0 && inferredLength < maxFragmentLength)
				{
					return true;
				}
			}
			
			if (alignmentIter2->strand == PlusStrand && alignmentIter1->strand == MinusStrand)
			{
				int inferredLength = alignmentIter1->region.end - alignmentIter2->region.start + 1;
				
				if (inferredLength >= 0 && inferredLength < maxFragmentLength)
				{
					return true;
				}
			}
		}
	}
	
	return false;
}

void CreateCompactAlignments(const RawAlignmentVec& alignments, CompAlignVec& compactAlignments, NameIndex& refNameIndex)
{
	for (RawAlignmentVecConstIter alignmentIter = alignments.begin(); alignmentIter != alignments.end(); alignmentIter++)
	{
		compactAlignments.push_back(CompactAlignment());
		compactAlignments.back().readID.fragmentIndex = SAFEPARSE(int, alignmentIter->fragment);
		compactAlignments.back().readID.readEnd = alignmentIter->readEnd;
		compactAlignments.back().refStrand.referenceIndex = refNameIndex.Index(alignmentIter->reference);
		compactAlignments.back().refStrand.strand = alignmentIter->strand;
		compactAlignments.back().region.start = alignmentIter->region.start;
		compactAlignments.back().region.end = alignmentIter->region.end;
		compactAlignments.back().alignProb = alignmentIter->alignProb;		
		compactAlignments.back().chimericProb = alignmentIter->chimericProb;		
		compactAlignments.back().validProb = alignmentIter->validProb;		
	}
}

void CreateMatePairs(const CompAlignVec& alignments1, const CompAlignVec alignments2, double fragmentLengthMean, double fragmentLengthStdDev, MatePairVec& matePairs)
{
	for (int alignPairIndex = 0; alignPairIndex < alignments1.size(); alignPairIndex++)
	{
		const CompactAlignment& alignment1 = alignments1[alignPairIndex];
		const CompactAlignment& alignment2 = alignments2[alignPairIndex];
		
		DebugCheck(alignment1.readID.fragmentIndex == alignment2.readID.fragmentIndex);
		DebugCheck(alignment1.readID.readEnd != alignment2.readID.readEnd);
		
		int length1 = alignment1.region.end - alignment1.region.start + 1;
		int length2 = alignment2.region.end - alignment2.region.start + 1;
		
		MatePair matePair;
		matePair.x = (alignment1.refStrand.strand == PlusStrand) ? alignment1.region.end : -alignment1.region.start;
		matePair.y = (alignment2.refStrand.strand == PlusStrand) ? alignment2.region.end : -alignment2.region.start;
		matePair.u = fragmentLengthMean - length1 - length2;
		matePair.s = fragmentLengthStdDev;
		
		matePairs.push_back(matePair);
	}
}

void OutputClusterMember(ostream& out, int clusterID, int clusterEnd, const CompactAlignment& alignment, const StringVec& referenceNames)
{
	out << clusterID << "\t";
	out << clusterEnd << "\t";
	out << alignment.readID.fragmentIndex << "\t";
	out << alignment.readID.readEnd << "\t";
	out << referenceNames[alignment.refStrand.referenceIndex] << "\t";
	out << ((alignment.refStrand.strand == PlusStrand) ? "+" : "-") << "\t";
	out << alignment.region.start << "\t";
	out << alignment.region.end << "\t";
	out << alignment.alignProb << "\t";
	out << alignment.chimericProb << "\t";
	out << alignment.validProb << endl;
}

void OutputClusters(ostream& out, int& clusterID, const NameIndex& refNameIndex, const CompAlignVec& alignments1, const CompAlignVec alignments2, const IntegerTable& clusters, const IntegerVec& alignPairIndices)
{
	for (int clusterIndex = 0; clusterIndex < clusters.size(); clusterIndex++)
	{
		const IntegerVec& cluster = clusters[clusterIndex];
		
		unordered_set<int> clusterFragmentIndices;
		for (int elementIndex = 0; elementIndex < cluster.size(); elementIndex++)
		{
			int alignPairIndex = alignPairIndices[cluster[elementIndex]];
			
			const CompactAlignment& alignment1 = alignments1[alignPairIndex];
			const CompactAlignment& alignment2 = alignments2[alignPairIndex];
			
			DebugCheck(alignment1.readID.fragmentIndex == alignment2.readID.fragmentIndex);
			DebugCheck(alignment1.readID.readEnd != alignment2.readID.readEnd);
			
			// For fragments with multiple alignments supporting the same cluster, select first alignment
			if (!clusterFragmentIndices.insert(alignment1.readID.fragmentIndex).second)
			{
				continue;
			}
			
			OutputClusterMember(out, clusterID, 0, alignment1, refNameIndex.Get());
			OutputClusterMember(out, clusterID, 1, alignment2, refNameIndex.Get());
		}
		
		clusterID++;
	}
}

int main(int argc, char* argv[])
{
	string samFilename;
	string alignFilename;
	string clustersFilename;
	double fragmentLengthMean;
	double fragmentLengthStdDev;
	int minClusterSize;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Pair Clustering Tool");
		TCLAP::ValueArg<string> samFilenameArg("s","sam","Read Sorted Sam Filename",false,"","string");
		TCLAP::ValueArg<string> alignFilenameArg("a","align","Read Sorted Compact Alignments Filename",false,"","string");
		TCLAP::ValueArg<string> clustersFilenameArg("o","output","Output Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> fragmentLengthMeanArg("u","fragmentmean","Fragment Length Mean",true,-1,"integer",cmd);
		TCLAP::ValueArg<double> fragmentLengthStdDevArg("d","fragmentstddev","Fragment Length Standard Deviation",true,-1,"integer",cmd);
		TCLAP::ValueArg<int> minClusterSizeArg("m","minclustersize","Minimum Cluster Size",true,-1,"integer",cmd);
		
		vector<TCLAP::Arg*> alignmentArgs;
		alignmentArgs.push_back(&samFilenameArg);
		alignmentArgs.push_back(&alignFilenameArg);
		cmd.xorAdd(alignmentArgs);
		
		cmd.parse(argc,argv);
		
		samFilename = samFilenameArg.getValue();
		alignFilename = alignFilenameArg.getValue();
		clustersFilename = clustersFilenameArg.getValue();
		fragmentLengthMean = fragmentLengthMeanArg.getValue();
		fragmentLengthStdDev = fragmentLengthStdDevArg.getValue();
		minClusterSize = minClusterSizeArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	AlignmentStream* alignmentStream = 0;
	if (!samFilename.empty())
	{
		alignmentStream = new SamAlignmentStream(samFilename);
	}
	else if (!alignFilename.empty())
	{
		alignmentStream = new CompactAlignmentStream(alignFilename);
	}
	
	FragmentAlignmentStream fragmentAlignmentStream(alignmentStream);
	
	cout << "Finding pairs of reference sequences connected by pairs of alignments" << endl;
	
	const int maxFragmentLength = (int)(fragmentLengthMean + 4 * fragmentLengthStdDev);
	
	DiscordantAlignments discordantAlignments;
	
	NameIndex refNameIndex;
	
	RawAlignmentVec alignments;
	while (fragmentAlignmentStream.GetNextAlignments(alignments))
	{
		if (IsConcordant(alignments, maxFragmentLength))
		{
			continue;
		}
		
		CompAlignVec compactAlignments;
		CreateCompactAlignments(alignments, compactAlignments, refNameIndex);
		
		discordantAlignments.AddFragmentAlignments(compactAlignments, fragmentLengthMean, fragmentLengthStdDev);
	}
	
	// Open output clusters file
	ofstream clustersFile(clustersFilename.c_str());
	if (!clustersFile)
	{
		cerr << "Error: unable to write to clusters file" << endl;		
		exit(1);
	}
	
	cout << "Creating clusters" << endl;
	
	int clusterID = 0;
	for (discordantAlignments.StartRefStrandIteration(); !discordantAlignments.FinishedRefStrandIteration(); discordantAlignments.NextRefStrandIteration())
	{
		CompAlignVec alignments1;
		CompAlignVec alignments2;
		discordantAlignments.RetrieveRefStrandAlignments(alignments1, alignments2);
		
		if (alignments1.size() < minClusterSize)
		{
			continue;
		}
		
		if (alignments1.size() == 1)
		{
			OutputClusterMember(clustersFile, clusterID, 0, alignments1.front(), refNameIndex.Get());
			OutputClusterMember(clustersFile, clusterID, 1, alignments2.front(), refNameIndex.Get());
			clusterID++;
			continue;
		}
		
		MatePairVec matePairs;
		CreateMatePairs(alignments1, alignments2, fragmentLengthMean, fragmentLengthStdDev, matePairs);
		
		MatePairDelauny delaunyClusterer;
		
		IntegerTable delaunyClusters;
		delaunyClusterer.DoClustering(matePairs, delaunyClusters);
		
		for (IntegerTableConstIter delaunyClusterIter = delaunyClusters.begin(); delaunyClusterIter != delaunyClusters.end(); delaunyClusterIter++)
		{
			MatePairVec delaunyMatePairs;
			IntegerVec alignPairIndices;
			for (IntegerVecConstIter elementIter = delaunyClusterIter->begin(); elementIter != delaunyClusterIter->end(); elementIter++)
			{
				delaunyMatePairs.push_back(matePairs[*elementIter]);
				alignPairIndices.push_back(*elementIter);
			}
			
			MatePairGibbs gibbsClusterer;
			
			IntegerTable gibbsClusters;
			gibbsClusterer.DoClustering(delaunyMatePairs, gibbsClusters);
			
			OutputClusters(clustersFile, clusterID, refNameIndex, alignments1, alignments2, gibbsClusters, alignPairIndices);
		}
	}
	
	clustersFile.close();
	
	cout << "Created " << clusterID << " clusters" << endl;
}

