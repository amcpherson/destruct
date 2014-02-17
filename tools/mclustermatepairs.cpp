/*
 *  mclustermatepairs.cpp
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
#include "Parsers.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


void ReadInstructions(const string& instructions, vector<string>& libNames, vector<double>& fragmentLengthMeans, vector<double>& fragmentLengthStdDevs, vector<vector<string> >& alignmentFilenames, vector<vector<string> >& alignmentFormats)
{
	istream* instructionsFile = 0;
	if (instructions == "-")
	{
		instructionsFile = &cin;
	}
	else
	{
		instructionsFile = new ifstream(instructions.c_str());
	}
	
	unordered_map<string,double> libFragmentLengthMeans;
	unordered_map<string,double> libFragmentLengthStdDevs;
	unordered_map<string,vector<string> > libAlignmentFilenames;
	unordered_map<string,vector<string> > libAlignmentFormats;
	
	StringVec fields;
	while (ReadTSV(*instructionsFile, fields))
	{
		if (fields.size() < 4)
		{
			cerr << "Must be 4 fields per line" << endl;
			cerr << "For an library entry fields are: 'library'{tab}library_name{tab}fragment_mean{tab}fragment_stddev" << endl;
			cerr << "For an alignment entry fields are: 'alignments'{tab}library_name{tab}alignment_filename{tab}alignment_format" << endl;
			exit(1);
		}
		
		const string& instype = fields[0];
		const string& library = fields[1];
		
		if (instype == "library")
		{
			libFragmentLengthMeans[library] = SAFEPARSE(double, fields[2]);
			libFragmentLengthStdDevs[library] = SAFEPARSE(double, fields[3]);
		}
		else if (instype == "alignments")
		{
			if (fields[3] != "sam" && fields[3] != "compact")
			{
				cerr << "Alignment format " << fields[3] << " is not supported" << endl;
				cerr << "Supported alignment formats are 'sam' and 'compact'" << endl;
				exit(1);
			}
			
			libAlignmentFilenames[library].push_back(fields[2]);
			libAlignmentFormats[library].push_back(fields[3]);
		}
		else
		{
			cerr << "First field must be either 'library' or 'alignments'" << endl;
			exit(1);
		}
	}
	
	if (instructions != "-")
	{
		delete instructionsFile;
	}
	
	for (unordered_map<string,double>::const_iterator libIter = libFragmentLengthMeans.begin(); libIter != libFragmentLengthMeans.end(); libIter++)
	{
		libNames.push_back(libIter->first);
		fragmentLengthMeans.push_back(libFragmentLengthMeans[libIter->first]);
		fragmentLengthStdDevs.push_back(libFragmentLengthStdDevs[libIter->first]);
		alignmentFilenames.push_back(libAlignmentFilenames[libIter->first]);
		alignmentFormats.push_back(libAlignmentFormats[libIter->first]);
	}
}

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

void CreateCompactAlignments(const RawAlignmentVec& alignments, CompAlignVec& compactAlignments, NameIndex& refNameIndex, int fragmentIndex)
{
	for (RawAlignmentVecConstIter alignmentIter = alignments.begin(); alignmentIter != alignments.end(); alignmentIter++)
	{
		compactAlignments.push_back(CompactAlignment());
		compactAlignments.back().readID.fragmentIndex = fragmentIndex;
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

void CreateMatePairs(const CompAlignVec& alignments1, const CompAlignVec alignments2, const vector<unsigned char>& fragmentLibIndices, const vector<double>& fragmentLengthMeans, const vector<double>& fragmentLengthStdDevs, MatePairVec& matePairs)
{
	for (int alignPairIndex = 0; alignPairIndex < alignments1.size(); alignPairIndex++)
	{
		const CompactAlignment& alignment1 = alignments1[alignPairIndex];
		const CompactAlignment& alignment2 = alignments2[alignPairIndex];
		
		DebugCheck(alignment1.readID.fragmentIndex == alignment2.readID.fragmentIndex);
		DebugCheck(alignment1.readID.readEnd != alignment2.readID.readEnd);
		
		int libIndex = fragmentLibIndices[alignment1.readID.fragmentIndex];
		
		int length1 = alignment1.region.end - alignment1.region.start + 1;
		int length2 = alignment2.region.end - alignment2.region.start + 1;
		
		MatePair matePair;
		matePair.x = (alignment1.refStrand.strand == PlusStrand) ? alignment1.region.end : -alignment1.region.start;
		matePair.y = (alignment2.refStrand.strand == PlusStrand) ? alignment2.region.end : -alignment2.region.start;
		matePair.u = fragmentLengthMeans[libIndex] - length1 - length2;
		matePair.s = fragmentLengthStdDevs[libIndex];
		
		matePairs.push_back(matePair);
	}
}

void OutputClusterMember(ostream& out, int clusterID, int clusterEnd, const CompactAlignment& alignment, const vector<string>& referenceNames, const vector<int>& fragmentIndices, const vector<unsigned char>& fragmentLibIndices, const vector<string>& libraryNames)
{
	out << clusterID << "\t";
	out << clusterEnd << "\t";
	out << fragmentIndices[alignment.readID.fragmentIndex] << "\t";
	out << alignment.readID.readEnd << "\t";
	out << referenceNames[alignment.refStrand.referenceIndex] << "\t";
	out << ((alignment.refStrand.strand == PlusStrand) ? "+" : "-") << "\t";
	out << alignment.region.start << "\t";
	out << alignment.region.end << "\t";
	out << alignment.alignProb << "\t";
	out << alignment.chimericProb << "\t";
	out << alignment.validProb << "\t";
	out << libraryNames[fragmentLibIndices[alignment.readID.fragmentIndex]] << endl;
}

int main(int argc, char* argv[])
{
	string instructions;
	string clustersFilename;
	int minClusterSize;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Pair Clustering Tool");
		TCLAP::ValueArg<string> instructionsArg("i","instructions","Instructions for Clustering Job",true,"","string",cmd);
		TCLAP::ValueArg<string> clustersFilenameArg("o","output","Output Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> minClusterSizeArg("m","minclustersize","Minimum Cluster Size",true,-1,"integer",cmd);
		cmd.parse(argc,argv);
		
		instructions = instructionsArg.getValue();
		clustersFilename = clustersFilenameArg.getValue();
		minClusterSize = minClusterSizeArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	vector<string> libraryNames;
	vector<double> fragmentLengthMeans;
	vector<double> fragmentLengthStdDevs;
	vector<vector<string> > alignmentFilenames;
	vector<vector<string> > alignmentFormats;
	ReadInstructions(instructions, libraryNames, fragmentLengthMeans, fragmentLengthStdDevs, alignmentFilenames, alignmentFormats);
	
	DiscordantAlignments discordantAlignments;
	
	NameIndex refNameIndex;
	
	vector<int> remapFragmentIndices;
	vector<unsigned char> fragmentLibIndices;
	
	for (unsigned char libIndex = 0; libIndex < (unsigned char)libraryNames.size(); libIndex++)
	{
		double fragmentLengthMean = fragmentLengthMeans[libIndex];
		double fragmentLengthStdDev = fragmentLengthStdDevs[libIndex];
		
		for (int filenameIndex = 0; filenameIndex < alignmentFilenames[libIndex].size(); filenameIndex++)
		{
			const string& alignmentFilename = alignmentFilenames[libIndex][filenameIndex];
			const string& alignmentFormat = alignmentFormats[libIndex][filenameIndex];
			
			cerr << "Reading alignments from " << alignmentFilename << endl;
			
			AlignmentStream* alignmentStream = 0;
			if (alignmentFormat == "sam")
			{
				alignmentStream = new SamAlignmentStream(alignmentFilename);
			}
			else if (alignmentFormat == "compact")
			{
				alignmentStream = new CompactAlignmentStream(alignmentFilename);
			}
			
			FragmentAlignmentStream fragmentAlignmentStream(alignmentStream);
			
			const int maxFragmentLength = (int)(fragmentLengthMean + 4 * fragmentLengthStdDev);
			
			RawAlignmentVec alignments;
			while (fragmentAlignmentStream.GetNextAlignments(alignments))
			{
				if (IsConcordant(alignments, maxFragmentLength))
				{
					continue;
				}
				
				int fragmentIndex = remapFragmentIndices.size();
				remapFragmentIndices.push_back(SAFEPARSE(int, alignments.front().fragment));
				fragmentLibIndices.push_back(libIndex);
				
				CompAlignVec compactAlignments;
				CreateCompactAlignments(alignments, compactAlignments, refNameIndex, fragmentIndex);
				
				discordantAlignments.AddFragmentAlignments(compactAlignments, fragmentLengthMean, fragmentLengthStdDev);
			}
			
			delete alignmentStream;
		}
	}

	cerr << "Read alignments for " << remapFragmentIndices.size() << " fragments" << endl;
	
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
			OutputClusterMember(clustersFile, clusterID, 0, alignments1.front(), refNameIndex.Get(), remapFragmentIndices, fragmentLibIndices, libraryNames);
			OutputClusterMember(clustersFile, clusterID, 1, alignments2.front(), refNameIndex.Get(), remapFragmentIndices, fragmentLibIndices, libraryNames);
			clusterID++;
			continue;
		}
		
		MatePairVec matePairs;
		CreateMatePairs(alignments1, alignments2, fragmentLibIndices, fragmentLengthMeans, fragmentLengthStdDevs, matePairs);
		
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
			
			for (int clusterIndex = 0; clusterIndex < gibbsClusters.size(); clusterIndex++)
			{
				const IntegerVec& cluster = gibbsClusters[clusterIndex];
				
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
					
					OutputClusterMember(clustersFile, clusterID, 0, alignment1, refNameIndex.Get(), remapFragmentIndices, fragmentLibIndices, libraryNames);
					OutputClusterMember(clustersFile, clusterID, 1, alignment2, refNameIndex.Get(), remapFragmentIndices, fragmentLibIndices, libraryNames);
				}
				
				clusterID++;
			}
		}
	}
	
	clustersFile.close();
	
	cout << "Created " << clusterID << " clusters" << endl;
}

