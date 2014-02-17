/*
 *  mpredictbreaks.cpp
 *
 */

#include "AlignmentStream.h"
#include "DebugCheck.h"
#include "Indexer.h"
#include "DiscordantAlignments.h"
#include "MatePairDelauny.h"
#include "MatePairGibbs.h"
#include "Parsers.h"
#include "Sequences.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


struct CompactPosition
{
	RefStrand refStrand;
	int position;
};

inline bool operator==(const CompactPosition& position1, const CompactPosition& position2)
{
	return position1.refStrand.id == position2.refStrand.id && position1.position == position2.position;
}

inline size_t hash_value(const CompactPosition& position)
{
	size_t seed = 0;
	hash_combine(seed, position.refStrand.id);
	hash_combine(seed, position.position);
    return seed;
}

struct CompactSplitAlignment
{
	CompactPosition position[2];
	string inserted;
	int score;
};

struct CompactSplitDefinition
{
	CompactPosition position[2];
	int numInserted;
};

inline bool operator==(const CompactSplitDefinition& split1, const CompactSplitDefinition& split2)
{
	return split1.position[0] == split2.position[0] && split1.position[1] == split2.position[1] && split1.numInserted == split2.numInserted;
}

inline size_t hash_value(const CompactSplitDefinition& split)
{
	size_t seed = 0;
	hash_combine(seed, split.position[0]);
	hash_combine(seed, split.position[1]);
	hash_combine(seed, split.numInserted);
    return seed;
}

const int cMateSearchLength = 1000;

bool Match(const LocationVec& clusterLocations, NameIndex& refNameIndex, const CompactSplitAlignment& splitAlignment, int flip)
{
	for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
	{
		int splitEnd = (clusterEnd + flip) % 2;
		
		int referenceIndex = refNameIndex.Index(clusterLocations[clusterEnd].refName);
		
		if (splitAlignment.position[splitEnd].refStrand.referenceIndex != referenceIndex)
		{
			return false;
		}
		
		if (splitAlignment.position[splitEnd].refStrand.strand != clusterLocations[clusterEnd].strand)
		{
			return false;
		}
		
		if (splitAlignment.position[splitEnd].position <= clusterLocations[clusterEnd].start - cMateSearchLength || splitAlignment.position[splitEnd].position >= clusterLocations[clusterEnd].end + cMateSearchLength)
		{
			return false;
		}
	}
	
	return true;
}

bool SplitScoreLessThan(const pair<CompactSplitDefinition,int>& a, const pair<CompactSplitDefinition,int>& b)
{
	return a.second < b.second;
}

void ReverseComplement(char& nucleotide)
{
	switch (nucleotide)
	{
		case 'A': nucleotide = 'T'; break;
		case 'C': nucleotide = 'G'; break;
		case 'T': nucleotide = 'A'; break;
		case 'G': nucleotide = 'C'; break;
		case 'a': nucleotide = 't'; break;
		case 'c': nucleotide = 'g'; break;
		case 't': nucleotide = 'a'; break;
		case 'g': nucleotide = 'c'; break;
	}
}

int CalculateOffset(int strand, int offset)
{
	int dir = (strand == PlusStrand) ? 1 : -1;
	return offset * dir;
}

int CalculateForwardHomology(int maxOffset, const NameIndex& refNameIndex, const Sequences& sequences, const CompactPosition& position1, const CompactPosition& position2)
{
	const char* seqPtr1 = sequences.Get(refNameIndex.Get(position1.refStrand.referenceIndex), position1.position);
	const char* seqPtr2 = sequences.Get(refNameIndex.Get(position2.refStrand.referenceIndex), position2.position);
	
	int homology = 0;
	for (int offset = 1; offset <= maxOffset; offset++)
	{
		char nt1 = *(seqPtr1 + CalculateOffset(position1.refStrand.strand, offset));
		char nt2 = *(seqPtr2 + CalculateOffset(position2.refStrand.strand, 1 - offset));
		
		if (position1.refStrand.strand != PlusStrand)
		{
			ReverseComplement(nt1);
		}
		
		if (position2.refStrand.strand != MinusStrand)
		{
			ReverseComplement(nt2);
		}
		
		if (nt1 != nt2)
		{
			break;
		}
		
		homology = offset;
	}
	
	return homology;
}

int PositionWithOffset(const CompactPosition& position, int offset)
{
	return position.position + CalculateOffset(position.refStrand.strand, offset);
}

void WriteBreakpoint(ostream& out, const NameIndex& refNameIndex, int clusterID, const CompactSplitDefinition& breakpoint, int splitCount, int offset1, int offset2)
{
	out << clusterID << "\t";
	out << refNameIndex.Get(breakpoint.position[0].refStrand.referenceIndex) << "\t";
	out << ((breakpoint.position[0].refStrand.strand == PlusStrand) ? "+" : "-") << "\t";
	out << PositionWithOffset(breakpoint.position[0], offset1) << "\t";
	out << refNameIndex.Get(breakpoint.position[1].refStrand.referenceIndex) << "\t";
	out << ((breakpoint.position[1].refStrand.strand == PlusStrand) ? "+" : "-") << "\t";
	out << PositionWithOffset(breakpoint.position[1], offset2) << "\t";
	out << breakpoint.numInserted << "\t";
	out << splitCount << endl;
}

struct ClusterSplitInfo
{
	unordered_map<CompactSplitDefinition,int> count;
	unordered_map<CompactSplitDefinition,int> score;
};

int main(int argc, char* argv[])
{
	string splitAlignmentFilename;
	string clustersFilename;
	string breakpointsFilename;
	string referenceFasta;
	
	try
	{
		TCLAP::CmdLine cmd("Breakpoint Prediction Tool");
		TCLAP::ValueArg<string> splitAlignmentFilenameArg("s","split","Split Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> breakpointsFilenameArg("b","breakpoints","Output Breakpoints Filename",true,"","integer",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		splitAlignmentFilename = splitAlignmentFilenameArg.getValue();
		clustersFilename = clustersFilenameArg.getValue();
		breakpointsFilename = breakpointsFilenameArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	const int cMaxBreakpointHomology = 25;
	
	NameIndex refNameIndex;
	
	cerr << "Reading clusters" << endl;
	
	unordered_map<int,unordered_set<int> > fragmentClusters;
	unordered_map<int,LocationVec> clusterLocations;
	
	ifstream clustersFile(clustersFilename.c_str());
	CheckFile(clustersFile, clustersFilename);
	
	ClusterReader clusterReader(clustersFile);
	
	while (clusterReader.Next())
	{
		int clusterID = clusterReader.FetchClusterID();
		
		clusterLocations[clusterID] = clusterReader.FetchLocations();
		IntegerVec fragmentIndices = clusterReader.FetchFragmentIndices();
		
		for (int fidx = 0; fidx < fragmentIndices.size(); fidx++)
		{
			for (int readEnd = 0; readEnd <= 1; readEnd++)
			{
				ReadID readID;
				readID.fragmentIndex = fragmentIndices[fidx];
				readID.readEnd = readEnd;
				
				fragmentClusters[readID.id].insert(clusterID);
			}
		}
	}
	
	cerr << "Reading split alignments" << endl;
	
	unordered_map<int,ClusterSplitInfo> clusterSplitInfo;
	
	ifstream splitAlignmentFile(splitAlignmentFilename.c_str());
	CheckFile(splitAlignmentFile, splitAlignmentFilename);
	
	int line = 1;
	StringVec fields;
	while (ReadTSV(splitAlignmentFile, fields))
	{
		if (fields.size() < 7)
		{
			cerr << "Invalid split alignments line, " << splitAlignmentFilename << ":" << line << endl;
			exit(1);
		}
		
		ReadID readID;
		readID.fragmentIndex = SAFEPARSE(int, fields[0]);
		readID.readEnd = SAFEPARSE(int, fields[1]);
		
		CompactSplitAlignment splitAlignment;
		
		splitAlignment.position[0].refStrand.referenceIndex = refNameIndex.Index(fields[2]);
		splitAlignment.position[0].refStrand.strand = (fields[3] == "+") ? PlusStrand : MinusStrand;
		splitAlignment.position[0].position = SAFEPARSE(int, fields[4]);
		
		splitAlignment.position[1].refStrand.referenceIndex = refNameIndex.Index(fields[5]);
		splitAlignment.position[1].refStrand.strand = (fields[6] == "+") ? PlusStrand : MinusStrand;
		splitAlignment.position[1].position = SAFEPARSE(int, fields[7]);
		
		splitAlignment.inserted = fields[8];
		
		splitAlignment.score = SAFEPARSE(int, fields[13]);
		
		unordered_map<int,unordered_set<int> >::const_iterator clustersIter = fragmentClusters.find(readID.id);
		if (clustersIter != fragmentClusters.end())
		{
			for (unordered_set<int>::const_iterator clusterIDIter = clustersIter->second.begin(); clusterIDIter != clustersIter->second.end(); clusterIDIter++)
			{
				int clusterID = *clusterIDIter;
				
				for (int flip = 0; flip <= 1; flip++)
				{
					if (Match(clusterLocations[clusterID], refNameIndex, splitAlignment, flip))
					{
						CompactSplitDefinition split;
						
						split.position[0] = splitAlignment.position[0];
						split.position[1] = splitAlignment.position[1];
						split.numInserted = (int)splitAlignment.inserted.size();
						
						if (flip == 1)
						{
							swap(split.position[0], split.position[1]);
						}
						
						clusterSplitInfo[clusterID].count.insert(make_pair(split, 0)).first->second++;
						clusterSplitInfo[clusterID].score.insert(make_pair(split, 0)).first->second += splitAlignment.score;
					}
				}
			}
		}
		
		line++;
	}
	
	cerr << "Reading reference fasta" << endl;
	
	Sequences referenceSequences(1000);
	referenceSequences.Read(referenceFasta);
	
	cerr << "Calculating breakpoint homology and outputting breakpoints" << endl;
	
	ofstream breakpointsFile(breakpointsFilename.c_str());
	CheckFile(breakpointsFile, breakpointsFilename);
	
	for (unordered_map<int,ClusterSplitInfo>::const_iterator splitInfoIter = clusterSplitInfo.begin(); splitInfoIter != clusterSplitInfo.end(); splitInfoIter++)
	{
		int clusterID = splitInfoIter->first;
		
		const CompactSplitDefinition& breakpoint = max_element(splitInfoIter->second.score.begin(), splitInfoIter->second.score.end(), SplitScoreLessThan)->first;
		int splitCount = splitInfoIter->second.count.find(breakpoint)->second;
		
		int maxOffset1 = CalculateForwardHomology(cMaxBreakpointHomology, refNameIndex, referenceSequences, breakpoint.position[0], breakpoint.position[1]);
		int maxOffset2 = CalculateForwardHomology(cMaxBreakpointHomology, refNameIndex, referenceSequences, breakpoint.position[1], breakpoint.position[0]);
		
		if (maxOffset1 + maxOffset2 + 1 > cMaxBreakpointHomology)
		{
			continue;
		}
		
		for (int offset = -maxOffset2; offset <= maxOffset1; offset++)
		{
			WriteBreakpoint(breakpointsFile, refNameIndex, clusterID, breakpoint, splitCount, offset, -offset);
		}
	}
}

