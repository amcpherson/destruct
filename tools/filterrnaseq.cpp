/*
 *  filterrnaseq.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "DebugCheck.h"
#include "Indexer.h"
#include "AlignmentStream.h"
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


bool LocationLessThan(const CompactLocation& loc1, const CompactLocation& loc2)
{
	if (loc1.refStrand.id < loc2.refStrand.id)
	{
		return true;
	}
	else if (loc1.refStrand.id > loc2.refStrand.id)
	{
		return false;
	}
	
	return loc1.region.start < loc2.region.start;
}

bool Overlap(const CompactLocation& r1, const CompactLocation& r2)
{
	return (r1.refStrand.id == r2.refStrand.id && !(r1.region.end < r2.region.start || r2.region.end < r1.region.start));
}

void Extend(CompactLocation& r1, const CompactLocation& r2)
{
	r1.region.start = min(r1.region.start, r2.region.start);
	r1.region.end = max(r1.region.end, r2.region.end);
}

void Merge(CompactLocationVec& locations)
{
	if (locations.empty())
	{
		return;
	}
	
	sort(locations.begin(), locations.end(), LocationLessThan);
	
	CompactLocationVec merged;
	merged.push_back(*locations.begin());
	for (CompactLocationVecConstIter locationIter = locations.begin(); locationIter != locations.end(); locationIter++)
	{
		if (Overlap(merged.back(), *locationIter))
		{
			Extend(merged.back(), *locationIter);
		}
		else
		{
			merged.push_back(*locationIter);
		}
	}
	
	locations = merged;
}

bool Overlaps(CompactLocationVec& locations1, CompactLocationVec& locations2)
{
	CompactLocationVec bothLocations;
	bothLocations.insert(bothLocations.end(), locations1.begin(), locations1.end());
	bothLocations.insert(bothLocations.end(), locations2.begin(), locations2.end());
	
	Merge(bothLocations);
	
	return bothLocations.size() < locations1.size() + locations2.size();
}

bool Intersects(const set<string>& set1, const set<string>& set2)
{
	set<string>::const_iterator iter1 = set1.begin();
	set<string>::const_iterator iter2 = set2.begin();
	
	while (iter1 != set1.end() && iter2 != set2.end())
	{
		if (*iter1 < *iter2)
		{
			iter1++;
		}
		else if (*iter2 < *iter1)
		{
			iter2++;
		}
		else
		{
			return true;
		}
	}
	
	return false;
}

bool IsFiltered(const RawAlignmentVec& alignments, const ExonRegions& exonRegions, int maxConcordantDistance, int maxGenomicLoci)
{
	NameIndex chromosomeNames;
	
	set<string> overlappedGenes[2];
	CompactLocationVec alignmentLocations[2];
	
	for (RawAlignmentVec::const_iterator alignmentIter = alignments.begin(); alignmentIter != alignments.end(); alignmentIter++)
	{
		string chromosome = alignmentIter->reference;
		int strand = alignmentIter->strand;
		Region region = alignmentIter->region;
		
		string alignGene;
		string alignTranscript;
		if (ParseTranscriptID(alignmentIter->reference, alignGene, alignTranscript) && exonRegions.IsTranscript(alignTranscript))
		{
			exonRegions.RemapTranscriptToGenome(alignTranscript, alignmentIter->strand, alignmentIter->region.start, chromosome, strand, region.start);
			exonRegions.RemapTranscriptToGenome(alignTranscript, alignmentIter->strand, alignmentIter->region.end, chromosome, strand, region.end);
			
			if (region.start > region.end)
			{
				swap(region.start, region.end);
			}
		}
		
		Region expandedRegion = region;
		expandedRegion.start -= maxConcordantDistance;
		expandedRegion.end += maxConcordantDistance;
		
		exonRegions.GetRegionGenes(chromosome, expandedRegion, overlappedGenes[alignmentIter->readEnd]);
		
		Region halfExpandedRegion = region;
		halfExpandedRegion.start -= maxConcordantDistance / 2;
		halfExpandedRegion.end += maxConcordantDistance / 2;
		
		CompactLocation location;
		location.refStrand.referenceIndex = chromosomeNames.Index(chromosome);
		location.refStrand.strand = 0;
		location.region = halfExpandedRegion;
		
		alignmentLocations[alignmentIter->readEnd].push_back(location);
	}
	
	if (Intersects(overlappedGenes[0], overlappedGenes[1]))
	{
		return true;
	}
		
	for (int readEnd = 0; readEnd <= 1; readEnd++)
	{
		Merge(alignmentLocations[readEnd]);
		
		if (alignmentLocations[readEnd].size() > maxGenomicLoci)
		{
			return true;
		}
	}
	
	if (Overlaps(alignmentLocations[0], alignmentLocations[1]))
	{
		return true;
	}
	
	return false;
}

void Write(ostream& out, const RawAlignmentVec& alignments)
{
	for (RawAlignmentVecConstIter alignmentIter = alignments.begin(); alignmentIter != alignments.end(); alignmentIter++)
	{
		out << alignmentIter->line << endl;
	}
}

int main(int argc, char* argv[])
{
	string samFilename;
	string alignFilename;
	string exonRegionsFilename;
	int maxConcordantDistance;
	int maxGenomicLoci;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Pair Filtering Tool For RNA-seq");
		TCLAP::ValueArg<string> samFilenameArg("s","sam","Read Sorted Sam Filename",false,"","string");
		TCLAP::ValueArg<string> alignFilenameArg("a","align","Read Sorted Compact Alignments Filename",true,"","string");
		TCLAP::ValueArg<string> exonRegionsFilenameArg("e","exons","Exon Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> maxConcordantDistanceArg("d","dist","Maximum Distance For Concordancy",true,0,"integer",cmd);
		TCLAP::ValueArg<int> maxGenomicLociArg("l","loci","Max Number of Genomic Loci Per Read",true,0,"integer",cmd);
		
		vector<TCLAP::Arg*> alignmentArgs;
		alignmentArgs.push_back(&samFilenameArg);
		alignmentArgs.push_back(&alignFilenameArg);
		cmd.xorAdd(alignmentArgs);
		
		cmd.parse(argc,argv);
		
		samFilename = samFilenameArg.getValue();
		alignFilename = alignFilenameArg.getValue();
		exonRegionsFilename = exonRegionsFilenameArg.getValue();
		maxConcordantDistance = maxConcordantDistanceArg.getValue();
		maxGenomicLoci = maxGenomicLociArg.getValue();
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
	
	ExonRegions exonRegions;
	exonRegions.Read(exonRegionsFilename);
	
	RawAlignmentVec alignments;
	while (fragmentAlignmentStream.GetNextAlignments(alignments))
	{
		if (IsFiltered(alignments, exonRegions, maxConcordantDistance, maxGenomicLoci))
		{
			continue;
		}
		
		Write(cout, alignments);
	}
}

