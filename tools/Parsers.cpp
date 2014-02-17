/*
 *  Parsers.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-09-02.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "Parsers.h"
#include "Common.h"
#include "DebugCheck.h"
#include "Indexer.h"

#include <map>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_set.hpp>

using namespace boost;
using namespace std;


class ClusterReaderRemapFragment : public ClusterReader
{
public:
	ClusterReaderRemapFragment(istream& clustersFile) : ClusterReader(clustersFile)
	{
		mNullLibIndex = mLibIndex.Index("");
	}
	
	bool Next()
	{
		if (!ClusterReader::Next())
		{
			return false;
		}
		
		mRemappedRowFragmentIndices.clear();
		mRemappedFragmentIndices.clear();
		for (StringTableConstIter rowIter = RowData().begin(); rowIter != RowData().end(); rowIter++)
		{
			const StringVec& fields = *rowIter;
			
			int fragmentIndex = SAFEPARSE(int, fields[2]);
			
			int libraryIndex = mNullLibIndex;
			if (fields.size() >= 12)
			{
				libraryIndex = mLibIndex.Index(fields[11]);
			}
			
			int remappedFragmentIndex = mLibFragmentIndex.Index(pair<int,int>(fragmentIndex, libraryIndex));
			
			mRemappedRowFragmentIndices.push_back(remappedFragmentIndex);
			
			if (fields[1] == "0")
			{
				mRemappedFragmentIndices.push_back(remappedFragmentIndex);
			}
		}
		
		return true;
	}
	
	const IntegerVec& RemappedRowFragmentIndices()
	{
		return mRemappedRowFragmentIndices;
	}
	
	const IntegerVec& RemappedFragmentIndices()
	{
		return mRemappedFragmentIndices;
	}
	
private:
	Indexer<string> mLibIndex;
	int mNullLibIndex;
	Indexer<pair<int,int> > mLibFragmentIndex;
	IntegerVec mRemappedRowFragmentIndices;
	IntegerVec mRemappedFragmentIndices;
};

void ClusterMembership::Read(IntegerVecMap& membership)
{
	ifstream clustersFile(mClustersFilename.c_str());
	CheckFile(clustersFile, mClustersFilename);
	
	ClusterReaderRemapFragment clusterReader(clustersFile);
	while (clusterReader.Next())
	{
		int clusterID = clusterReader.FetchClusterID();
		membership[clusterID] = clusterReader.RemappedFragmentIndices();
	}
	
	clustersFile.close();
}

void ClusterMembership::Read(IntegerVecMap& membership, IntegerMap& distances)
{
	ifstream clustersFile(mClustersFilename.c_str());
	CheckFile(clustersFile, mClustersFilename);
	
	ClusterReaderRemapFragment clusterReader(clustersFile);
	while (clusterReader.Next())
	{
		int clusterID = clusterReader.FetchClusterID();
		
		membership[clusterID] = clusterReader.RemappedFragmentIndices();
		
		LocationVec clusterLocations = clusterReader.FetchLocations();
		
		if (clusterLocations[0].refName != clusterLocations[1].refName)
		{
			continue;
		}
		
		int positions[2];
		for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
		{
			positions[clusterEnd] = (clusterLocations[clusterEnd].strand == PlusStrand) ? clusterLocations[clusterEnd].end : clusterLocations[clusterEnd].start;
		}
		
		distances[clusterID] = abs(positions[0] - positions[1]);
	}
	
	clustersFile.close();
}

void ClusterMembership::Write(const string& filename, const IntegerVecMap& membership)
{
	ifstream clustersFile(mClustersFilename.c_str());
	CheckFile(clustersFile, mClustersFilename);
	
	ofstream outputFile(filename.c_str());
	CheckFile(outputFile, filename);
	
	ClusterReaderRemapFragment clusterReader(clustersFile);
	while (clusterReader.Next())
	{
		int clusterID = clusterReader.FetchClusterID();
		
		IntegerVecMapConstIter clusterIter = membership.find(clusterID);
		
		if (clusterIter == membership.end())
		{
			continue;
		}
		
		IntegerSet fragmentFilter(clusterIter->second.begin(), clusterIter->second.end());
		
		const StringTable& rowData = clusterReader.RowData();
		IntegerVec rowFragmentIndices = clusterReader.RemappedRowFragmentIndices();
		
		for (int rowIndex = 0; rowIndex < (int)rowData.size(); rowIndex++)
		{
			int fragmentIndex = rowFragmentIndices[rowIndex];
			
			if (fragmentFilter.find(fragmentIndex) != fragmentFilter.end())
			{
				WriteTSV(outputFile, rowData[rowIndex]);
			}
		}
	}
	
	clustersFile.close();
}

void ClusterMembership::Write(const string& filename, const IntegerSet& clusterIDs)
{
	ifstream clustersFile(mClustersFilename.c_str());
	CheckFile(clustersFile, mClustersFilename);
	
	ofstream outputFile(filename.c_str());
	CheckFile(outputFile, filename);
	
	ClusterReader clusterReader(clustersFile);
	while (clusterReader.Next())
	{
		int clusterID = clusterReader.FetchClusterID();
		
		if (clusterIDs.find(clusterID) == clusterIDs.end())
		{
			continue;
		}
		
		const StringTable& rowData = clusterReader.RowData();
		
		for (int rowIndex = 0; rowIndex < (int)rowData.size(); rowIndex++)
		{
			WriteTSV(outputFile, rowData[rowIndex]);
		}
	}
	
	clustersFile.close();
}

void ReadClusters(const string& clustersFilename, CompactLocationVecMap& clusters, NameIndex& references)
{
	// Open clusters file
	ifstream clustersFile(clustersFilename.c_str());
	CheckFile(clustersFile, clustersFilename);
	
	ClusterReader clusterReader(clustersFile);
	while (clusterReader.Next())
	{
		int clusterID = clusterReader.FetchClusterID();
		LocationVec clusterLocations = clusterReader.FetchLocations();
		
		for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
		{
			ClusterEndID clusterEndID;
			clusterEndID.clusterID = clusterID;
			clusterEndID.clusterEnd = clusterEnd;
			
			CompactLocation location;
			location.refStrand.referenceIndex = references.Index(clusterLocations[clusterEnd].refName);
			location.refStrand.strand = clusterLocations[clusterEnd].strand;
			location.region.start = clusterLocations[clusterEnd].start;
			location.region.end = clusterLocations[clusterEnd].end;
			
			clusters[clusterID].push_back(location);
		}
	}
	
	clustersFile.close();
}

void ReadCorroboration(const string& corroborationFilename, IntegerVecPairVec& corroboration)
{
	// Open corroboration file
	ifstream corroborationFile(corroborationFilename.c_str());
	if (!corroborationFile)
	{
		cerr << "Error: unable to read from corroboration file " << corroborationFilename << endl;		
		exit(1);
	}
	
	// Parse file contents
	string line;
	int lineNumber = 0;
	while (getline(corroborationFile, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty corroboration line " << lineNumber << " of " << corroborationFilename << endl;
			exit(1);
		}
		
		vector<string> corroborationFields;
		split(corroborationFields, line, is_any_of("\t"));
		
		if (corroborationFields.size() < 1)
		{
			cerr << "Error: Format error for corroboration line " << lineNumber << " of " << corroborationFilename << endl;
			exit(1);
		}
		
		IntegerVec rnaClusterIDs;
		if (corroborationFields[0] != "")
		{
			int rnaClusterID = SAFEPARSE(int, corroborationFields[0]);
			rnaClusterIDs.push_back(rnaClusterID);
		}
		
		IntegerVec dnaClusterIDs;
		for (int fieldIndex = 1; fieldIndex < (int)corroborationFields.size(); fieldIndex++)
		{
			dnaClusterIDs.push_back(SAFEPARSE(int, corroborationFields[fieldIndex]));
		}
		
		corroboration.push_back(IntegerVecPair());
		swap(corroboration.back().first, rnaClusterIDs);
		swap(corroboration.back().second, dnaClusterIDs);		
	}
	
	corroborationFile.close();
}

void IntepretAlignString(const string& alignString, Location& alignRegion)
{
	string::size_type colonDividerPos = alignString.find_first_of(":");
	if (colonDividerPos == string::npos || colonDividerPos == 0)
	{
		cerr << "Error: Unable to interpret strand for " << alignString << endl;
		exit(1);
	}
	
	char strand = alignString[colonDividerPos - 1];
	
	if (strand == '+')
	{
		alignRegion.strand = PlusStrand;
	}
	else if (strand == '-')
	{
		alignRegion.strand = MinusStrand;
	}
	else
	{
		cerr << "Error: Unable to interpret strand for " << alignString << endl;
		exit(1);
	}
	
	vector<string> alignFields;
	split(alignFields, alignString, is_any_of("+-:"));
	
	if (alignFields.size() != 4)
	{
		cerr << "Error: Unable to interpret alignment string " << alignString << endl;
		exit(1);
	}
	
	alignRegion.refName = alignFields[0];
	alignRegion.start = SAFEPARSE(int, alignFields[2]);
	alignRegion.end = SAFEPARSE(int, alignFields[3]);
}

void ReadAlignRegionPairs(const string& filename, LocationVecMap& alignRegionPairs)
{
	ifstream alignRegionPairsFile(filename.c_str());
	if (!alignRegionPairsFile.good())
	{
		cerr << "Error: Unable to open align region pairs file " << filename << endl;
		exit(1);
	}
	
	string line;
	int lineNumber = 0;
	
	while (getline(alignRegionPairsFile, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			continue;
		}
		
		vector<string> alignRegionFields;
		split(alignRegionFields, line, is_any_of("\t"));
		
		if (alignRegionFields.size() < 5)
		{
			continue;
		}
		
		int pairID = SAFEPARSE(int, alignRegionFields[0]);
		int pairEnd = SAFEPARSE(int, alignRegionFields[1]);
		
		DebugCheck(pairEnd == 0 || pairEnd == 1);
		
		Location alignRegion;
		alignRegion.refName = alignRegionFields[2];
		alignRegion.strand = InterpretStrand(alignRegionFields[3]);
		alignRegion.start = SAFEPARSE(int, alignRegionFields[4]);
		alignRegion.end = SAFEPARSE(int, alignRegionFields[5]);
		
		alignRegionPairs[pairID].resize(2);
		alignRegionPairs[pairID][pairEnd] = alignRegion;
	}
	
	alignRegionPairsFile.close();
}

void ReadDoubleMap(const string& filename, DoubleMap& values)
{
	ifstream file(filename.c_str());
	if (!file.good())
	{
		cerr << "Error: Unable to open file " << filename << endl;
		exit(1);
	}
	
	string line;
	int lineNumber = 0;
	
	while (getline(file, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			continue;
		}
		
		vector<string> fields;
		split(fields, line, is_any_of("\t"));
		
		if (fields.size() < 2)
		{
			continue;
		}
		
		int id = SAFEPARSE(int, fields[0]);
		double value = SAFEPARSE(long double, fields[1]);
		
		values[id] = value;
	}
	
	file.close();
}

void ReadDoubleVec(const string& filename, DoubleVec& values)
{
	ifstream file(filename.c_str());
	if (!file.good())
	{
		cerr << "Error: Unable to open file " << filename << endl;
		exit(1);
	}
	
	string line;
	int lineNumber = 0;
	
	while (getline(file, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			continue;
		}
		
		values.push_back(SAFEPARSE(double, line));
	}
	
	file.close();
}

void ReadClusterProbabilities(const string& filename, DoubleMap& probabilities)
{
	ifstream file(filename.c_str());
	if (!file.good())
	{
		cerr << "Error: Unable to open file " << filename << endl;
		exit(1);
	}
	
	string line;
	int lineNumber = 0;
	
	while (getline(file, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			continue;
		}
		
		vector<string> fields;
		split(fields, line, is_any_of("\t"));
		
		if (fields.size() < 4)
		{
			continue;
		}
		
		int id = SAFEPARSE(int, fields[0]);
		double alignProb = SAFEPARSE(long double, fields[1]);
		double chimericProb = SAFEPARSE(long double, fields[2]);
		double validProb = SAFEPARSE(long double, fields[3]);
		
		probabilities[id] = alignProb * chimericProb * validProb;
	}
	
	file.close();
}

void ReadIntegerVecMap(const string& filename, IntegerVecMap& values)
{
	ifstream file(filename.c_str());
	if (!file.good())
	{
		cerr << "Error: Unable to open file " << filename << endl;
		exit(1);
	}
	
	string line;
	int lineNumber = 0;
	
	while (getline(file, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			continue;
		}
		
		vector<string> fields;
		split(fields, line, is_any_of("\t"));
		
		if (fields.size() < 2)
		{
			continue;
		}
		
		int id = SAFEPARSE(int, fields[0]);
		int value = SAFEPARSE(int, fields[1]);
		
		values[id].push_back(value);
	}
	
	file.close();
}

bool ReadTSV(istream& file, StringVec& fields)
{
	string line;
	if (!getline(file, line))
	{
		return false;
	}
	
	split(fields, line, is_any_of("\t"));
	
	return true;
}

void WriteTSV(ostream& file, const StringVec& fields)
{
	for (StringVecConstIter fieldIter = fields.begin(); fieldIter != fields.end(); fieldIter++)
	{
		if (fieldIter != fields.begin())
		{
			file << "\t";
		}
		file << *fieldIter;
	}
	file << "\n";
}

void ReadStringPairs(const string& filename, vector<pair<string,string> >& values)
{
	ifstream file(filename.c_str());
	if (!file.good())
	{
		cerr << "Error: Unable to open file " << filename << endl;
		exit(1);
	}
	
	string line;
	int lineNumber = 0;
	
	while (getline(file, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			continue;
		}
		
		vector<string> fields;
		split(fields, line, is_any_of("\t"));
		
		if (fields.size() < 2)
		{
			continue;
		}
		
		values.push_back(make_pair(fields[0],fields[1]));
	}
	
	file.close();
}

ClusterReader::ClusterReader(istream& clustersFile) : mClustersFile(clustersFile), mLineNumber(0)
{
	if (!mClustersFile)
	{
		cerr << "Error: unable to read clusters" << endl;		
		exit(1);
	}
	
	mGetResult = ReadTSV(mClustersFile, mRow);
	
	if (mGetResult)
	{
		mLineNumber++;
		mClusterID = mRow[0];
	}
}

bool ClusterReader::Next()
{
	if (!mGetResult)
	{
		return false;
	}
	
	mClusterRows.clear();
	while (mGetResult && mRow[0] == mClusterID)
	{
		if (mRow.size() < 8)
		{
			cerr << "Error: Invalid clusters line " << mLineNumber << endl;
			exit(1);
		}
		
		mClusterRows.push_back(mRow);
		mGetResult = ReadTSV(mClustersFile, mRow);
		mLineNumber++;
	}
	
	if (mGetResult)
	{
		mClusterID = mRow[0];
	}
	
	return true;
}

int ClusterReader::FetchClusterID() const
{
	return SAFEPARSE(int, mClusterRows.front()[0]);
}

LocationVec ClusterReader::FetchLocations() const
{
	LocationVec locations(2);
	
	bool initialized[2] = {false, false};
	for (StringTableConstIter rowIter = mClusterRows.begin(); rowIter != mClusterRows.end(); rowIter++)
	{
		const StringVec& fields = *rowIter;
		
		int clusterEnd = SAFEPARSE(int, fields[1]);
		
		const string& refName = fields[4];
		int strand = InterpretStrand(fields[5]);
		
		int start = SAFEPARSE(int, fields[6]);
		int end = SAFEPARSE(int, fields[7]);
		
		locations[clusterEnd].refName = refName;
		locations[clusterEnd].strand = strand;
		
		if (!initialized[clusterEnd])
		{
			locations[clusterEnd].start = start;
			locations[clusterEnd].end = end;
		}
		else
		{
			locations[clusterEnd].start = min(locations[clusterEnd].start, start);
			locations[clusterEnd].end = max(locations[clusterEnd].end, end);
		}
		
		initialized[clusterEnd] = true;
	}
	
	if (!initialized[0] || !initialized[1])
	{
		cerr << "Error: cluster " << mClusterRows.front()[0] << " is invalid" << endl;		
		exit(1);
	}
	
	return locations;
}

IntegerVec ClusterReader::FetchFragmentIndices() const
{
	IntegerVec fragmentIndices;
	
	for (StringTableConstIter rowIter = mClusterRows.begin(); rowIter != mClusterRows.end(); rowIter++)
	{
		const StringVec& fields = *rowIter;
		
		if (fields[1] != "0")
		{
			continue;
		}
		
		fragmentIndices.push_back(SAFEPARSE(int, fields[2]));
	}
	
	return fragmentIndices;
}

DoubleTable ClusterReader::FetchProbabilities(int fieldIndex) const
{
	IntegerVec fragmentIndices;
	DoubleMap readProbabilities[2];
	
	for (StringTableConstIter rowIter = mClusterRows.begin(); rowIter != mClusterRows.end(); rowIter++)
	{
		const StringVec& fields = *rowIter;
		
		int fragmentIndex = SAFEPARSE(int, fields[2]);
		int readEnd = SAFEPARSE(int, fields[3]);
		
		readProbabilities[readEnd][fragmentIndex] = SAFEPARSE(double, fields[fieldIndex]);
		
		if (fields[1] == "0")
		{
			fragmentIndices.push_back(fragmentIndex);
		}
	}
	
	DoubleTable fragmentProbabilities;
	
	for (IntegerVecConstIter fragmentIter = fragmentIndices.begin(); fragmentIter != fragmentIndices.end(); fragmentIter++)
	{
		fragmentProbabilities.push_back(DoubleVec(2));
		fragmentProbabilities.back()[0] = readProbabilities[0].find(*fragmentIter)->second;
		fragmentProbabilities.back()[1] = readProbabilities[1].find(*fragmentIter)->second;
	}
	
	return fragmentProbabilities;
}

DoubleTable ClusterReader::FetchAlignProbabilities() const
{
	return FetchProbabilities(8);
}

DoubleTable ClusterReader::FetchChimericProbabilities() const
{
	return FetchProbabilities(9);
}

DoubleTable ClusterReader::FetchValidProbabilities() const
{
	return FetchProbabilities(10);
}

StringVec ClusterReader::FetchLibraryNames() const
{
	StringVec libraryNames;
	
	for (StringTableConstIter rowIter = mClusterRows.begin(); rowIter != mClusterRows.end(); rowIter++)
	{
		const StringVec& fields = *rowIter;
		
		if (fields[1] != "0")
		{
			continue;
		}
		
		if (fields.size() >= 12)
		{
			libraryNames.push_back(fields[11]);
		}
		else
		{
			libraryNames.push_back("");
		}
	}
	
	return libraryNames;
}

const StringTable& ClusterReader::RowData() const
{
	return mClusterRows;
}

int ClusterReader::FetchFragmentCount() const
{
	return (int)mClusterRows.size() / 2;
}

void ReadFAI(const string& faiFilename, vector<string>& referenceNames, vector<long>& referenceLengths)
{
	ifstream faiFile(faiFilename.c_str());
	if (!faiFile)
	{
		cerr << "Error: unable to fai file" << endl;
		exit(1);
	}
	
	string line;
	int lineNumber = 0;
	while (getline(faiFile, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty fai line " << lineNumber << " of " << faiFilename << endl;
			exit(1);
		}
		
		vector<string> faiFields;
		split(faiFields, line, is_any_of("\t"));
		
		if (faiFields.size() < 2)
		{
			cerr << "Error: Format error for fai line " << lineNumber << " of " << faiFilename << endl;
			exit(1);
		}
		
		string referenceName = faiFields[0];
		long referenceLength = SAFEPARSE(long, faiFields[1]);
		
		referenceNames.push_back(referenceName);
		referenceLengths.push_back(referenceLength);
	}
	
	faiFile.close();
}


