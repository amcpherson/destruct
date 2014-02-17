/*
 *  Parsers.h
 */

#ifndef PARSERS_H_
#define PARSERS_H_

#include "Common.h"
#include "Indexer.h"

#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_set.hpp>

using namespace boost;
using namespace std;


void ReadClusters(const string& clustersFilename, CompactLocationVecMap& clusters, NameIndex& references);
void ReadCorroboration(const string& corroborationFilename, IntegerVecPairVec& corroboration);
void IntepretAlignString(const string& alignString, Location& alignRegion);
void ReadAlignRegionPairs(const string& filename, LocationVecMap& alignRegionPairs);
void ReadDoubleMap(const string& filename, DoubleMap& values);
void ReadDoubleVec(const string& filename, DoubleVec& values);
void ReadClusterProbabilities(const string& filename, DoubleMap& probabilities);
void ReadIntegerVecMap(const string& filename, IntegerVecMap& values);
bool ReadTSV(istream& file, StringVec& fields);
void WriteTSV(ostream& file, const StringVec& fields);
void ReadStringPairs(const string& filename, vector<pair<string,string> >& values);
void ReadFAI(const string& faiFilename, vector<string>& referenceNames, vector<long>& referenceLengths);

class ClusterReader
{
public:
	ClusterReader(istream& clustersFile);
	
	virtual bool Next();
	
	int FetchClusterID() const;
	LocationVec FetchLocations() const;
	IntegerVec FetchFragmentIndices() const;
	DoubleTable FetchAlignProbabilities() const;
	DoubleTable FetchChimericProbabilities() const;
	DoubleTable FetchValidProbabilities() const;
	StringVec FetchLibraryNames() const;
	int FetchFragmentCount() const;
	
	const StringTable& RowData() const;
	
private:
	DoubleTable FetchProbabilities(int fieldIndex) const;
	
	istream& mClustersFile;
	
	bool mGetResult;
	StringVec mRow;
	string mClusterID;
	StringTable mClusterRows;
	int mLineNumber;
};

class ClusterMembership
{
public:
	explicit ClusterMembership(const string& filename) : mClustersFilename(filename) {}
	
	void Read(IntegerVecMap& membership);
	void Read(IntegerVecMap& membership, IntegerMap& distances);
	void Write(const string& filename, const IntegerVecMap& membership);
	void Write(const string& filename, const IntegerSet& membership);

private:
	const string& mClustersFilename;
};

#endif


