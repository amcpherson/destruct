/*
 *  RegionDB.h
 *
 *  Created by Andrew McPherson on 12-01-12.
 *
 */

#ifndef REGIONDB_H_
#define REGIONDB_H_

#include "Common.h"
#include "DebugCheck.h"
#include "Indexer.h"

#include <string>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;


struct CompactSimpleRegion
{
	string ref;
	int start;
	int end;
};

typedef vector<CompactSimpleRegion> CompactSimpleRegionVec;

class RegionDB
{
public:
	explicit RegionDB(int binSpacing = 10000) : mBinSpacing(binSpacing) {}
	
	void Add(const CompactSimpleRegion& region);
	void Add(const CompactSimpleRegionVec& regions);
	void Add(const string& regionsFilename);
	
	bool Contained(string ref, int start, int end) const;
	bool Overlapped(string ref, int start, int end) const;
	
	void Clear();
	
private:
	typedef pair<string,int> RefBin;
	typedef pair<int,int> StartEnd;
	typedef vector<StartEnd> StartEndVec;
	
	inline bool Contains(const StartEnd& region1, const StartEnd& region2) const
	{
		return (region2.first >= region1.first && region2.second <= region1.second);
	}
	
	inline bool Overlaps(const StartEnd& region1, const StartEnd& region2) const
	{
		return !(region1.second < region2.first || region2.second < region1.first);
	}
		
	int mBinSpacing;
	unordered_map<RefBin,StartEndVec> mBinned;
};

#endif

