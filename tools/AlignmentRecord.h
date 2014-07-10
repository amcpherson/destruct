/*
 *  AlignmentRecord.h
 *
 */

#ifndef ALIGNMENTRECORD_H_
#define ALIGNMENTRECORD_H_

#include "Common.h"

#include <ostream>


struct AlignmentPairKey
{
    int libID;
    int readID;
    int alignID[2];
};

struct SpanningAlignmentRecord
{
    int libID;
    int readID;
    int readEnd;
    int alignID;
    string chromosome;
    string strand;
    int start;
    int end;
    int score;

    int GetOuterPosition() const;
};

std::ostream & operator<<(std::ostream &os, const SpanningAlignmentRecord& record);

std::istream & operator>>(std::istream &is, SpanningAlignmentRecord& record);

struct SplitAlignmentRecord
{
    int libID;
    int readID;
    int readEnd;
    int alignID[2];
    string chromosome[2];
    string strand[2];
    int position[2];
    string inserted;
    int score;

    AlignmentPairKey GetAlignmentPairKey() const;
};

std::ostream & operator<<(std::ostream &os, const SplitAlignmentRecord& record);

std::istream & operator>>(std::istream &is, SplitAlignmentRecord& record);

struct ClusterMemberRecord
{
    int clusterID;
    int clusterEnd;
    int libID;
    int readID;
    int readEnd;
    int alignID;
};

std::ostream & operator<<(std::ostream &os, const ClusterMemberRecord& record);

std::istream & operator>>(std::istream &is, ClusterMemberRecord& record);

template<typename TRecordType>
std::ostream & operator<<(std::ostream &os, const vector<TRecordType>& records)
{
    for (typename vector<TRecordType>::const_iterator recordIter = records.begin(); recordIter != records.end(); recordIter++)
    {
        os << *recordIter;
    }

    return os;
}

template<typename TRecordType>
bool ReadEqual(const TRecordType& a, const TRecordType& b)
{
    return (a.libID == b.libID &&
            a.readID == b.readID);
}

template<typename TRecordType>
bool ClusterReadEqual(const TRecordType& a, const TRecordType& b)
{
    return (a.clusterID == b.clusterID &&
            a.libID == b.libID &&
            a.readID == b.readID);
}

template<typename TRecordType>
class GroupedRecordsStream
{
public:
    GroupedRecordsStream(std::istream& is) : mStream(is)
    {
        mGood = (is >> mNextRecord);
    }
    
    template<typename TEqualFunc>
    bool Next(vector<TRecordType>& records, TEqualFunc eq)
    {
        if (!mGood)
        {
            return false;
        }
        
        records.clear();
        records.push_back(mNextRecord);
        
        while ((mGood = (mStream >> mNextRecord)))
        {
            if (eq(records.front(), mNextRecord))
            {
                break;
            }
            else
            {
                records.push_back(mNextRecord);
            }
        }
        
        return true;
    }
    
protected:
    std::istream& mStream;
    TRecordType mNextRecord;
    bool mGood;
};


#endif
