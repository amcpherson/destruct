/*
 *  AlignmentRecord.cpp
 *
 */

#include "AlignmentRecord.h"
#include "DebugCheck.h"

#include <ostream>

int SpanningAlignmentRecord::GetOuterPosition() const
{
    DebugCheck(strand == "+" || strand == "-");

    if (strand == "+")
    {
        return start;
    }
    else
    {
        return end;
    }
}

std::ostream & operator<<(std::ostream &os, const SpanningAlignmentRecord& record)
{
    os << record.libID << "\t";
    os << record.readID << "\t";
    os << record.readEnd << "\t";
    os << record.alignID << "\t";
    os << record.chromosome << "\t";
    os << record.strand << "\t";
    os << record.start << "\t";
    os << record.end << "\t";
    os << record.score << std::endl;
    return os;
}

std::istream & operator>>(std::istream &is, SpanningAlignmentRecord& record)
{
    is >> record.libID;
    is >> record.readID;
    is >> record.readEnd;
    is >> record.alignID;
    is >> record.chromosome;
    is >> record.strand;
    is >> record.start;
    is >> record.end;
    is >> record.score;
    return is;
}

std::ostream & operator<<(std::ostream &os, const SplitAlignmentRecord& record)
{
    os << record.libID << "\t";
    os << record.readID << "\t";
    os << record.readEnd << "\t";
    for (int readEnd = 0; readEnd < 2; readEnd++)
    {
        os << record.alignID[readEnd] << "\t";
        os << record.chromosome[readEnd] << "\t";
        os << record.strand[readEnd] << "\t";
        os << record.position[readEnd] << "\t";
    }
    os << record.inserted << "\t";
    os << record.score << std::endl;
    return os;
}

std::istream & operator>>(std::istream &is, SplitAlignmentRecord& record)
{
    is >> record.libID;
    is >> record.readID;
    is >> record.readEnd;
    for (int readEnd = 0; readEnd < 2; readEnd++)
    {
        is >> record.alignID[readEnd];
        is >> record.chromosome[readEnd];
        is >> record.strand[readEnd];
        is >> record.position[readEnd];
    }
    is >> record.inserted;
    is >> record.score;
    return is;
}

AlignmentPairKey SplitAlignmentRecord::GetAlignmentPairKey() const
{
    AlignmentPairKey alignPairKey;
    alignPairKey.libID = libID;
    alignPairKey.readID = readID;
    for (int readEnd = 0; readEnd < 2; readEnd++)
    {
        alignPairKey.alignID[readEnd] = alignID[readEnd];
    }
    return alignPairKey;
}

std::ostream & operator<<(std::ostream &os, const ClusterMemberRecord& record)
{
    os << record.clusterID << "\t";
    os << record.clusterEnd << "\t";
    os << record.libID << "\t";
    os << record.readID << "\t";
    os << record.readEnd << "\t";
    os << record.alignID << std::endl;
    return os;
}

std::istream & operator>>(std::istream &is, ClusterMemberRecord& record)
{
    is >> record.clusterID;
    is >> record.clusterEnd;
    is >> record.libID;
    is >> record.readID;
    is >> record.readEnd;
    is >> record.alignID;
    return is;
}

