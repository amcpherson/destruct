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
    os << record.alignID1 << "\t";
    os << record.alignID2 << "\t";
    os << record.chromosome1 << "\t";
    os << record.strand1 << "\t";
    os << record.position1 << "\t";
    os << record.chromosome2 << "\t";
    os << record.strand2 << "\t";
    os << record.position2 << "\t";
    os << record.score << std::endl;
    return os;
}

std::istream & operator>>(std::istream &is, SplitAlignmentRecord& record)
{
    is >> record.libID;
    is >> record.readID;
    is >> record.readEnd;
    is >> record.alignID1;
    is >> record.alignID2;
    is >> record.chromosome1;
    is >> record.strand1;
    is >> record.position1;
    is >> record.chromosome2;
    is >> record.strand2;
    is >> record.position2;
    is >> record.score;
    return is;
}

