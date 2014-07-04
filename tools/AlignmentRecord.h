/*
 *  AlignmentRecord.h
 *
 */

#ifndef ALIGNMENTRECORD_H_
#define ALIGNMENTRECORD_H_

#include "Common.h"

#include <ostream>

struct SpanningAlignmentRecord
{
    int libID;
    int readID;
    int readEnd;
    int alignID;
    string chromosome;
    char strand;
    int start;
    int end;
    int score;
};

std::ostream & operator<<(std::ostream &os, const SpanningAlignmentRecord& record)
{
    os << record.libID << "\t";
    os << record.readID << "\t";
    os << record.readEnd << "\t";
    os << record.alignID << "\t";
    os << record.chromosome << "\t";
    os << ((record.strand == PlusStrand) ? '+' : '-') << "\t";
    os << record.start << "\t";
    os << record.end << "\t";
    os << record.score << std::endl;
    return os;
}

struct SplitAlignmentRecord
{
    int libID;
    int readID;
    int readEnd;
    int alignID1;
    int alignID2;
    string chromosome1;
    char strand1;
    int position1;
    string chromosome2;
    char strand2;
    int position2;
    int score;
};

std::ostream & operator<<(std::ostream &os, const SplitAlignmentRecord& record)
{
    os << record.libID << "\t";
    os << record.readID << "\t";
    os << record.readEnd << "\t";
    os << record.alignID1 << "\t";
    os << record.alignID2 << "\t";
    os << record.chromosome1 << "\t";
    os << ((record.strand1 == PlusStrand) ? '+' : '-') << "\t";
    os << record.position1 << "\t";
    os << record.chromosome2 << "\t";
    os << ((record.strand2 == PlusStrand) ? '+' : '-') << "\t";
    os << record.position2 << "\t";
    os << record.score << std::endl;
    return os;
}


#endif
