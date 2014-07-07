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

std::ostream & operator<<(std::ostream &os, const SplitAlignmentRecord& record);

std::istream & operator>>(std::istream &is, SplitAlignmentRecord& record);


#endif
