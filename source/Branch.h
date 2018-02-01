/*  Copyright 2016, 2017 Andrew Sigorskih, Dmitry Penzar, Sergei Spirin

    This file is part of UMAST, taken from our another program PQ.

    UMAST is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    UMAST is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with UMAST (a file named "COPYING.txt").
    If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "add.h"

#ifndef __BRANCH__H
#define __BRANCH__H

typedef struct
{
    unsigned size;
    int leavesNum; // -1 for uncalculated
    INT* branch;
}Branch;

static const char intSize = sizeof(INT) * 8;

size_t getBitMaskIntSize(unsigned size);

size_t branchGetIntSize(Branch* br);

Branch* branchCreate(unsigned size);

void branchDelete(Branch* branch);

unsigned countZeroRightNum_(INT p);

unsigned countZeroRightNum(INT p);

size_t* branchGetLeavesPos(Branch* br, size_t* leavesNum, size_t maxNum);

void branchNormalize(Branch* br);

char* branchToString(Branch* br);

char branchContradict(Branch* br1, Branch* br2);

void branchPrint(Branch* br);

char branchIsSubset(Branch* br1, Branch* br2);

Branch* branchAnd(Branch* br1, Branch* br2);

void branchAndDest(Branch* br1, Branch* br2, Branch* andBranch);

Branch* branchOr(Branch* br1, Branch* br2);

void branchOrSumDest(Branch* br1, Branch* br2, Branch* orBranch);

void branchOrDest(Branch* br1, Branch* br2, Branch* orBranch);

Branch* branchReverse(Branch* br);

Branch* branchCopy(Branch* br);

Branch* branchCopyToDest(Branch* br, Branch* dest);

char branchIsZero(Branch* br);

void branchAddLeafUnsafe(Branch* br, int leafPos);

int branchCompare(Branch* br1, Branch* br2);

int vBranchCompare(const void* branch1, const void* branch2);

#endif
