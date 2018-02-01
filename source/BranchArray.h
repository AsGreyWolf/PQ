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

#ifndef __BRANCH_ARRAY_H
#define __BRANCH_ARRAY_H

#include "Branch.h"

typedef struct
{
    Branch** array;
    unsigned size;
    unsigned maxSize;
}BranchArray;

BranchArray* branchArrayCreate(unsigned startSize);
void branchArrayDelete(BranchArray* ba);
void branchArrayAdd(BranchArray* ba, Branch* branch);
void branchArraySort(BranchArray* ba);
void branchArrayExtend(BranchArray* dest, BranchArray* source);

#endif
