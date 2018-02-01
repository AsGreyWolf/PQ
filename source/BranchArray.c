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
#include "BranchArray.h"

BranchArray* branchArrayCreate(unsigned startSize)
{
    BranchArray* ba = malloc(sizeof(BranchArray));
    ba->array = calloc(startSize, sizeof(Branch*));
    ba->maxSize = startSize;
    ba->size = 0;
    return ba;
}

void branchArrayDelete(BranchArray* ba)
{
    int i = 0;
    free(ba->array);
    free(ba);
}

void branchArrayAdd(BranchArray* ba, Branch* branch)
{
    if (ba->size && (ba->array[0]->size != branch->size))
    {
        fprintf(stderr, "branchArrayAdd: Branches are not of the same size\n");
        exit(1);

/*        raiseError("Branches are not of the same size", __FILE__, __FUNCTION__, __LINE__);*/
    }

    if (ba->size == ba->maxSize)
    {
        ba->maxSize = ba->maxSize * 3 / 2 + 1;
        ba->array = realloc(ba->array, sizeof(Branch*) * ba->maxSize);
    }
    ba->array[ba->size++] = branch;
}



void branchArrayExtend(BranchArray* dest, BranchArray* source)
{
    unsigned newSize = dest->size + source->size;
    if (newSize >= dest->maxSize)
    {
        dest->maxSize = newSize * 3 / 2 + 1;
        dest->array = realloc(dest->array, sizeof(Branch*) * dest->maxSize);
    }

    memcpy(dest->array + dest->size, source->array, source->size * sizeof(Branch*));
    dest->size += source->size;
}

void branchArraySort(BranchArray* ba)
{
    qsort(ba->array, ba->size, sizeof(Branch*), vBranchCompare);
    return;
}