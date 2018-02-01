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

#include "BranchAllocator.h"
#include <assert.h>


size_t getBranchFullSize(unsigned treeLeavesNum){
    return sizeof(Branch) + sizeof(INT) * getBitMaskIntSize(treeLeavesNum);
}

BranchAllocator* branchAllocatorCreate(unsigned capacity, unsigned treeLeavesNum){
    size_t branchFullSize = getBranchFullSize(treeLeavesNum);
    BranchAllocator* allocator = malloc(sizeof(BranchAllocator));
    allocator->length = 0;
    allocator->capacity = capacity;
    allocator->treeLeavesNum = treeLeavesNum;
    allocator->store = calloc(capacity, branchFullSize);
    allocator->branchFullSize = branchFullSize;
    return allocator;
}

void branchAllocatorDelete(BranchAllocator* allocator){
    free(allocator->store);
    free(allocator);
}

Branch* branchAllocatorGetBranch(BranchAllocator* allocator){
    if (allocator->length == allocator->capacity){
        perror("Memory limit for branch allocator");
        exit(1);
    }

    Branch* br = (Branch*)(allocator->store + allocator->length * (allocator->branchFullSize));
    ++allocator->length;
    br->branch = (INT*)( (char*)br + sizeof(Branch));
    br->leavesNum = 0;
    br->size = allocator->treeLeavesNum;
    return br;
}
