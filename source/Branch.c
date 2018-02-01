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

#include "Branch.h"

size_t getBitMaskIntSize(unsigned size){
    if (size % intSize == 0){
        return size / intSize;
    }
    return size / intSize + 1;
}


size_t branchGetIntSize(Branch* br){
    return getBitMaskIntSize(br->size);
}

Branch* branchCreate(unsigned size)
{
    Branch* branch;
    branch = (Branch*)malloc(sizeof(Branch));
    branch->size = size;
    branch->branch = (INT*)calloc(sizeof(INT), branchGetIntSize(branch));
    branch->leavesNum = 0;
    return branch;
}


void branchDelete(Branch* branch)
{
    free(branch->branch);
    free(branch);
}

Branch* branchOr(Branch* br1, Branch* br2)
{
    int i = 0;
    Branch* orBranch = NULL;
    if (br1->size != br2->size)
    {
      fprintf(stderr, "branchOR: Branches are not of the same size\n");
      exit(1);
    }
    orBranch = branchCreate(br1->size);
    for(i = 0; i < branchGetIntSize(br1); ++i)
    {
        orBranch->branch[i] = br1->branch[i] | br2->branch[i];
    }
    orBranch->leavesNum = -1;
    return orBranch;
}

void branchOrDest(Branch* br1, Branch* br2, Branch* orBranch){

  int i = 0;
  if (br1->size != br2->size || orBranch->size != br2->size)
  {
    fprintf(stderr, "branchORDest: Branches are not of the same size\n");
    exit(1);
  }
  for(i = 0; i < branchGetIntSize(br1); ++i)
  {
      orBranch->branch[i] = br1->branch[i] | br2->branch[i];
  }
  orBranch->leavesNum = -1;
  return;
}

/*
 for this function we require that there is no intersection between br1 and br2
 so we can set result branch leaves num to the sum of br1 and br2 leaves nums
*/
void branchOrSumDest(Branch* br1, Branch* br2, Branch* orBranch){
  int i = 0;
  if (br1->size != br2->size || orBranch->size != br2->size){
        fprintf(stderr, "branchORSumDest: Branches are not of the same size\n");
        exit(1);
  }
  for(i = 0; i < branchGetIntSize(br1); ++i)
  {
      orBranch->branch[i] = br1->branch[i] | br2->branch[i];
  }
  orBranch->leavesNum = br1->leavesNum + br2->leavesNum;
  return;
}

Branch* branchAnd(Branch* br1, Branch* br2)
{
    int i = 0;
    Branch* andBranch = NULL;
    if (br1->size != br2->size)
    {
        fprintf(stderr, "branchAnd: Branches are not of the same size\n");
        exit(1);
    }
    andBranch = branchCreate(br1->size);

    for(i = 0; i < branchGetIntSize(br1); ++i)
    {
        andBranch->branch[i] = br1->branch[i] & br2->branch[i];
    }
    andBranch->leavesNum = -1;
    return andBranch;
}

void branchAndDest(Branch* br1, Branch* br2, Branch* andBranch){
    int i = 0;
    if (br1->size != br2->size || andBranch->size != br2->size)
    {
        fprintf(stderr, "branchAndDest: Branches are not of the same size\n");
        exit(1);
    }

    for(i = 0; i < branchGetIntSize(br1); ++i)
    {
        andBranch->branch[i] = (br1->branch[i] & br2->branch[i]);
    }
    andBranch->leavesNum = -1;
}

unsigned countZeroRightNum_(INT p)
{

    unsigned int v = p;  // find the number of trailing zeros in v
    int r = 0;           // put the result in r
    static const int Mod37BitPosition[] = // map a bit value mod 37 to its position
    {
      32, 0, 1, 26, 2, 23, 27, 0, 3, 16, 24, 30, 28, 11, 0, 13, 4,
        7, 17, 0, 25, 22, 31, 15, 29, 10, 12, 6, 0, 21, 14, 9, 5,
          20, 8, 19, 18
    };

    //r = Mod37BitPosition[(-v & v) % 37];
    r = Mod37BitPosition[(- ((int)v) & ((int)v) ) % 37];
    return r;
}

unsigned countZeroRightNum(INT p)
{
    if (p == 0)
    {
        return 64;
    }
    uint64_t y = p;
    //find number of trailing zeros
    int r;            // result goes here
    static const char MultiplyDeBruijnBitPosition[64] =
    {
            0, 1, 48, 2, 57, 49, 28, 3, 61, 58, 50, 42, 38, 29, 17, 4,
            62, 55, 59, 36, 53, 51, 43, 22, 45, 39, 33, 30, 24, 18, 12, 5,
            63, 47, 56, 27, 60, 41, 37, 16, 54, 35, 52, 21, 44, 32, 23, 11,
            46, 26, 40, 15, 34, 20, 31, 10, 25, 14, 19, 9, 13, 8, 7, 6
    };
    r = MultiplyDeBruijnBitPosition[((uint64_t)((y & -y) * 0x03F79D71B4CB0A89U)) >> 58];
    //printf("%lu %llu %d\n", p, y, r);
    return r;
}

size_t* branchGetLeavesPos(Branch* br, size_t* leavesNum, size_t maxNum)
{
    unsigned i = 0;
    unsigned j = 0;
    unsigned k = 0;
    size_t* positions;
    unsigned curSize = 0;

    positions = malloc(sizeof(size_t) * br->size);
    for(i = 0; i < branchGetIntSize(br); ++i)
    {
        j = 0;
        while(j < intSize)
        {
            k = countZeroRightNum((br->branch[i]) >> j);
            if (k != intSize)
            {
                positions[curSize++] = k + j +  i * intSize;
            }
            j += k + 1;
            if (curSize > maxNum){
                free(positions);
                *leavesNum = 0;
                return NULL;
            }

        }
    }
    if (!curSize){
        free(positions);
        return 0;
    }
    positions = realloc(positions, sizeof(size_t) * curSize);
    *leavesNum = curSize;
    return positions;
}


void branchNormalize(Branch* br)
{
    int i = 0;
    INT p = 0;
    if ((br->branch[0] & 1) == 1)
    {
        for(i = 0; i < branchGetIntSize(br); ++i)
        {
            br->branch[i] = ~(br->branch[i]);
        }
        --i;
        if (br->size % intSize)
        {
            p = (INT)1 << ((INT)(br->size) % (INT)(intSize) );
            br->branch[i] = br->branch[i] % p;
        }
    }
    return;
}


char* branchToString(Branch* br)
{
    int j = 0;
    int k = 0;
    INT p = 0;
    int curSize = 0;
    char* str = malloc(sizeof(char) * (br->size + 1));
    str[br->size] = '\0';
    for(j = 0; j < branchGetIntSize(br); ++j)
    {
        p = 1;
        curSize = br->size - intSize * j;
        curSize = curSize > intSize ? intSize : curSize;

        for(k = 0; k < curSize; ++k)
        {
            sprintf(str + j * intSize + k, "%lu", (br->branch[j] & p) >> k);
            p <<= 1;
        }
    }
    str[br->size] = '\0';
    return str;
}

void branchPrint(Branch* br)
{
    int j = 0;
    int k = 0;
    INT p = 0;
    int curSize = br->size;

    for(j = 0; j < branchGetIntSize(br); ++j)
    {
        p = 1;
        curSize = curSize > intSize ? intSize : curSize;

        for(k = 0; k < curSize; ++k)
        {
            printf("%lu", (br->branch[j] & p) >> k);
            p <<= 1;
        }
        curSize -= intSize;
    }
    printf("\n");
}

char branchIsSubset(Branch* br1, Branch* br2) // 0 - not 1 - br2 is subset br1,
    //-1 - br1 is subset br2
{
    int i = 0;
    char isSubset = 0;

    if (br1->size != br2->size)
    {
        fprintf(stderr, "branchIsSubset: Branches are not of the same size\n");
        exit(1);

/*        raiseError("Branches are not of the same size", __FILE__, __FUNCTION__, __LINE__); */
    }

    for(i = 0; i < branchGetIntSize(br1); ++i)
    {
        if ((br1->branch[i] & br2->branch[i]) != 0)
        {
            isSubset = 1;
            break;
        }
    }

    if (!isSubset) return 0;

    for (i = 0; i < branchGetIntSize(br2); ++i)
    {
        if ((br1->branch[i] & (~ br2->branch[i])) != 0) return 1;
    }

    return -1;
}


char branchContradict(Branch* br1, Branch* br2)
{
    int i = 0;
    char isEmpty11 = 1;
    char isEmpty00 = 1;
    char isEmpty10 = 1;
    char isEmpty01 = 1;

    if (br1->size != br2->size)
    {
        fprintf(stderr, "branchContradict: Branches are not of the same size\n");
        exit(1);
    }


    for(i = 0; i < branchGetIntSize(br1); ++i)
    {
        if (br1->branch[i] & br2->branch[i]) isEmpty11 = 0;
        if ((~ br1->branch[i]) & (~ br2->branch[i])) isEmpty00 = 0;
        if ((~ br1->branch[i]) & br2->branch[i]) isEmpty01 = 0;
        if (br1->branch[i] & (~ br2->branch[i])) isEmpty10 = 0;
    }


    if (! (isEmpty11 || isEmpty00 || isEmpty01 || isEmpty10)) return 1;
    return 0;
}

Branch* branchReverse(Branch* br)
{
    int i = 0;
    Branch* revBranch = branchCreate(br->size);
    for(i = 0; i < branchGetIntSize(br); ++i)
    {
        revBranch->branch[i] = ~br->branch[i];
    }
    revBranch->leavesNum = -1;
    return revBranch;
}

Branch* branchCopy(Branch* br)
{
    int i = 0;
    Branch* copy = branchCreate(br->size);
    for(i = 0; i < branchGetIntSize(br); ++i)
    {
        copy->branch[i] = br->branch[i];
    }
    copy->leavesNum = br->leavesNum;
    return copy;
}

Branch* branchCopyToDest(Branch* br, Branch* dest){
  int i = 0;
  if (br->size != dest->size){
      fprintf(stderr, "branchCopyToDest: Branches are not of the same size\n");
      exit(1);
  }


  for(i = 0; i < branchGetIntSize(br); ++i){
      dest->branch[i] = br->branch[i];
  }
  dest->leavesNum = br->leavesNum;
  return dest;
}

void branchAddLeafUnsafe(Branch* br, int leafPos){
    INT p = 1;
    p = p << (leafPos & (intSize - 1));
    br->branch[leafPos / intSize] |= p;
    ++br->leavesNum;
}

char branchIsZero(Branch* br)
{
    int i = 0;
    for(i = 0; i < branchGetIntSize(br); ++i)
    {
        if(br->branch[i] != 0) return false;
    }
    return true;
}

int branchCompare(Branch* br1, Branch* br2)
{
    int i = 0;

    if (br1->size != br2->size)
    {
        fprintf(stderr, "branchComplare: Branches are not of the same size\n");
        exit(1);
/*        raiseError("Branches are not of the same size", __FILE__, __FUNCTION__, __LINE__); */
    }
    for(; i < branchGetIntSize(br1); ++i)
    {
        if (br1->branch[i] > br2->branch[i])
        {
            return 1;
        }
        else if (br1->branch[i] < br2->branch[i])
        {
            return -1;
        }
    }
    return 0;
}

int vBranchCompare(const void* branch1, const void* branch2)
{
    int i = 0;
    Branch* brPtr1;
    Branch* brPtr2;

    brPtr1 = *(Branch**)branch1;
    brPtr2 = *(Branch**)branch2;

    return branchCompare(brPtr1, brPtr2);
}
