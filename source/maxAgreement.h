#ifndef _MAX_AGREEMENT_H_
#define _MAX_AGREEMENT_H_
#include "Tree.h"
#include "add.h"
#include "Branch.h"
#include "BranchArray.h"
#include "BranchAllocator.h"
#include <stdlib.h>
#include <string.h>

//const int nodeTypeCount = 6;

BranchArray* treeRootedToBranchArray(Tree* tree, int* permutation);
void stackCycle(Node* node, NodeStack* stack);
Tree* deleteLeaves(Tree* tree, char** leavesToDelete, int leavesToDeleteAmount);
Tree** treesPrune(Tree* tree1, Tree* tree2);
Tree* treeRoot(Tree* tree, unsigned nodeID, unsigned neighbourID, char newTree);
Tree* treeUnRoot(Tree* tree, char newTree);
unsigned* treeRootAndTopSort(Tree* tree, unsigned nodeID, unsigned neighbourID, unsigned* setPermutation);
int* calculateLeavesPermutation(Tree* tree1, Tree* tree2);
void branchCalculateLeavesPosNum(Branch* br);
int branchGetLeavesPosNum(Branch* br);
int* getTreeLeavesPos(Tree* tree);
unsigned* branchToLeavesArr(Branch* br, unsigned leavesNum);
int** getAllRoots(Tree* tree);
int* getRandMaxBranch(Branch*** TAB, int rows, int cols);
Tree* makeMAST(Branch* br, Tree* tree1);
int find3MaxPos(int a0, int a1, int a2);
int findParent(int elementID, Tree* tree, int* setPermutation);
int findChilds(int elementID, Tree* tree, int* setPermutation, int* child1,
    int* child2);
void countVariants(Branch*** TAB, int a, int w, int b, int c, int x, int y,
    int* bestValue, int* bestPos);    
Branch* MAST(Tree* tree1, Tree* tree2, unsigned* set1, unsigned* set2, unsigned* setPermutation1, unsigned* setPermutation2);
Branch* UMASTStep(Tree* intree1, Tree* intree2, unsigned* set1, unsigned* set2,
    unsigned* setPermutation1, unsigned* setPermutation2, Branch**** TAB, Branch** rootRow,
    BranchAllocator* brAllocator);
Tree* makeUMASTTree(Branch* br, Tree* tree1);
Tree* UMAST(Tree* tree1, Tree* tree2);
#endif