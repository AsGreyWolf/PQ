#include "consensus.h"

BranchOcc* branchOccCreate(Branch* branch, unsigned occurence)
{
    BranchOcc* brc = malloc(sizeof(BranchOcc));
    brc->branch = branch;
    brc->occurence = occurence;
    return brc;
}

void branchOccDelete(BranchOcc* brc)
{
    free(brc);
}


int branchOccCompare(BranchOcc* brc1, BranchOcc* brc2)
{
    return brc1->occurence - brc2->occurence;
}


int vBranchOccCompare(const void* brc1, const void* brc2)
{
    BranchOcc* brcPtr1 = *(BranchOcc**)brc1;
    BranchOcc* brcPtr2 = *(BranchOcc**)brc2;
    return branchOccCompare(brcPtr1, brcPtr2);
}

int branchOccCompareByBranch(BranchOcc* brc1, BranchOcc* brc2)
{
    return branchCompare(brc1->branch, brc2->branch);
}

int vBranchOccCompareByBranch(const void* brc1, const void* brc2)
{
    BranchOcc* brcPtr1 = *(BranchOcc**)brc1;
    BranchOcc* brcPtr2 = *(BranchOcc**)brc2;
    return branchOccCompareByBranch(brcPtr1, brcPtr2);
}


BranchCounter* branchCounterCreate(size_t startMaxSize)
{
    BranchCounter* bc = malloc(sizeof(BranchCounter));
    bc->size = 0;
    bc->maxSize = startMaxSize;
    bc->array = calloc(bc->maxSize, sizeof(BranchOcc*));
    return bc;
}


void branchCounterDelete(BranchCounter* bc)
{
    free(bc->array);
    free(bc);
}

void branchCounterAdd(BranchCounter* bc, Branch* br, size_t branchOccurence)
{
    if (bc->size == bc->maxSize)
    {
        bc->maxSize = bc->maxSize * 3 / 2 + 1;
        bc->array = realloc(bc->array, bc->maxSize * sizeof(BranchOcc*));
    }
    bc->array[bc->size++] = branchOccCreate(br, branchOccurence);
}

void branchCounterInc(BranchCounter* bc, unsigned pos, unsigned diff)
{
    bc->array[pos]->occurence += diff;
}

BranchCounter* branchCount(BranchArray* ba)
{
    BranchCounter* bc = branchCounterCreate(10);
    int i = 0;
    branchCounterAdd(bc, ba->array[i],1);
    for(i = 1; i < ba->size; ++i)
    {
        if (branchCompare(ba->array[i - 1], ba->array[i]))
        {
            branchCounterAdd(bc, ba->array[i], ba->array[i]->weight);
        }
        else
        {
            branchCounterInc(bc, bc->size - 1, ba->array[i]->weight);
        }
    }

    return bc;
}

void branchCounterSort(BranchCounter* brc)
{
    qsort(brc->array, brc->size, sizeof(BranchOcc*), vBranchOccCompare);
    return;
}

void branchCounterSortByBranch(BranchCounter* brc)
{
    qsort(brc->array, brc->size, sizeof(BranchOcc*), vBranchOccCompareByBranch);
}

BranchArray* treeToBranch(Tree* tree, int* permutation, unsigned int weight)
{
    INT p = 1;
    int i = 0;
    int j = 0;
    unsigned branchNum = tree->nodesNum;
    BranchArray* ba;
    BranchArray* result;
    NodeStack* stack;
    Node* curNode = 0;
    Node* nextNode = 0;
    char* isTrivial;

    ba = branchArrayCreate(branchNum);
    result = branchArrayCreate(tree->nodesNum - 1);
    stack = nodeStackCreate(tree->nodesNum);

    if (tree->leavesNum == 0)
    {
        return 0;
    }
    if (tree == 0)
    {
        perror("Function treeToBranch: null Tree pointer");
        exit(1);
/*        raiseError("null Tree pointer", __FILE__, __FUNCTION__, __LINE__); */
    }

    treeWash(tree);

    for(i = 0; i < branchNum; ++i)
    {
        branchArrayAdd(ba, branchCreateWeighted(tree->leavesNum, weight));
    }

    for(i = 0; i < tree->leavesNum; ++i)
    {
        p = 1;
        p = p << (permutation[i] & (intSize - 1));
        ba->array[tree->leaves[i]->pos]->branch[permutation[i] / intSize] |= p;
        tree->leaves[i]->color = BLACK;
    }

    if (tree->nodes[0]->name == 0)
    {
        curNode = tree->nodes[0];
    }
    else
    {
        curNode = tree->nodes[1];
    }

    nodeStackPush(stack, curNode);
    curNode->color = GREY;

    while (stack->curSize != 0)
    {
        curNode = nodeStackPeek(stack);
        nextNode = 0;
        for(i = 0; i < curNode->neiNum; ++i)
        {
            if (curNode->neighbours[i]->color == WHITE)
            {
                nextNode = curNode->neighbours[i];
                break;
            }
        }
        if (nextNode)
        {
            nodeStackPush(stack, nextNode);
            nextNode->color = GREY;
        }
        else
        {
            for(i = 0; i < curNode->neiNum; ++i)
            {
                for(j = 0; j < tree->leavesNum / 64 + 1; ++j)
                {
                    ba->array[curNode->pos]->branch[j] |= \
                            ba->array[curNode->neighbours[i]->pos]->branch[j];
                }
            }
            nodeStackPop(stack);
            curNode->color = BLACK;
        }
    }


    isTrivial = (char*)calloc(sizeof(char), branchNum);
    isTrivial[curNode->pos] = 1;
    for (i = 0; i < tree->leavesNum; ++i)
    {
        isTrivial[tree->leaves[i]->pos] = 1;
    }

    
    for(i = 0; i < ba->size; ++i)
    {
        if (! isTrivial[i]) branchArrayAdd(result, ba->array[i]);
        else branchDelete(ba->array[i]);
    }


    for(i = 0; i < result->size; ++i)
    {
        branchNormalize(result->array[i]);
    }

    free(isTrivial);
    branchArrayDelete(ba);
    nodeStackDelete(stack);
    return result;
}

ParserNode* parserNodeCreate(BranchOcc* branchOcc)
{
    ParserNode* nd = malloc(sizeof(ParserNode));
    nd->branchOcc = branchOcc;
    nd->parent = NULL;
    nd->right = NULL;
    nd->left = NULL;
    nd->treeNode = NULL;
    nd->color = WHITE;
    return nd;
}


void parserNodeDelete(ParserNode* nd)
{
    free(nd);
}

ParserTree* parserTreeCreate()
{
    ParserTree* tree = malloc(sizeof(ParserTree));
    tree->size = 0;
    tree->root = NULL;
    return tree;
}

void parserTreeDelete(ParserTree* tree)
{
    ParserNode** nodes = malloc(sizeof(ParserNode*) * tree->size);
    size_t curStart = 0;
    size_t curEnd = 1;
    size_t newEnd = 0;
    int i = 0;
    nodes[0] = tree->root;

    while(curEnd != tree->size)
    {
        newEnd = curEnd;
        for(i = curStart; i < curEnd; ++i)
        {
            if (nodes[i]->left != NULL)
            {
                nodes[newEnd++] = nodes[i]->left;
            }    
            if (nodes[i]->right != NULL)
            {
                nodes[newEnd++] = nodes[i]->right;
            }
            free(nodes[i]);
        } 
        curStart = curEnd;
        curEnd = newEnd;
    }
    free(nodes);
}


void parserTreePrint(ParserTree* tree)
{
    ParserNode** nodes = malloc(sizeof(ParserNode*) * tree->size * 2);
    size_t curStart = 0;
    size_t curEnd = 1;
    size_t newEnd = 0;
    int i = 0;
    int k = 0;
    nodes[0] = tree->root;

    while(curEnd != curStart)
    {
        newEnd = curEnd;
        printf("Level %d\n", k++);

        for(i = curStart; i < curEnd; ++i)
        {
            if (nodes[i]->left != NULL)
            {
                nodes[newEnd++] = nodes[i]->left;
            }    
            if (nodes[i]->right != NULL)
            {
                nodes[newEnd++] = nodes[i]->right;
            }
            branchPrint(nodes[i]->branchOcc->branch);
        } 
        curStart = curEnd;
        curEnd = newEnd;
        printf("End\n");
    }
    free(nodes);
}

void parserTreeAdd(ParserTree* tree, BranchOcc* branchOcc,
        char* name) // 0 for inner, other for leaf
{
    ParserNode* curNode;
    ParserNode* temp = NULL;
    char isSubset = 1;
    
    ++(tree->size);
    if (tree->root == NULL)
    {
            tree->root = parserNodeCreate(branchOcc);
        if (name) tree->root->treeNode = leafCreate(name);
        return;
    }
    curNode = tree->root;
    temp = NULL;
    isSubset = 1;
    while(true)
    {
        if( (isSubset = 
                        branchIsSubset(branchOcc->branch,
                            curNode->branchOcc->branch)) == -1)
        {
            if (curNode->left != NULL)
            {
                curNode = curNode->left;
            }
            else
            {
                curNode->left = 
                    parserNodeCreate(branchOcc);
                curNode->left->parent = curNode;
                curNode = curNode->left;
                break;
            }
        }
        else if (isSubset == 0)
        {
            if (curNode->right != NULL)
            {
                curNode = curNode->right;
            }
            else
            {
                curNode->right = parserNodeCreate(branchOcc);
                curNode->right->parent = curNode;
                curNode = curNode->right;
                break;
            }
        }
        else // isSubset == 1
        {
            //printf("is subset = 1\n");
            //fprintf(stderr,"Perhaps something've gone wrong. Please,\
            //        make sure you've give branches in reverse order,\
            //        if yes - inform me");
            //exit(1);
            temp = parserNodeCreate(branchOcc);
            if (curNode == tree->root)
            {
                tree->root = temp;
                tree->root->left = curNode;
                curNode->parent = tree->root;
            }
            else
            {
                if (curNode->parent->right == curNode) 
                {
                    curNode->parent->right = temp;
                }
                else
                { 
                    curNode->parent->left = temp;
                }
                temp->parent = curNode->parent;
                curNode->parent = temp;
                temp->left = curNode;
            }
            if ( curNode->right != NULL && (!branchIsSubset(temp->branchOcc->branch, 
                        curNode->right->branchOcc->branch)))
            {
                temp->right = curNode->right;
                curNode->right->parent = temp;
                curNode->right = NULL;
            }

            curNode = temp;
            break;
        }
    }
    if (name) curNode->treeNode = leafCreate(name);

    return;

}

Tree* branchCounterToTree(BranchCounter* bc, char** names)
{
    int i = 0;
    char isSubset = 0;
    ParserNode* temp = NULL;
    ParserNode* curNode = NULL;
    Branch* branch = NULL;
    INT p = 1;
    size_t leavesNum;
    ParserTree* parser;
    Tree* tree;
    size_t curNodesNum = 0;
    size_t curLeavesNum = 0;
    BranchOcc** trivialBranches;
    Branch* orBranch = NULL;
    ParserNode** stack;
    size_t stackSize = 0;
    char ready;

    leavesNum = bc->array[0]->branch->size;
    parser = parserTreeCreate();
    tree = treeCreate(); 
    tree->nodes = (Node**)malloc(sizeof(Node*) * (leavesNum * 2 - 2));
    tree->leaves = (Node**)malloc(sizeof(Node*) * (leavesNum));
    tree->leavesNum = leavesNum;

    branch = branchCreate(leavesNum);
    branch->branch[0] = 1;
    branchNormalize(branch);

    trivialBranches = (BranchOcc**)malloc(sizeof(BranchOcc*) * leavesNum); 
    trivialBranches[0] = branchOccCreate(branch, 0);
    parserTreeAdd(parser, trivialBranches[0], names[0]);

    branchCounterSortByBranch(bc);


    for(i = bc->size - 1; i >= 0; --i)
    {
        parserTreeAdd(parser, bc->array[i], NULL);
    }

    for(i = 1; i < leavesNum; ++i)
    {// add trivial leaves to simplify algorithm
        p = 1;
        p = p << (i & (intSize - 1));
        branch = branchCreate(leavesNum);   
        branch->branch[i / intSize] |= p;
        
        trivialBranches[i] = branchOccCreate(branch, 0);
        parserTreeAdd(parser, trivialBranches[i], names[i]);
    }


    curNode = parser->root;

    stack = malloc(sizeof(ParserNode*) * (leavesNum * 2 - 2));
    stackSize = 0;
    stack[stackSize++] = parser->root;
    ready = 1;
    while(stackSize > 0)
    {
        ready = 1;
        curNode = stack[stackSize - 1];
        if (curNode->left != NULL && curNode->left->color != BLACK)
        {
            stack[stackSize++] = curNode->left;
            ready = 0;
        }
        if (curNode->right != NULL && curNode->right->color != BLACK)
        {
            stack[stackSize++] = curNode->right;
            ready = 0;
        }
        if (ready)
        {
            curNode->color = BLACK;

            if (curNode->left == NULL && curNode->right == NULL)
            {
/*
                      001100
                       /
                 001000 
                       \
                         000100  we are here 
*/
                tree->nodes[curNodesNum++] = curNode->treeNode;
                curNode->treeNode->pos = curNodesNum - 1;
                tree->leaves[curLeavesNum++] = curNode->treeNode;
            }
            else if (curNode->left == NULL )
            {   
/*
                      001100 or 001110
                       /
                 001000 we are here
                       \
                         000100  or 0001100 
*/

                tree->nodes[curNodesNum++] = curNode->treeNode;
                tree->leaves[curLeavesNum++] = curNode->treeNode;
                curNode->treeNode->pos = curNodesNum - 1;
                
                if (curNode->right->right == NULL)
                {
                    tree->nodes[curNodesNum++] = nodeCreate();              
                    tree->nodes[curNodesNum - 1]->pos = curNodesNum - 1;

                    nodeAddNeighbour(tree->nodes[curNodesNum - 1],
                        curNode->right->treeNode, 
                        curNode->right->branchOcc->occurence);
                    nodeAddNeighbour(curNode->right->treeNode,
                        tree->nodes[curNodesNum - 1], 
                        curNode->right->branchOcc->occurence);
                    nodeAddNeighbour(curNode->treeNode,
                        tree->nodes[curNodesNum - 1], 
                        0);
                    nodeAddNeighbour(tree->nodes[curNodesNum - 1],
                        curNode->treeNode, 
                        0);
                    curNode->treeNode = tree->nodes[curNodesNum - 1];
                }
                else // curNode->right->right != NULL
                {
                    if (curNode->treeNode != NULL)
                    {
                        nodeAddNeighbour(curNode->treeNode,
                                curNode->right->treeNode, 0);
                        nodeAddNeighbour(curNode->right->treeNode,
                                curNode->treeNode, 0);
                        curNode->treeNode = curNode->right->treeNode;
                    }
                    else
                    {
                        perror("Something've gone wrong");
                        exit(1);
/*                      raiseError("Something've gone wrong\n", __FILE__, __FUNCTION__, __LINE__); */
                    }
                }
                
            }
            else if (curNode->right == NULL)
            {
                if (curNode == parser->root)
                {
                    nodeAddNeighbour(curNode->left->treeNode, curNode->treeNode,
                        curNode->branchOcc->occurence);
                    nodeAddNeighbour(curNode->treeNode, curNode->left->treeNode,
                        curNode->branchOcc->occurence);

                    tree->nodes[curNodesNum++] = curNode->treeNode;
                    tree->leaves[curLeavesNum++] = curNode->treeNode;
                    curNode->treeNode->pos = curNodesNum - 1;
                }
                else
                {
                    curNode->treeNode = curNode->left->treeNode;
                }
            }
            else //(curNode->right != NULL && curNode->left != NULL)
            {
                /*orBranch = branchOr(curNode->branchOcc->branch,
                        curNode->right->branchOcc->branch);*/
                curNode->treeNode = nodeCreate();
                nodeAddNeighbour(curNode->treeNode, curNode->left->treeNode,
                    curNode->branchOcc->occurence);
                nodeAddNeighbour(curNode->left->treeNode, curNode->treeNode,
                    curNode->branchOcc->occurence);
                nodeAddNeighbour(curNode->treeNode, curNode->right->treeNode,
                    curNode->right->branchOcc->occurence);
                nodeAddNeighbour(curNode->right->treeNode, curNode->treeNode,
                    curNode->right->branchOcc->occurence);
                tree->nodes[curNodesNum++] = curNode->treeNode;
                curNode->treeNode->pos = curNodesNum - 1;
            }
                /*
                else
                {
                    curNode->treeNode = curNode->left->treeNode;
                    nodeAddNeighbour(curNode->treeNode, 
                                     curNode->right->treeNode,
                                     curNode->right->branchOcc->occurence);
                    nodeAddNeighbour(curNode->right->treeNode,
                                     curNode->treeNode,
                                     curNode->right->branchOcc->occurence);
                    
                }
                branchDelete(orBranch);*/
            --stackSize;
        }
    }

    for(i = 0; i < leavesNum; ++i)
    {
        branchDelete(trivialBranches[i]->branch);
        branchOccDelete(trivialBranches[i]);
    }
    
    tree->nodesNum = curNodesNum;

    free(stack);
    parserTreeDelete(parser);
    return tree;
}

// branches sorted by occurence in BranchCounter are required
BranchCounter* majorityRule(BranchCounter* bc, unsigned threshold)
{
    BranchCounter* consensus = branchCounterCreate(10);
    int i = 0;
    int j = 0;
    int k = 0;
    char take = 1;
    for(i = bc->size - 1; i >= 0; --i)
    {
        if (bc->array[i]->occurence > threshold)
        {
            branchCounterAdd(consensus, bc->array[i]->branch, bc->array[i]->occurence);
        }
        else break;
    }
    return consensus;
}

// branches sorted by occurence in BranchCounter are required
BranchCounter* majorityExtendedRule(BranchCounter* bc, unsigned threshold)
{
    BranchCounter* consensus = branchCounterCreate(10);
    int i = 0;
    int j = 0;
    int k = 0;
    char take = 1;
    for(i = bc->size - 1; i >= 0; --i)
    {
        if (bc->array[i]->occurence > threshold)
        {
            branchCounterAdd(consensus, bc->array[i]->branch, bc->array[i]->occurence);
        }
        else break;
    }
    j = i;
    // select non-contradicting branches
    for(i = j; i >= 0; --i)
    {
        take = 1;
        for (k = 0; k < consensus->size; ++k)
        {
            if (branchContradict(consensus->array[k]->branch, bc->array[i]->branch))
            {
                take = 0;
            }   
        }
        if (take)
        {
            branchCounterAdd(consensus, bc->array[i]->branch, bc->array[i]->occurence);
        }
    }
    
    return consensus;
}


Tree* makeConsensus(Tree** treeArray, unsigned int* treesWeight, size_t treeNum, double threshold,
        char extended)  
{
    int i = 0;
    char** initTreeNames = treeGetNames(treeArray[0]);
    int* permutation = getRange(0, treeArray[0]->leavesNum);
    BranchArray* ba = treeToBranch(treeArray[0], permutation, treesWeight[0]);
    BranchArray* temp = NULL;
    char** treeNames = NULL;
    Tree* consensusTree = NULL;
    BranchCounter* consensus = NULL;
    BranchCounter* bc = NULL;
    free(permutation);
    for(i = 1; i < treeNum; ++i)
    {
        treeNames = treeGetNames(treeArray[i]); 
        permutation = calculatePermutation(treeNames, initTreeNames,
                treeArray[i]->leavesNum);
        
        temp = treeToBranch(treeArray[i], permutation, treesWeight[i]);
        
        branchArrayExtend(ba, temp);
        free(permutation);
        branchArrayDelete(temp);
    }
    

    branchArraySort(ba);
    bc = branchCount(ba);
    branchCounterSort(bc);
    if (extended)
    {
        printf("Applying extended majority rule\n");
        consensus = majorityExtendedRule(bc, treeNum * threshold);
    }
    else
    {
        printf("Applying majority rule with threshold %.02lf\n", threshold);
        consensus = majorityRule(bc, treeNum * threshold);
    }

    printf("Chosen branches:\n");
    for (i = 0; i < consensus->size; ++i)
    {
        printf("Branch: %s ", branchToString(consensus->array[i]->branch));
        printf("Occurence: %u\n", consensus->array[i]->occurence);
    }

    consensusTree = branchCounterToTree(consensus, initTreeNames);

    free(initTreeNames);
    free(treeNames);
    for(i = 0; i < consensus->size; ++i) branchOccDelete(consensus->array[i]);
    for(i = 0; i < bc->size; ++i) branchOccDelete(bc->array[i]);
    for(i = 0; i < ba->size; ++i) branchDelete(ba->array[i]);
    branchCounterDelete(consensus);
    branchCounterDelete(bc);
    branchArrayDelete(ba);
    
    return consensusTree;
}
