#include "nniTree.h"
#define GRADUS 12000000

TreeWithScore* simpleNNI(Tree* inTree, HashAlignment* alignment, PWM* pwmMatrix, \
        int alpha, GapOpt gapOpt, INT**** hashScore)
{   
    int i, j, variant;
    int isUpdate;
    char** treeNames;
    char** seqNames;
    int* permutation;
    Tree* resultTree;
    INT score;
    TreeWithScore* result;
    TreeWithScore* curBestTree;
    int counter;

    treeNames = treeGetNames(inTree);
    seqNames = hashAlignmentGetSeqNames(alignment);
    permutation = calculatePermutation(treeNames, seqNames,
                                       alignment->alignmentSize);
    free(treeNames);
    free(seqNames);

    resultTree = treeCopy(inTree, 1);

    score = countScoreHash(alignment, resultTree, pwmMatrix,
                           alpha, gapOpt, hashScore, permutation);
    printf("Performing simple NNI hill climbing, initial score is %lu\n", score);
    result = treeWithScoreCreate(resultTree, score);

    isUpdate = TRUE;
    while (isUpdate)
    {
        isUpdate = FALSE;
        curBestTree = treeWithScoreCreate(treeCopy(result->tree, 0), result->score);
        //fprintf(stdout, "Next round\n");
        counter = 0;
        for(i = 0; (i < result->tree->nodesNum) && (!isUpdate); ++i)
        {
            if (result->tree->nodes[i]->neiNum == 3)
            {
                for(j = 0; (j < 3)  && (!isUpdate); ++j)
                {
                    if (result->tree->nodes[i]->neighbours[j]->pos > i && \
                        result->tree->nodes[i]->neighbours[j]->neiNum == 3)
                    {
                        for (variant = 1; (variant < 3)  && (!isUpdate); ++variant)
                        {
                            /*fprintf(stdout, "Watch branch %d %d %d\n", 
                             *                 i,  result->tree->nodes[i]->neighbours[j]->pos, variant); */
                            result->tree = treeNNIMove(result->tree, i, j, variant, 0, 1);
                            score = countScoreHash(alignment, result->tree, pwmMatrix, \
                                                   alpha, gapOpt, hashScore, permutation);
                            counter++;
                            printf("%d NNIs regarded\r", counter);
                            if (score > curBestTree->score)
                            {
                                isUpdate = TRUE;
                                printf("\nTree with score %lu found\n", score);
                                treeDelete(curBestTree->tree);
                                curBestTree->tree = treeCopy(result->tree, 0);
                                curBestTree->score = score;
                            }
                            result->tree = treeNNIMove(result->tree, i, j, variant, 0, 0); /* backward */
                        } /* for variant */
                    } /* if inner branch */
                } /* for adjacent branches */
            } /* if node */
        } /* for vertices */
        treeWithScoreDelete(result);
        result = curBestTree;
    } /* while new score received */
    printf("\nNo better trees found\n");
    free(permutation);
    return result;
} /* simpleNNI */

TreeWithScore* gradientNNI(Tree* inTree, HashAlignment* alignment, PWM* pwmMatrix, \
        int alpha, GapOpt gapOpt, INT**** hashScore)
{   
    int i, j, variant;
    int isUpdate;
    char** treeNames;
    char** seqNames;
    int* permutation;
    Tree* resultTree;
    INT score;
    TreeWithScore* result;
    TreeWithScore* curBestTree;

    treeNames = treeGetNames(inTree);
    seqNames = hashAlignmentGetSeqNames(alignment);
    permutation = calculatePermutation(treeNames, seqNames,
                                       alignment->alignmentSize);
    free(treeNames);
    free(seqNames);

    resultTree = treeCopy(inTree, 1);

    score = countScoreHash(alignment, resultTree, pwmMatrix,
                           alpha, gapOpt, hashScore, permutation);
    printf("Performing gradient NNI hill climbing, initial score is %lu\n", score);
    result = treeWithScoreCreate(resultTree, score);

    isUpdate = TRUE;
    while (isUpdate)
    {
        isUpdate = FALSE;
        curBestTree = treeWithScoreCreate(treeCopy(result->tree, 0), result->score);
        //fprintf(stdout, "Next round\n");
        for(i = 0; i < result->tree->nodesNum; ++i)
        {
            if (result->tree->nodes[i]->neiNum != 1)
            {
                for(j = 0; j < 3; ++j)
                {
                    if (result->tree->nodes[i]->neighbours[j]->pos > i && \
                        result->tree->nodes[i]->neighbours[j]->neiNum != 1)
                    {
                        for (variant = 1; variant < 3; ++variant)
                        {
                            /*fprintf(stdout, "Watch branch %d %d %d\n", 
                             *                 i,  result->tree->nodes[i]->neighbours[j]->pos, variant); */
                            result->tree = treeNNIMove(result->tree, i, j, variant, 0, 1);
                            score = countScoreHash(alignment, result->tree, pwmMatrix,
                                                   alpha, gapOpt, hashScore, permutation);
                            if (score > curBestTree->score)
                            {
                                isUpdate = TRUE;
                                printf("Tree with score %lu found\n", score);
                                treeDelete(curBestTree->tree);
                                curBestTree->tree = treeCopy(result->tree, 0);
                                curBestTree->score = score;
                            }
                            result->tree = treeNNIMove(result->tree, i, j, variant, 0, 0); /* backward */
                        } /* for variant */
                    } /* if inner branch */
                } /* for adjacent branches */
            } /* if node */
        } /* for vertices */
        treeWithScoreDelete(result);
        result = curBestTree;
    } /* while new score received */
    printf("No better trees found\n");
    free(permutation);
    return result;
} /* gradientNNI */

TrajectoryElement* trajectoryElementCreate(TreeWithScore* treeWs, unsigned long int time)
{
    TrajectoryElement* trajectoryElement;
    Tree* newTree;
    TreeWithScore* newTreeWS;

    newTree = treeCopy(treeWs->tree, 0);
    trajectoryElement = (TrajectoryElement*)malloc(sizeof(TrajectoryElement));
    trajectoryElement->treeWS = treeWithScoreCreate(newTree, treeWs->score);
    trajectoryElement->time = time;
    trajectoryElement->next = NULL;
    return trajectoryElement;
}  /* trajectoryElementCreate */

void trajectoryElementDelete(TrajectoryElement* elem)
{
    treeWithScoreDelete(elem->treeWS);
    free(elem);
    return;
}  /* trajectoryElementDelete */
 
Trajectory* trajectoryCreate(unsigned long int trajectoryTime, unsigned int temperature)
{
    Trajectory* trajectory;

    trajectory = (Trajectory*)malloc(sizeof(Trajectory));
    trajectory->bestPoint = NULL;
    trajectory->head = NULL;
    trajectory->tail = NULL;
    trajectory->trajectoryTime = trajectoryTime;
    trajectory->temperature = temperature;
    trajectory->size = 0;
    return trajectory;
}  /* trajectoryCreate */

void trajectoryDelete(Trajectory* trajectory)
{
    TrajectoryElement* curElement;
    TrajectoryElement* nextElement;

    curElement = trajectory->head;
    while (curElement != NULL)
    {
        nextElement = curElement->next;
        trajectoryElementDelete(curElement);
        curElement = nextElement;
    }
    free(trajectory);
    return;
} /* trajectoryDelete */

void trajectoryAdd(Trajectory* trajectory, TreeWithScore* treeWS, unsigned int time)
{
    TrajectoryElement* newElement;
    TrajectoryElement* tmpElement;

    newElement = trajectoryElementCreate(treeWS, time);
    if (trajectory->size == 0)
    {
        trajectory->head = newElement;
        trajectory->bestPoint = newElement;
    }
    else
    {
        trajectory->tail->next = newElement;
        if (treeWS->score > trajectory->bestPoint->treeWS->score)
        { 
            /* printf("New best point: old %lu, new %lu\n", trajectory->bestPoint->treeWS->score, treeWS->score); */
            trajectory->bestPoint = newElement;
        }
    }

    trajectory->tail = newElement;
    trajectory->tail->next = NULL;
    ++(trajectory->size);
    return;
} /* trajectoryAdd */


char isNextPoint(int prevScore, int newScore, unsigned long revTemp)
// 0 if not take to traectory, 1 if take;
{
    double minProbability;
    double coeff;
    double probability;

    if (newScore >= prevScore)
    {
        return 1;
    }
    else
    {
        minProbability = (double) rand();
        minProbability /= RAND_MAX; 
        coeff = ((double)(newScore - prevScore))/prevScore;
        coeff *= revTemp;
        probability = exp(coeff);
        if ( probability > minProbability)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}  /* isNextPoint */

char isNextPoint2(int prevScore, int newScore, unsigned long revPrevTemp, unsigned long revNewTemp)
// 0 if not take to traectory, 1 if take;
{
    double minProbability;
    double coeff;
    double probability;
    int maxScore = prevScore < newScore ? newScore : prevScore;

        minProbability = (double) rand();
        minProbability /= RAND_MAX; 
        coeff = ((double)(newScore - prevScore))/maxScore;
        coeff = coeff*revPrevTemp - coeff*revNewTemp;
        probability = exp(coeff);
        if ( probability > minProbability)
        {
            return 1;
        }
        else
        {
            return 0;
        }
}  /* isNextPoint */

int trajectoryNNIStep(
		unsigned int mcStyle, 
		unsigned int newMcStyle, 
		HashAlignment* alignment,
                PWM* pwmMatrix,
		int alpha,
	       	GapOpt gapOpt,
	       	INT**** hashScore,
		int* permutation,
		TreeWithScore* curPoint, 
		double* curTemperature, double startTemperature, double stepTemperature, 
		unsigned long int* curTime, unsigned long int trajectoryTime, 
		double*** p,
		Trajectory* resultTrajectory) {
    int trLen = 0;
    if ( (mcStyle == 1) || (newMcStyle == 1) || (mcStyle == 4))
    { 
        for(int i = 0; (i < curPoint->tree->nodesNum) && (*curTime < trajectoryTime); ++i)
        {
            if (curPoint->tree->nodes[i]->neiNum != 1) /* if vertice i is a node */
            {
                for(int j = 0; (j < 3) && (*curTime < trajectoryTime); ++j)
                {
                    if (curPoint->tree->nodes[i]->neighbours[j]->pos > i && /* to regard each branch only once */
                        curPoint->tree->nodes[i]->neighbours[j]->neiNum != 1) /* branch is not trivial */
                    {
                        for(int variant = 1; (variant < 3) && (*curTime < trajectoryTime); ++variant)
                        {
                            curPoint->tree = treeNNIMove(curPoint->tree, i, j, variant, 0, 1);
                            INT score = countScoreHash(alignment, curPoint->tree, pwmMatrix,
                                                   alpha, gapOpt, hashScore, permutation);
                            //printf("Second: %lu Temperature : %lf, %lu, %lu\n",  *curTime, *curTemperature, curPoint->score, score);

                            if (isNextPoint(curPoint->score, score, (unsigned long)(GRADUS/ *curTemperature)))
                            { 
                                trLen++;
                                trajectoryAdd(resultTrajectory,
                                              treeWithScoreCreate(treeCopy(curPoint->tree,0), score),
                                              *curTime); 
				if ((int)*curTemperature!=0)
                                	printf("%lu\t%d\t%lu\t%lu\n", *curTime, (int)*curTemperature, score,
                                        	resultTrajectory->bestPoint->treeWS->score);
                                curPoint->score = score;
                            }
                            else
                            {
                                // return to prev state
                                curPoint->tree = treeNNIMove(curPoint->tree, i, j, variant, 0, 0);
                            }
                            ++*curTime;
                            *curTemperature -= stepTemperature;
                        } /* for variants of NNI */
                    } /* if branch is inner and not yet regarded */
                } /* for adjacent branches */
            } /* if vertice i is a node */
        }  /* for all vertices */
    } /* while maximum number of steps is not reached, mcStyle 1*/

    if (mcStyle == 2)
    { 
        int new = FALSE;
        int* nodePermutation = getPermutation(curPoint->tree->nodesNum);
        for(int ii = 0; (ii < curPoint->tree->nodesNum) && (!new) && (*curTime < trajectoryTime); ++ii)
        { 
            int i = nodePermutation[ii];        
            if (curPoint->tree->nodes[i]->neiNum != 1) /* if vertice i is a node */
            {
                for(int j = 0; (j < 3) && (!new) && (*curTime < trajectoryTime); ++j)
                {
                    if (curPoint->tree->nodes[i]->neighbours[j]->pos > i && /* to regard each branch only once */
                        curPoint->tree->nodes[i]->neighbours[j]->neiNum != 1) /* branch is not trivial */
                    {
                        for(int variant = 1; variant < 3 && (!new) && (*curTime < trajectoryTime); ++variant)
                        {
                            /* printf("Second: %lu Temperature : %u\n",  *curTime, *curTemperature); */
                            curPoint->tree = treeNNIMove(curPoint->tree, i, j, variant, 0, 1);
                            INT score = countScoreHash(alignment, curPoint->tree, pwmMatrix,
                                                   alpha, gapOpt, hashScore, permutation);

                            if (isNextPoint(curPoint->score, score, (unsigned long)(GRADUS/ *curTemperature)))
                            { 
                                trLen++;
                                new = TRUE;
                                trajectoryAdd(resultTrajectory,
                                              treeWithScoreCreate(treeCopy(curPoint->tree,0), score),
                                              *curTime); 
                                printf("%lu\t%d\t%lu\t%lu\n", *curTime, (int)*curTemperature, score,
                                        resultTrajectory->bestPoint->treeWS->score);
                                curPoint->score = score;
                            }
                            else
                            {
                                // return to prev state
                                curPoint->tree = treeNNIMove(curPoint->tree, i, j, variant, 0, 0);
                            }
                            ++*curTime;
                            *curTemperature -= stepTemperature;
                        } /* for variants of NNI */
                    } /* if branch is inner and not yet regarded */
                } /* for adjacent branches */
            } /* if vertice i is a node */
        }  /* for all vertices */
        free(nodePermutation);
    } /* while maximum number of steps is not reached, mcStyle 2 */

    if (mcStyle == 3) 
    { 

            double Z = 1.0; /* Statistical sum; 1 is for the current state */
            for(int i = 0; i < curPoint->tree->nodesNum; ++i)
            { 
                if (curPoint->tree->nodes[i]->neiNum != 1) /* if vertice i is a node */
                {
                    for(int j = 0; j < 3; ++j)
                    {
                        if (curPoint->tree->nodes[i]->neighbours[j]->pos > i && /* to regard each branch only once */
                            curPoint->tree->nodes[i]->neighbours[j]->neiNum != 1) /* branch is not trivial */
                        {
                            for(int variant = 1; variant < 3; ++variant)
                            {
                                curPoint->tree = treeNNIMove(curPoint->tree, i, j, variant, 0, 1);
                                INT score = countScoreHash(alignment, curPoint->tree, pwmMatrix,
                                                       alpha, gapOpt, hashScore, permutation);
                                double coeff = ((double)score - curPoint->score)/curPoint->score;
                                coeff /= *curTemperature;
                                coeff *= GRADUS;
                                if (coeff > 10) coeff = 10.0;
                                p[i][j][variant - 1] = exp(coeff);
                                Z += p[i][j][variant - 1];
                                ++*curTime;
                                curPoint->tree = treeNNIMove(curPoint->tree, i, j, variant, 0, 0); /* return to previous */
                            } /* for variants of NNI */
                        } /* if branch is inner and not yet regarded */
                    } /* for adjacent branches */
                } /* if vertice i is a node */
            }  /* for all vertices */

            /**** choose i, j, variant to move ****/
            int choice = rand();
            double partialStatSum = 0.0;
            int randVariant = 0; /* means current state */
	    int randI, randJ;
            for(int i = 0; (i < curPoint->tree->nodesNum) && (randVariant == 0); ++i)
            { 
                if (curPoint->tree->nodes[i]->neiNum != 1) /* if vertice i is a node */
                {
                    for(int j = 0; (j < 3) && (randVariant == 0); ++j)
                    {
                        if (curPoint->tree->nodes[i]->neighbours[j]->pos > i && /* to regard each branch only once */
                            curPoint->tree->nodes[i]->neighbours[j]->neiNum != 1) /* branch is not trivial */
                        {
                            for(int variant = 1; (variant < 3) && (randVariant == 0); ++variant)
                            {
                                partialStatSum += p[i][j][variant - 1];
                                if (partialStatSum*RAND_MAX/Z > choice) 
                                {
                                    randVariant = variant; /* randVariant > 0 now! */
                                    randI = i;
                                    randJ = j;
                                } /* if */
                            } /* for variant */
                        } /* if */
                    } /* for j */
                } /* if */
            } /* for i */ 
            if (randVariant > 0) /* a new tree is chosen */
            {
                curPoint->tree = treeNNIMove(curPoint->tree, randI, randJ, randVariant, 0, 1);
                INT score = countScoreHash(alignment, curPoint->tree, pwmMatrix,
                                       alpha, gapOpt, hashScore, permutation);
                trLen++;
                Tree* newTree = treeCopy(curPoint->tree, 0);
                TreeWithScore* newTreeWS = treeWithScoreCreate(newTree, score);
                trajectoryAdd(resultTrajectory, newTreeWS, *curTime);
                treeWithScoreDelete(newTreeWS);
                printf("%lu\t%d\t%lu\t%lu\n", *curTime, (int)*curTemperature, score,
                        resultTrajectory->bestPoint->treeWS->score);
                curPoint->score = score;
            }
            *curTemperature = startTemperature - stepTemperature * *curTime;
    } /* if mcStyle is 3 */
    return trLen;
}

Trajectory* trajectoryNNI(Tree* inTree, HashAlignment* alignment,
                          PWM* pwmMatrix,  int alpha, GapOpt gapOpt, INT**** hashScore,
                          unsigned long trajectoryTime, unsigned int temperature, unsigned int mc3chains, unsigned int mcStyle)
    // traectoryTime means how many iterations should algorithm do 
{
    int treeNum = mc3chains;    
    char** treeNames;
    char** seqNames;
    int* permutation;
    int* trLen = malloc(sizeof(int) * treeNum);
    double* curTemperature = malloc(sizeof(double) * treeNum);
    Tree* curTree;
    INT score;
    TreeWithScore** curPoint = malloc(sizeof(TreeWithScore*) * treeNum);
    Trajectory** resultTrajectory = malloc(sizeof(Trajectory*) * treeNum);
    unsigned long int* curTime = malloc(sizeof(unsigned long int) * treeNum);
    double step;
    unsigned int newMcStyle = 0;
    double*** p = NULL;

    srand(time(NULL)); // need for correct run isNextPoint;

    treeNames = treeGetNames(inTree);
    seqNames = hashAlignmentGetSeqNames(alignment);
    permutation = calculatePermutation(treeNames, seqNames, alignment->alignmentSize);
    curTree = treeCopy(inTree, 1);
    treeWash(inTree);
    score = countScoreHash(alignment, curTree, pwmMatrix, 
                           alpha, gapOpt, hashScore, permutation);
    for (int i = 0; i < treeNum; i++) {
	trLen[i] = 0;
	curTime[i] = 0;
        curPoint[i] = treeWithScoreCreate(treeCopy(curTree, 1), score);
    	resultTrajectory[i] = trajectoryCreate(trajectoryTime, temperature);
    	trajectoryAdd(resultTrajectory[i], curPoint[i], 0);
    }
    step = (double)temperature / trajectoryTime;

    if (!((mcStyle == 1) || (mcStyle == 2) || (mcStyle == 3) || (mcStyle == 4))) 
    {
        fprintf(stderr, "Warning: unknown mcStyle %u, set to default (1)\n", mcStyle);
        newMcStyle = 1;
    }

    printf("Peforming Monte-Carlo NNI search\n");
    printf("Initial temperature is %u\n%lu steps, style %u\n", temperature, trajectoryTime, mcStyle); 
    printf("Initial score is %lu\n", score);
    printf("Step\tTemp.\tScore\tBest score\n");
    if (mcStyle == 3) {
        p = (double***) malloc(sizeof(double**)*(curPoint[0]->tree->nodesNum));
        for(int i = 0; i < curPoint[0]->tree->nodesNum; ++i)
        {
            p[i] = (double**) malloc(sizeof(double*)*3);
            for(int j = 0; j < 3; ++j)
            {
                p[i][j] = (double*) malloc(sizeof(double)*2);
                for(int variant = 1; variant < 3; ++variant)
                {
                    p[i][j][variant - 1] = 0.0;
                }
            }
        }
    }
    for (int i = 0; i < treeNum; i++) {
        curTemperature[i] = (double) temperature;
    }
    if (mcStyle == 4) {
        for (int i = 0; i < treeNum; i++) {
            curTemperature[i] = temperature * 1.0 / (treeNum-1) * i;
	    if (curTemperature[i] == 0) curTemperature[i]=1;
	}
	step = 0;
    }
    while(curTime[0] < trajectoryTime) {
#pragma omp parallel for
        for (int i = 0; i < treeNum; i++) {
        	trLen[i] += trajectoryNNIStep(mcStyle, newMcStyle, 
			alignment,
			pwmMatrix,
			alpha,
			gapOpt,
			hashScore,
			permutation,
			curPoint[i], 
			&curTemperature[i], temperature, step, 
			&curTime[i], trajectoryTime, 
			p,
			resultTrajectory[i]);
	}
	if (mcStyle == 4) /*for(int i=0;i<20;i++)*/{
	    int leftTree = 0;
	    int rightTree = 0;
	    while (leftTree == rightTree) {
	    	leftTree = rand() % treeNum;
	    	rightTree = rand() % treeNum;
	    }
            if (isNextPoint2(
				    curPoint[leftTree]->score,
				    curPoint[rightTree]->score,
				    (unsigned long)(GRADUS/ curTemperature[leftTree]),
				    (unsigned long)(GRADUS/ curTemperature[rightTree]))) {
		TreeWithScore* tmp = curPoint[leftTree];
		curPoint[leftTree] = curPoint[rightTree];
		curPoint[rightTree] = tmp;
		printf("Swapping %d and %d\n", leftTree, rightTree);
                trajectoryAdd(resultTrajectory[leftTree],
                    treeWithScoreCreate(treeCopy(curPoint[leftTree]->tree,0), curPoint[leftTree]->score),
                    curTime[leftTree]); 
                trajectoryAdd(resultTrajectory[rightTree],
                    treeWithScoreCreate(treeCopy(curPoint[rightTree]->tree,0), curPoint[rightTree]->score),
                    curTime[rightTree]); 
                printf("%lu\t%d\t%lu\t%lu\n", curTime[leftTree], (int)curTemperature[leftTree], curPoint[leftTree]->score,
                                        resultTrajectory[leftTree]->bestPoint->treeWS->score);
                printf("%lu\t%d\t%lu\t%lu\n", curTime[rightTree], (int)curTemperature[rightTree], curPoint[rightTree]->score,
                                        resultTrajectory[rightTree]->bestPoint->treeWS->score);
	    }
	}
    } /* while maximum number of steps is not reached */
    printf("Trajectory length is %u\n", trLen[0]);
    return resultTrajectory[0];
}  /* trajectoryNNI */
