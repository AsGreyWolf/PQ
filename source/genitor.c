#include "genitor.h"
//#define MAX_GA_ITER 18
//#define MAX_LEADER_ITER 10001

unsigned treeIsUnique(TreeWithScore* intree, int pos, TreeWithScore** trees, unsigned treeNum)
{
	int i;
	unsigned unique = TRUE;
	Tree* intersection;

	for (i = 0; i < treeNum; i++)
	{
		if (i != pos)
		{
			if (trees[i]->score == intree->score)
			{
				intersection = UMAST(trees[i]->tree, intree->tree);

				if (intersection->leavesNum == intree->tree->leavesNum)
				{
					unique = FALSE;
				}

				treeDelete(intersection);
			}
		}
	}
	return unique;
}

TreeWithScore* crossover(TreeWithScore* tree1, TreeWithScore* tree2,\
			HashAlignment* alignment, int alpha,\
			 GapOpt gapOpt, PWM* pwmMatrix, INT**** hashScore)
{
	int* indexes;
	int* permutation;
	unsigned i, j;
	TreeWithScore* result;
	char** seqNames;
	char** treeNames;
	char** leavesToGrow;

	indexes = (int*)calloc(sizeof(int), tree1->tree->leavesNum);

	for (i = 0; i < tree1->tree->leavesNum; i++)
	{
		indexes[i] = 0;
	}

	result = treeWithScoreCreate(UMAST(tree1->tree, tree2->tree), 0);
	leavesToGrow = (char**)calloc(sizeof(char*), tree1->tree->leavesNum);

	for (i = 0; i < tree1->tree->leavesNum; i++)
	{
		for (j = 0; j < result->tree->leavesNum; j++)
		{
			if (strcmp(tree1->tree->leaves[i]->name, result->tree->leaves[j]->name) == 0)
			{
				indexes[i] = 1;
				break;
			}
		}
	}

	for (i = 0; i < result->tree->leavesNum; i++)
	{
		leavesToGrow[i] = result->tree->leaves[i]->name;
	}

	i = result->tree->leavesNum;

	for(j = 0; j < tree1->tree->leavesNum; j++)
	{
		if (i > tree1->tree->leavesNum) 
		{
			fprintf(stderr, "failed to copy leaves, genitor:crossover\n");
			exit(1);
		}
		if (indexes[j] == 0) 
		{
			leavesToGrow[i] = tree1->tree->leaves[j]->name;
			i++;
		}
	}

	seqNames = hashAlignmentGetSeqNames(alignment);
	permutation = calculatePermutation(leavesToGrow, seqNames, alignment->alignmentSize); //leavesNum == alignmensSise?

	for (i = result->tree->leavesNum; i < tree1->tree->leavesNum; i++)
	{
		result = getBestChild(alignment, result, leavesToGrow[i], alpha, gapOpt,\
        					pwmMatrix, hashScore, permutation);
	}

	treeLCAFinderCalculate(result->tree);
	result->score = 0; //maybe not needed
	result->score = countScoreHash(alignment, result->tree, pwmMatrix, alpha, gapOpt,\
								hashScore, permutation);
	free(indexes); 
	free(permutation); 
	free(seqNames); 
	free(treeNames); 
	free(leavesToGrow);
	return result;
}

TreeWithScore* genitor(TreeWithScore** trees, unsigned treeNum, HashAlignment* alignment,\
					int alpha, GapOpt gapOpt, PWM* pwmMatrix, INT**** hashScore, unsigned iterNum, \
					unsigned iterNew, unsigned iterLim)
{
	INT leaderScore;
	int* permutation;
	char** seqNames;
	char** treeNames;
	unsigned i, j, k;
	unsigned t = 0;
	unsigned failure = 0;
	unsigned leaderIter = 0;
	unsigned maxFailureStreak = 0;
	unsigned tmpFailureStreak = 0;
	unsigned maxLeaderStreak = 0;
	unsigned* uniqueSet;
	TreeWithScore* offspring;
	TreeWithScore** initPop;
	TreeWithScore** population;

	srand(time(NULL));
	seqNames = hashAlignmentGetSeqNames(alignment);
	uniqueSet = (unsigned*)calloc(sizeof(unsigned), treeNum);
	offspring = (TreeWithScore*)malloc(sizeof(TreeWithScore));
	initPop = (TreeWithScore**)malloc(sizeof(TreeWithScore*) * treeNum);
	population = (TreeWithScore**)malloc(sizeof(TreeWithScore*) * treeNum);

	for (i = 0;  i < treeNum; i++)
	{
		initPop[i] = treeWithScoreCreate(treeCopy(trees[i]->tree, 1), trees[i]->score);
	}

	k = 0;
	for (i = 0; i < treeNum; i++)
	{
		uniqueSet[i] = treeIsUnique(initPop[i], i, initPop, treeNum);
		k += uniqueSet[i];
	}
	printf("%d unique trees\n", k);

	if (k == 0)
	{
		for (i = 1; i < treeNum; i++)
		{
			if (initPop[i]->score != initPop[0]->score)
			{
				k++;
				break;
			}
		}
		if (k == 0)
		{
			offspring = treeWithScoreCreate(treeCopy(initPop[0]->tree, 1), initPop[0]->score);
			printf("Whole population consists of a single tree, no need for genitor run;\n");
			printf("Returning the only tree..\n");
			free(uniqueSet);
			free(population);
			for (i = 0; i < treeNum; i++)
			{
				treeWithScoreDelete(initPop[i]);
			}
			free(initPop);
			return offspring;
		}
	}

	population[0] = treeWithScoreCreate(treeCopy(initPop[0]->tree, 1), initPop[0]->score);
	k = 1;

	for (i = 1; i < treeNum; i++)
	{
		if (treeIsUnique(initPop[i], -1, population, k))
		{
			population[k] = treeWithScoreCreate(treeCopy(initPop[i]->tree, 1), initPop[i]->score);
			k++;
		}
	}

	t = 0;
	failure = 0;
	treeWithScoreSort(population, k); 
	leaderScore = population[k - 1]->score;
	
	while ((t < iterNum) && (k < treeNum))
	{
		i = rand() % k;
		j = rand() % k;
		while (i == j) 
		{
			j = rand() % k;
		}

		offspring = crossover(population[i], population[j], alignment, alpha,\
							gapOpt, pwmMatrix, hashScore);

		treeNames = treeGetNames(offspring->tree);
        permutation = calculatePermutation(treeNames, seqNames, alignment->alignmentSize);
		offspring->score = countScoreHash(alignment, offspring->tree, pwmMatrix, alpha, gapOpt,\
								hashScore, permutation);

		printf("Iter: %4d, Score: %ld, parents: %3d and %3d, ", t + 1, offspring->score, i + 1, j + 1);

		if (!(treeIsUnique(offspring, -1, population, k)))
		{
			failure++;
			printf("unique: 0 ");
		}
		if (treeIsUnique(offspring, -1, population, k))
		{
			if (tmpFailureStreak < failure)
			{
				tmpFailureStreak = failure;
			}
			failure = 0;
			population[k] = offspring;
			treeWithScoreSort(population, k);
			printf("unique: 1 ");
			k++;
		}

		if (population[k - 1]->score != leaderScore)
		{
			leaderScore = population[k - 1]->score;
			if (maxLeaderStreak < leaderIter)
			{
				maxLeaderStreak = leaderIter;
			}
			leaderIter = 0;
			if (maxFailureStreak < tmpFailureStreak)
			{
				maxFailureStreak = tmpFailureStreak;
			}
			tmpFailureStreak = 0;
			printf("new leader: yes\n");
		}
		else
		{
			leaderIter++;
			printf("new leader: no\n");
		}


		if (failure == iterLim)
		{
			printf("No new trees for %d iterations\n", iterLim);
			printf("Longest failure streak: %d, longest leader streak: %d\n", maxFailureStreak, maxLeaderStreak);
			return population[k - 1];
		}

		if (leaderIter == iterNew) 
		{
			printf("Leader remains for %d iterations\n", iterNew);
			printf("Longest failure streak: %d, longest leader streak: %d\n", maxFailureStreak, maxLeaderStreak);
			return population[k - 1];
		}

		t++;

		if (k == treeNum)
		{
			printf("unique population successfully generated, starting genitor\n");
			for (i = 0; i < treeNum; i++)
			{
				treeWithScoreDelete(initPop[i]);
			}
			free(initPop);
			break;
		}
	}

	treeWithScoreSort(population, treeNum); 
	leaderScore = population[treeNum - 1]->score; 

	while (t < iterNum)
	{
		i = rand() % treeNum;
		j = rand() % treeNum;
		while (i == j) 
		{
			j = rand() % treeNum;
		}

		offspring = crossover(population[i], population[j], alignment, alpha,\
							gapOpt, pwmMatrix, hashScore);

		treeNames = treeGetNames(offspring->tree);
        permutation = calculatePermutation(treeNames, seqNames, alignment->alignmentSize);
		offspring->score = countScoreHash(alignment, offspring->tree, pwmMatrix, alpha, gapOpt,\
								hashScore, permutation);

		printf("Iter: %4d, Score: %ld, parents: %3d and %3d, ", t + 1, offspring->score, i + 1, j + 1);

		if (offspring->score <= population[0]->score)
		{
			failure++;
			printf("bad score ");
		}
		else if (!(treeIsUnique(offspring, -1, population, treeNum)))
		{
			failure++;
			printf("unique: 0 ");
		}
		if ((offspring->score > population[0]->score) && 
                   (treeIsUnique(offspring, -1, population, treeNum)))
		{
			if (tmpFailureStreak < failure)
			{
				tmpFailureStreak = failure;
			}
			failure = 0;
			population[0] = offspring;
			treeWithScoreSort(population, treeNum);
			printf("unique: 1 ");
		}

		if (population[treeNum - 1]->score != leaderScore)
		{
			leaderScore = population[treeNum - 1]->score;
			if (maxLeaderStreak < leaderIter)
			{
				maxLeaderStreak = leaderIter;
			}
			leaderIter = 0;
			if (maxFailureStreak < tmpFailureStreak)
			{
				maxFailureStreak = tmpFailureStreak;
			}
			tmpFailureStreak = 0;
			printf("new leader: yes\n");
		}
		else
		{
			leaderIter++;
			printf("new leader: no\n");
		}


		if (failure == iterLim)
		{
			printf("No new trees for %d iterations\n", iterLim);
			printf("Longest failure streak: %d, longest leader streak: %d\n", maxFailureStreak, maxLeaderStreak);
			return population[treeNum - 1];
		}

		if (leaderIter == iterNew) 
		{
			printf("Leader remains for %d iterations\n", iterNew);
			printf("Longest failure streak: %d, longest leader streak: %d\n", maxFailureStreak, maxLeaderStreak);
			return population[treeNum - 1];
		}

		t++;
	}

	free(uniqueSet);
	printf("Genitor finished, returning leader\n");
	printf("Longest failure streak: %d, longest leader streak: %d\n", maxFailureStreak, maxLeaderStreak);
	return population[treeNum - 1];
}