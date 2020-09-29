#ifndef RCM_RCM_H
#define RCM_RCM_H

#include <stdlib.h>
#include <stdio.h>
#include "../quicksort/quicksort.h"
#include "../graph/graph.h"

/*
 * Function Signatures
 */
int* rcm(int vertices, int* colIDs, int* rowDel);
void printPermutation(int* permutation, char name[50], int length);

#endif //RCM_RCM_H
