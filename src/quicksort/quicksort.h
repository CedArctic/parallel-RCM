#ifndef RCM_QUICKSORT_H
#define RCM_QUICKSORT_H

#include "../graph/graph.h"
/*
 * Quicksort for nodes based on their degree. Used to order children array during placement of nodes in permutation.
 *
 * arr: Array of node pointers to be sorted based on degrees
 * low: Starting index
 * high: End index (= length of arr - 1 to sort the complete array)
 */
void quickSortNodes(node** arr, int low, int high);
int partitionNodes(node** arr, int low, int high);
void swap_nodes(node **a, node **b);

#endif //RCM_QUICKSORT_H
