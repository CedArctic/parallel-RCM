#ifndef RCM_GRAPH_H
#define RCM_GRAPH_H

#include <stdlib.h>
#include <limits.h>
#include <omp.h>
#include <math.h>
#include <string.h>
#include "../../config.h"

// Node Structure
typedef struct {
    // Node ID, its level and number of neighbors (degrees)
    int id, level, neighborsNum;
    // Array of neighbors
    int *neighbors;
} node;

/*
 * Returns an array of graph nodes built from a csr format matrix
 *
 * vertices: dimension of adjacency matrix
 * nnzCount: length of nnz and colIDs
 * colIDs: ontains the column ids of the nnz elements
 * rowDel: row delimiters of the nnz and colIDs arrays. It's length is vertices+1
 */
node* buildGraph(int vertices, const int* colIDs, const int* rowDel);
node* buildGraphSequential(int vertices, const int* colIDs, const int* rowDel);

/*
 * Return a list of nodes-neighbors, children, and the count of children of a node that are on a specified level
 *
 * graph: array of nodes
 * n: node of which the children at level will be returned
 * level: level
 * children: pointer to an array of node pointers where the results (children) will be placed
 * int returned: number of neighbors at level
 */
int neighborsAtLevel(node* graph, node* n, int level, node** children);

/*
 * Returns the height, width (maximum number of nodes on a level) and number of nodes on each level of the graph
 *
 * graph: array of nodes
 * vertices: the number of nodes
 * height, width: pointers to place the results
 * levelCounts: pointer that will point to an array containing the number of nodes on each level. Level number is the
 *              array index. Length of the array is height+1 (since height starts counting from 0 for the root)
 */
void graphHeightWidthCounts(node* graph, int vertices, int* height, int* width, int** levelCounts);

/*
 * Returns a pointer to an array of pointers to the nodes at a specific level.
 * Requires running graphHeightWidthCounts() first to count the number of nodes on the desired level. The array returned
 * has the length that occurs from levelCounts of graphHeightWidthCounts().
 *
 * graph: array of nodes
 * vertices: the number of nodes
 * level: the level of which we want to find the nodes
 * levelNodes: number of the nodes that will be returned. This can be found by first running graphHeightWidthCounts() and
 *             looking at levelCounts[level]
 * WARNING: The nodes in the return array are copies of the ones given in graph, however each node contains a pointer
 * to a list of neighbors. This means that both original and copy nodes point to the same space in memory for their neighbors.
 * This shouldn't be an issue since the neighbors array is constructed once, before graphVertices() is run, and never altered.
 */
node** graphVerticesAt(node* graph, int vertices, int level, int levelNodes);

/*
 * Returns an array with the prefix sums of the provided counts array
 *
 * counts: input array to the prefix sum algorithm
 * height: the maximum level in counts (levels start from 0 so the length of counts is height + 1)
 */
int* prefixSums(int* counts, int height);

/*
 * Sequential prefix sums on array
 *
 * array: input array
 * elements: length of array
 */
void prefixSumsSeq(int* array, int elements);

#endif //RCM_GRAPH_H
