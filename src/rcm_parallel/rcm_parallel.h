#ifndef RCM_RCM_PARALLEL_H
#define RCM_RCM_PARALLEL_H

#include <omp.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "../rcm/rcm.h"
#include "../graph/graph.h"
#include "../../config.h"

/*
 * Breadth first search
 *
 * graph: Array of node structs on which bfs will be conducted
 * vertices: The number of elements in the graph array
 * root: The index number in the graph array of the node to be considered as the root node
 * CHUNK: Double number used to chunk up the work in a multithread environment
 * Returns: the number of nodes in the graph visited by the BFS algorithm using the provided root node
 *
 * complete_bfs does the same breadth first search but makes sure that all nodes are reached. Used in disconnected graphs
*/
void complete_bfs(node* graph, int vertices, int initialRoot, double chunk);
int bfs(node* graph, int vertices, int root, double chunk);


/*
 * Pseudo-optimal way of selecting two nodes further apart optimal for RCM roots. Returns an array with two node ids / indexes
 * on the provided graph
 *
 * graph: Array of node structs that are the base graph
 * vertices: The number of elements in the graph array
 * CHUNK: Double number used to chunk up the work in a multithread environment
 *
 */
int* pseudo_diameter(node* graph, int vertices, double chunk);

/*
 * Shrinking strategy: Accepts a pointer to the candidate node array and reallocates it to point to a reduced set
 * according to the shrinking strategy and then returns the size of the new set
 *
 * candSet: pointer to the candidates array - pointer to an array of node pointers
 * candSetSize: length of the above mentioned array
 *
 * Shrinking strategy: Keep one node of each degree.
 *
 */
int shrinkStrategy(node*** candSet, int candSetSize);

/*
 * Placement algorithm to form the final permutation array
 *
 * graph: node array where bfs has been performed
 * vertices: number of nodes in the graph
 * source: index of the root node in the graph found using pseudo_diameter()
 * levels: number of levels of the graph (graph height + 1)
 * levelCounts: prefix sums array of levels. Is of length 'levels'
 *
 */
int* placement(node* graph, int vertices, int source, int levels, int width, int* levelCounts);

/*
 * Unified RCM, contains all the steps in one function. Returns a permutation array with length 'vertices'
 *
 * graph: input node array
 * vertices: number of nodes in the graph
 * root: The index number in the graph array of the node to be considered as the root node
 * CHUNK: Double number used to chunk up the work in a multithread environment
 */
int* rcmParallelFull(node* graph, int vertices, int root, double chunk);

/*
 * Comparator function for qsort.
 *
 * Receives two pointers to nodes and returns which node has the highest degree.
 */
int cmpNodes (const void * a, const void * b);

#endif //RCM_RCM_PARALLEL_H
