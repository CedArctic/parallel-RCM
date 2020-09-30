#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include "src/rcm/rcm.h"
#include "src/rcm_parallel/rcm_parallel.h"
#include "src/csr/csr.h"
#include "src/graph/graph.h"
#include "config.h"

int main() {

    // Configure OpenMP threads number
    omp_set_num_threads(2);

    // Load matrix in CSR format
    int *nnz, *colIDs, *rowDel;
    int nnzCount = NNZ_COUNT;
    loadcsr("csr.csv", GRAPH_DIM, nnzCount, &nnz, &colIDs, &rowDel);

    // Variables for benchmarking
    clock_t start, end, secStart;
    double sections[5];
    double timings[3] = {0};
    double avg = 0;

    // Variable declarations for parallel approaches
    int height, width, *levelCounts;

    int* permutation2;
    int* permutation_seq;

    // Run 5 tests
    for(int i = 0; i < 5; i++) {

        // Clear memory from previous loop
        if(i > 0) {
            free(permutation_seq);
        }

        // === Sequential Approach ===
        // Get permutation using Sequential RCM
        start = clock();
        permutation_seq = rcm(GRAPH_DIM, colIDs, rowDel);
        end = clock();
        timings[0] += ((double) (end - start)) * 1000 / CLOCKS_PER_SEC;


        // === Parallel Approaches
        // Build graph
        start = clock();
        node *graph = buildGraph(GRAPH_DIM, colIDs, rowDel);
        end = clock();
        avg = ((double) (end - start)) * 1000 / CLOCKS_PER_SEC;
        sections[0] = avg;

        // Determine root node
        int *startEndNodes = pseudo_diameter(graph, GRAPH_DIM, CHUNK);

        // === Unified Approach ===
        // Create a copy of graph to pass to rcmParallel()
        node *inputGraph = malloc(GRAPH_DIM * sizeof(node));
        memcpy(inputGraph, graph, GRAPH_DIM * sizeof(node));
        start = clock();
        permutation2 = rcmParallelFull(inputGraph, GRAPH_DIM, startEndNodes[0], CHUNK);
        end = clock();
        timings[2] += avg + ((double) (end - start)) * 1000 / CLOCKS_PER_SEC;

        // Free loop memory
        for(int j = 0; j < GRAPH_DIM; j ++)
            free(graph[j].neighbors);
        free(graph);
        free(startEndNodes);
        free(levelCounts);

    }

    printPermutation(permutation2,"unified.csv", GRAPH_DIM);
    printPermutation(permutation_seq,"sequential.csv", GRAPH_DIM);
    printf("Sequential Execution Time: %f ms\n", timings[0]/5);
    printf("Unified Execution Time: %f ms\n", timings[2]/5);

    // Memory freeing
    free(nnz);
    free(colIDs);
    free(rowDel);
    free(permutation2);
    free(permutation_seq);

    return 0;
}
