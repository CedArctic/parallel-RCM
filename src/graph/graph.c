#include "graph.h"

// Returns an array of graph nodes built from a csr format matrix
node* buildGraph(int vertices, const int* colIDs, const int* rowDel){

    // Allocate memory for graph array
    node* graph = calloc(vertices, sizeof(node));

    // Go through nodes and initialize them
    #pragma omp parallel shared(graph, vertices, colIDs, rowDel) default(none)
    {
            #pragma omp single
            for (int i = 0; i < vertices; i++) {
                graph[i].id = i;
                graph[i].level = INT_MAX;

                // Count neighbors
                graph[i].neighborsNum = rowDel[i + 1] - rowDel[i];

                // Add neighbors to the neighbors list
                graph[i].neighbors = malloc(graph[i].neighborsNum * sizeof(int));
            }
            #pragma omp for
            for (int i = 0; i < vertices; i++) {
                for (int k = rowDel[graph[i].id], elCounter = 0; k < rowDel[graph[i].id + 1]; k++) {
                    if (colIDs[k] != graph[i].id) {
                        graph[i].neighbors[elCounter] = colIDs[k];
                        elCounter++;
                    }else{
                        graph[i].neighborsNum--;
                    }
                }

            }
    }
    return graph;
}

node* buildGraphSequential(int vertices, const int* colIDs, const int* rowDel){

    // Allocate memory for graph array
    node* graph = calloc(vertices, sizeof(node));

    // Go through nodes and initialize them
    for(int i = 0; i < vertices; i++){
        graph[i].id = i;
        graph[i].level = INT_MAX;

        // Count neighbors
        graph[i].neighborsNum = rowDel[i+1] - rowDel[i];

        // Add neighbors to the neighbors list
        graph[i].neighbors = malloc(graph[i].neighborsNum * sizeof(int));
        int elCounter = 0;
        for(int k = rowDel[graph[i].id]; k <  rowDel[graph[i].id + 1]; k++){
            if(colIDs[k] == graph[i].id){
                graph[i].neighborsNum--;
            }
            if(colIDs[k] != graph[i].id){
                graph[i].neighbors[elCounter] = colIDs[k];
                elCounter++;
            }
        }

    }

    return graph;

}

// Return a list of nodes-neighbors, children, and the count of children of a node that are on a specified level
int neighborsAtLevel(node* graph, node* n, int level, node** children){

    // Count children
    int childrenCount = 0;

    // Place children in the children array
    for(int j = 0; j < n->neighborsNum; j++){
        if(graph[n->neighbors[j]].level == level){
            children[childrenCount++] = &graph[n->neighbors[j]];
        }
    }

    // Return children count
    return childrenCount;
}

// Get the height, width and number of nodes on each level of the graph
void graphHeightWidthCounts(node* graph, int vertices, int* height, int* width, int** levelCounts){

    // Allocate shared arrays
    int maxThreads = omp_get_max_threads();
    int* localcount = calloc(maxThreads * vertices, sizeof(int));
    int* localmax = calloc(maxThreads, sizeof(int));
    #pragma omp parallel shared(localcount, localmax, graph, vertices, height, levelCounts) default(none)
    {
        int tid = omp_get_thread_num();
        int teamThreads = omp_get_team_size(omp_get_level());
        #pragma omp for schedule(guided) nowait
        for(int i = 0; i < vertices; i++){
            localcount[tid * vertices + graph[i].level]++;
            if(localmax[tid] < graph[i].level){
                localmax[tid] = graph[i].level;
            }
        }

        #pragma omp single
        {
            // Determine height (max_level)
            *height = 0;
            for (int k = 0; k < teamThreads; k++) {
                if (*height < localmax[k]) {
                    *height = localmax[k];
                }
            }
            // Form counts
            *levelCounts = calloc(*height + 1, sizeof(int));
        }

        #pragma omp for schedule(guided) nowait
        for(int i = 0; i < *height + 1; i++){
            for(int id = 0; id < teamThreads; id++){
                (*levelCounts)[i] += localcount[id * vertices + i];
            }
        }
    }
    // Determine width (maximum number of nodes on a level)
    *width = *levelCounts[0];
    for(int j = 1; j < *height + 1; j++){
        if((*levelCounts)[j] > *width){
            *width = (*levelCounts)[j];
        }
    }

    // Free memory
    free(localcount);
    free(localmax);
}

// Returns a pointer to an array with the nodes at a specific level.
node** graphVerticesAt(node* graph, int vertices, int level, int levelNodes){

    // Allocate returned array
    node** verticesAtLevel = calloc(levelNodes, sizeof(node*));
    int verticesAtLevelIdx = 0;

    // Scan graph adding vertices at level to the results array
    for(int i = 0; i < vertices; i++){
        if(graph[i].level == level){
            verticesAtLevel[verticesAtLevelIdx] = &graph[i];
            verticesAtLevelIdx++;
        }
    }
    return verticesAtLevel;
}

// Returns an array with the prefix sums of the provided counts array
int* prefixSums(int* counts, int height){

    // Create a copy of the provided array to act upon and return
    // maxLevel + 1 because levels start from 0
    int* prefix_sums;
    prefix_sums = malloc((height + 1) * sizeof(int));
    memcpy(prefix_sums, counts, (height + 1) * sizeof(int));

    // Check if work is worth done sequentially
    if(height < PREFIX_SUMS_SEQ_LIMIT){
        prefixSumsSeq(prefix_sums, height + 1);
        return prefix_sums;
    }

    // Get and set maximum number of threads and calculate work per thread and number of switching stages
    int maxThreads = PREFIX_SUMS_THREADS;
//    omp_set_num_threads(maxThreads);
    // We assume a system with threads at a power of 2
    int num_changes = (int)log2(maxThreads);

    // Chunk up the work
    int chunk = (height + 1) / maxThreads;
    int remainder = (height + 1) % maxThreads;
    int workIndexes[maxThreads + 1];
    workIndexes[0] = 0;
    for(int i = 1, sum = 0; i < maxThreads + 1; i++){
        sum += chunk;
        if(i <= remainder)
            sum++;
        workIndexes[i] = sum;
    }

    // Parallel calculate local prefix sums
    int cPrefix[maxThreads], cTotal[maxThreads], lPrefix[maxThreads], lTotal[maxThreads];
    #pragma omp parallel shared(prefix_sums, workIndexes, cPrefix, cTotal, lPrefix, lTotal, maxThreads, num_changes) \
    default(none) num_threads(PREFIX_SUMS_THREADS)
    {
        int tid = omp_get_thread_num();
        int elementsNum = workIndexes[tid+1] - workIndexes[tid];
        prefixSumsSeq(prefix_sums + workIndexes[tid], elementsNum);
        cPrefix[tid] = cTotal[tid] = prefix_sums[workIndexes[tid+1] - 1];
        lPrefix[tid] = lTotal[tid] = prefix_sums[workIndexes[tid+1] - 1];

        #pragma omp barrier
        for(int i = 0; i < num_changes; i++){
            int tidn = tid ^ ((int)pow(2, i));
            if((tidn < maxThreads) && (tidn != tid)){
                if(tidn > tid){
                    lPrefix[tidn] += cTotal[tid];
                    lTotal[tidn] += cTotal[tid];
                }else{
                    lTotal[tidn] += cTotal[tid];
                }
            }
            // Sync threads before writing to the cross arrays (cTotal/cPrefix)
            #pragma omp barrier
            cPrefix[tid] = lPrefix[tid];
            cTotal[tid] = lTotal[tid];
            #pragma omp barrier
        }

        if(tid != 0){
            for(int j = workIndexes[tid]; j < workIndexes[tid+1]; j++){
                prefix_sums[j] += cPrefix[tid - 1];
            }
        }
    }

    return prefix_sums;
}

// Sequential prefix sums on array
void prefixSumsSeq(int* array, int elements){
    int sum = 0;
    for(int i = 0; i < elements; i++){
        sum += array[i];
        array[i] = sum;
    }
}