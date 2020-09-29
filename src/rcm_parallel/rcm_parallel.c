#include "rcm_parallel.h"

// Shrinking strategy
int shrinkStrategy(node*** candSet, int candSetSize){

    // Find the maximum degrees of the initial candidate set
    int maxDegrees = 0;
    for(int i = 0; i < candSetSize; i++){
        if(maxDegrees < (*candSet)[i]->neighborsNum){
            maxDegrees = (*candSet)[i]->neighborsNum;
        }
    }

    // This is an integer array with a cell for each degree. If a cell is marked with a non zero or one number it means that
    // a node with that degree has already been added to the new candidates set. Length is (maxDegrees + 1) because degree
    // numbering starts from 0
    int* degreeAdded = calloc(maxDegrees + 1, sizeof(int));

    // Count how many nodes will be added to the new set
    int newCandSetSize = 0;
    for(int j = 0; j < candSetSize; j++){
        if(degreeAdded[(*candSet)[j]->neighborsNum] == 0){
            degreeAdded[(*candSet)[j]->neighborsNum] = 1;
            newCandSetSize++;
        }
    }

    // Allocate memory for new shrinked candidate set
    node** newCandSet = calloc(newCandSetSize, sizeof(node*));
    int newCandSetIdx = 0;

    // Scan nodes adding one of each degree
    for(int k = 0; k < candSetSize; k++){
        if(degreeAdded[(*candSet)[k]->neighborsNum] == 1){
            degreeAdded[(*candSet)[k]->neighborsNum] = 0;
            newCandSet[newCandSetIdx] = (*candSet)[k];
            newCandSetIdx++;
        }
    }

    // Free old memory and swap with new
    free(degreeAdded);
    free(*candSet);
    *candSet = newCandSet;
    return newCandSetSize;
}

// Pseudo-optimal way of selecting two nodes further apart optimal for RCM roots
int* pseudo_diameter(node* graph, int vertices, double chunk){

    // Allocate result arrays for bfs searches
    node *forwardBFS, *reverseBFS;

    // Initialize pointers to start and end nodes
    int diameterEnd = -1;
    int diameterStart = graph[0].id;
    // Initialize start node to the node with the smallest degree
    for(int i = 0; i < vertices; i++){
        if(graph[i].neighborsNum < graph[diameterStart].neighborsNum){
            diameterStart = graph[i].id;
        }
    }

    // Declare variables used in the main loop
    int minWidth, localDiameter, forwardWidth, reverseHeight, reverseWidth;
    int *forwardLevelCounts, *reverseLevelCounts;
    node** candSet;

    do {
        forwardBFS = malloc(GRAPH_DIM * sizeof(node));
        memcpy(forwardBFS, graph, GRAPH_DIM * sizeof(node));
        complete_bfs(forwardBFS, vertices, diameterStart, chunk);
        graphHeightWidthCounts(forwardBFS, vertices, &localDiameter, &forwardWidth, &forwardLevelCounts);
        candSet = graphVerticesAt(forwardBFS, vertices, localDiameter, forwardLevelCounts[localDiameter]);
        // Print candidate set
        int shrinkedCandSetSize = shrinkStrategy(&candSet, forwardLevelCounts[localDiameter]);
        minWidth = INT_MAX;
        for(int j = 0; j < shrinkedCandSetSize; j++){
            reverseBFS = malloc(GRAPH_DIM * sizeof(node));
            memcpy(reverseBFS, graph, GRAPH_DIM * sizeof(node));
            complete_bfs(reverseBFS, vertices, candSet[j]->id, chunk);
            graphHeightWidthCounts(reverseBFS, vertices, &reverseHeight, &reverseWidth, &reverseLevelCounts);
            if(reverseWidth < minWidth){
                if(reverseHeight > localDiameter){
                    diameterStart = candSet[j]->id;
                    diameterEnd = -1;
                    break;
                }else{
                    minWidth = reverseWidth;
                    diameterEnd = candSet[j]->id;
                }
            }

            // Memory freeing for loop
            free(reverseLevelCounts);
            free(reverseBFS);
        }

        // Memory freeing for loop
        free(forwardLevelCounts);
        if(diameterEnd == -1){
            free(reverseLevelCounts);
            free(reverseBFS);
        }
        free(candSet);
        free(forwardBFS);

    }while (diameterEnd == -1);

    // Allocate results array and place the start and end nodes inside
    int* results = calloc(2, sizeof(int));
    results[0] = diameterStart;
    results[1] = diameterEnd;

    return results;
}

void complete_bfs(node* graph, const int vertices, int initialRoot, double chunk){
    int relaxedNodes = 0;
    int root = initialRoot;
    do {
        relaxedNodes += bfs(graph, vertices, root, chunk);
        // If BFS hasn't reached all nodes yet, choose a new root to continue the search
        if(relaxedNodes < vertices){
            root = -1;
            for(int j = 0; j < vertices; j++){
//                if(((root == -1) || (graph[root].neighborsNum > graph[j].neighborsNum)) \
//                        && (graph[j].level == INT_MAX)){
                if(graph[j].level == INT_MAX){
                    root = j;
                    break;
                }
            }
        }
    }while(relaxedNodes < vertices);
}

int bfs(node* graph, const int vertices, int root, double chunk){

    // Shared worklist variables
    node** wl;
    wl = calloc(100*vertices, sizeof(node*));
    int wlHead = 0, wlTail = -1;
    int localHead, localTail;
    int sizeChunk;

    // Counter of nodes that have been reached
    int reachedNodes = 0;

    // Set root level
    graph[root].level = 0;

    // Count nodes relaxed
    int relaxedNodes = 1;

    // Enqueue Root node
    wl[++wlTail] = &graph[root];

    #pragma omp parallel \
    private(localHead, localTail, sizeChunk) \
    shared(graph, wl, wlTail, wlHead, reachedNodes, vertices, chunk, relaxedNodes) default(none)
    {
        // Counter for number of nodes relaxed for the first time
        int localRelaxedNodes = 0;

        node** relaxedwl = malloc((int)(chunk * vertices) * sizeof(node*));
        while(((wlTail - wlHead) >= 0)){

            // Shift shared worklist head and increment reached nodes counter
            #pragma omp critical
            {
                localHead = wlHead;
                localTail = wlTail;
                sizeChunk = ceil(chunk * (double)(localTail - localHead + 1));
                wlHead += sizeChunk;
                reachedNodes += sizeChunk;
            }

            // Work chunking. Initialize local worklist and copy node* chunk to it
            node** localwl = wl + localHead;

            // Fixed Point iteration
            int relaxedWlIdx = 0;
            for(int i = 0; i < sizeChunk; i++){
                int level = localwl[i]->level + 1;
                for(int j = 0; j < localwl[i]->neighborsNum; j++){
                    if (level < graph[localwl[i]->neighbors[j]].level) {
                        if(graph[localwl[i]->neighbors[j]].level == INT_MAX){
                            localRelaxedNodes++;
                        }
                        #pragma omp atomic write
                        graph[localwl[i]->neighbors[j]].level = level;
                        relaxedwl[relaxedWlIdx++] = &graph[localwl[i]->neighbors[j]];
                    }

                }
            }

            // Relaxing nodes & Checking for multiple components
            #pragma omp critical
            {
                // Relax Nodes
                for(int k = 0; k < relaxedWlIdx; k++){
                    wl[++wlTail] = relaxedwl[k];
                }
            }

        }
        // Free locally allocated memory before next iteration
        free(relaxedwl);

        #pragma omp atomic update
        relaxedNodes += localRelaxedNodes;
    }

    // Free worklist
    free(wl);

    return relaxedNodes;

}

int cmpNodes (const void * a, const void * b) {
    return ( ((node*)a)->neighborsNum - ((node*)b)->neighborsNum );
}

int* placement(node* graph, int vertices, int source, int levels, int width, int* levelCounts){

    // Allocate results array
    volatile int* permutation;

    // Allocate array indicating if a node has already been added to the permutation
    volatile int* nodeAdded;

    // Allocate and initialize readOffset and writeOffset arrays
    volatile int* readOffset;
    volatile int* writeOffset;

    // Shared children array
    node** children;

    // Non-parallelized initialization (deprecated)
//    // Allocate results array
//    volatile int* permutation = malloc(vertices * sizeof(int));
//
//    // Allocate array indicating if a node has already been added to the permutation
//    volatile int* nodeAdded = calloc(vertices, sizeof(int));
//
//    // Allocate and initialize readOffset and writeOffset arrays
//    volatile int* readOffset = malloc((levels+1) * sizeof(int));
//    volatile int* writeOffset = malloc((levels+1) * sizeof(int));
//    memcpy(readOffset + 1, levelCounts, levels * sizeof(int));
//    memcpy(writeOffset + 1, levelCounts, levels * sizeof(int));
//    readOffset[0] = 0;
//
//    // Add root to the permutation
//    permutation[0] = source;
//    writeOffset[0] = 1;
//    nodeAdded[source] = 1;
//
//    // Also add all other level 0 nodes - for graphs with multiple components
//    // Disabling this might yield better performance on connected graphs
//    for(int i = 0; i < vertices; i++){
//        if((graph[i].level == 0) && (i!=source)){
//            permutation[writeOffset[0]++] = i;
//            nodeAdded[i] = 1;
//        }
//    }
//
//    // Get maximum number of threads
//    int maxThreads = omp_get_max_threads();
//
//    // Shared children array
//    node** children = malloc(maxThreads * width * sizeof(node*));

    // Initialize waiting locks
//    omp_lock_t locks[maxThreads];
//    for (int i = 0; i < maxThreads; i++)
//        omp_init_lock(&locks[i]);

    #pragma omp parallel shared(readOffset, writeOffset, permutation, levels, levelCounts, graph, nodeAdded, width, children, vertices, source) \
    default(none)
    {
        int threadsNum = omp_get_num_threads();

        #pragma omp single nowait
        {
            children = malloc(threadsNum * width * sizeof(node*));
        }

        #pragma omp single nowait
        {
            // Allocate and initialize readOffset and writeOffset arrays
            readOffset = malloc((levels+1) * sizeof(int));
            memcpy(readOffset + 1, levelCounts, levels * sizeof(int));
            readOffset[0] = 0;
        };

        #pragma omp single nowait
        {
            permutation = malloc(vertices * sizeof(int));
            permutation[0] = source;
        };

        #pragma omp single
        {
            nodeAdded = malloc(vertices * sizeof(int));
            // Add root to the permutation
            nodeAdded[source] = 1;
        };

        #pragma omp single
        {
            writeOffset = malloc((levels+1) * sizeof(int));
            memcpy(writeOffset + 1, levelCounts, levels * sizeof(int));
            writeOffset[0] = 1;

            // Also add all other level 0 nodes
            // Disabling this might yield better performance on connected graphs
            for(int i = 0; i < vertices; i++){
                if((graph[i].level == 0) && (i!=source)){
                    permutation[writeOffset[0]++] = i;
                    nodeAdded[i] = 1;
                }
            }
        };

        int tid = omp_get_thread_num();
//        omp_set_lock(&locks[tid]);
        node** localChildren = children + tid*width;
//        #pragma omp barrier
        for(int l = tid; l < levels; l += threadsNum){
            while(readOffset[l] != levelCounts[l]){

//                #pragma omp flush (writeOffset)

                while(readOffset[l] == writeOffset[l]){
//                    omp_set_lock(&locks[(l-1) % threadsNum]);
                }

//                if (l > 0) omp_unset_lock(&locks[(l-1) % threadsNum]);

                node N = graph[permutation[readOffset[l]]];
                ++readOffset[l];

                int childCount = neighborsAtLevel(graph, &N, N.level + 1, localChildren);

//              quickSortNodes(localChildren, 0, childCount - 1);
                qsort(localChildren, childCount, sizeof(node*), cmpNodes);

                for(int c = 0; c < childCount; c++){
                    // Check if child already added
                    if(nodeAdded[localChildren[c]->id] != 1){
//                        omp_test_lock(&locks[l % threadsNum]);

                        nodeAdded[localChildren[c]->id] = 1;
                        permutation[writeOffset[l+1]] = localChildren[c]->id;
                        ++writeOffset[l+1];

//                        omp_unset_lock(&locks[l % threadsNum]);
                    }
                }
            }
        }
    };

    // Memory freeing
    free(children);
    free(readOffset);
    free(writeOffset);
    free(nodeAdded);

    return permutation;

}

int* rcmParallelFull(node* graph, int vertices, int root, double chunk){

    // OpenMP Configuration
    int maxThreads = MAX_THREADS;

    // === BFS ===
    // Shared worklist variables
    node** wl;
    wl = calloc(100*vertices, sizeof(node*));
    int wlHead = 0, wlTail = -1;
    // Set root level
    graph[root].level = 0;
    // Enqueue Root node
    wl[++wlTail] = &graph[root];
    // Total number of nodes that have been relaxed in BFS
    int relaxedNodes = 0;
    // Temporary variable used for multi-component graph bfs searches
    int currentRoot = root;

    // === Graph Height, Width, Counts
    // Allocate shared arrays
    int* localcount = calloc(maxThreads * vertices, sizeof(int));
    int* localmax = calloc(maxThreads, sizeof(int));
    int height, width;
    int* levelCounts;

    // === Prefix Sums ===
    // Parallel calculate local prefix sums
    int workIndexes[maxThreads + 1];
    int cPrefix[maxThreads];
    int cTotal[maxThreads];
    int lPrefix[maxThreads];
    int lTotal[maxThreads];
    workIndexes[0] = 0;

    // === Placement ===
    // Allocate results array
    volatile int* permutation = malloc(vertices * sizeof(int));
    // Allocate array indicating if a node has already been added to the permutation
    volatile int* nodeAdded = calloc(vertices, sizeof(int));
    volatile int *readOffset, *writeOffset;

    #pragma omp parallel \
    shared(graph, wl, wlTail, wlHead, vertices, chunk, localcount, localmax, height, width, maxThreads, \
    levelCounts, cPrefix, cTotal, lPrefix, lTotal, permutation, readOffset, writeOffset, nodeAdded, root, workIndexes, \
    currentRoot, relaxedNodes) \
    default(none) num_threads(MAX_THREADS)
    {

        // === BFS ===
        do {

            #pragma omp single
            {
                wlHead = 0;
                wlTail = -1;
                // Set root level
                graph[currentRoot].level = 0;
                // Enqueue Root node
                wl[++wlTail] = &graph[currentRoot];
                relaxedNodes++;
            }

            int localHead, localTail, sizeChunk;

            // Counter for number of nodes relaxed for the first time
            int localRelaxedNodes = 0;

            node** relaxedwl = malloc((int)(chunk * vertices) * sizeof(node*));

            while(((wlTail - wlHead) >= 0)){

                // Shift shared worklist head and increment reached nodes counter
                #pragma omp critical
                {
                    localHead = wlHead;
                    localTail = wlTail;
                    sizeChunk = ceil(chunk * (double)(localTail - localHead + 1));
                    wlHead += sizeChunk;
                }

                // Work chunking. Initialize local worklist and copy node* chunk to it
                node** localwl = wl + localHead;

                // Fixed Point iteration
                int relaxedWlIdx = 0;
                for(int i = 0; i < sizeChunk; i++){
                    int level = localwl[i]->level + 1;
                    for(int j = 0; j < localwl[i]->neighborsNum; j++){
                        if(level < graph[localwl[i]->neighbors[j]].level){
                            if(graph[localwl[i]->neighbors[j]].level == INT_MAX){
                                localRelaxedNodes++;
                            }
                            #pragma omp atomic write
                            graph[localwl[i]->neighbors[j]].level = level;
                            relaxedwl[relaxedWlIdx++] = &graph[localwl[i]->neighbors[j]];
                        }
                    }
                }

                // Relaxing nodes
                #pragma omp critical
                {
                    // Relax Nodes
                    for(int k = 0; k < relaxedWlIdx; k++){
                        wl[++wlTail] = relaxedwl[k];
                    }
                }

            }
            // Free locally allocated memory before next iteration
            free(relaxedwl);

            #pragma omp atomic update
            relaxedNodes += localRelaxedNodes;
            //TODO: Remove debug?
#pragma omp barrier
            // If BFS hasn't reached all nodes yet, choose a new root to continue the search
            #pragma omp single
            {
                if(relaxedNodes < vertices){
                    currentRoot = -1;
                    for(int j = 0; j < vertices; j++){
                        if(graph[j].level == INT_MAX){
                            currentRoot = j;
                            break;
                        }
                    }
                }
            };
        }while(relaxedNodes < vertices);


        // === Height, Width, Counts ===
        //TODO: Remove debug?
#pragma omp barrier
        // Warning: Segfaults occur if bellow barrier is removed
        int tid = omp_get_thread_num();
        #pragma omp barrier
        #pragma omp for schedule(auto) nowait
        for(int i = 0; i < vertices; i++){
            localcount[tid * vertices + graph[i].level]++;
            if(localmax[tid] < graph[i].level){
                localmax[tid] = graph[i].level;
            }
        }

        #pragma omp single
        {
            // Determine height (max_level)
            height = 0;
            for (int k = 0; k < maxThreads; k++) {
                if (height < localmax[k]) {
                    height = localmax[k];
                }
            }
            // Form counts
            levelCounts = calloc(height + 1, sizeof(int));
        }

        // TODO: Decide on a schedule - guided seems good, dynamic is bad
        //TODO: Remove debug?
#pragma omp flush(height)
        #pragma omp for schedule(auto)
        for(int i = 0; i < height + 1; i++){
            for(int id = 0; id < maxThreads; id++){
                (levelCounts)[i] += localcount[id * vertices + i];
            }
        }

        #pragma omp single nowait
        {
            // Determine width (maximum number of nodes on a level)
            width = levelCounts[0];
            for(int j = 1; j < height + 1; j++){
                if((levelCounts)[j] > width){
                    width = (levelCounts)[j];
                }
            }
        };

        // === Prefix Sums ===
        int num_changes = (int)log2(maxThreads);
        //TODO: Remove debug?
#pragma omp flush(height)
#pragma omp barrier
        // Distribute work to threads
        #pragma omp single
        {
            int chunk = (height + 1) / maxThreads;
            int remainder = (height + 1) % maxThreads;
            for(int i = 1, sum = 0; i < maxThreads + 1; i++){
                sum += chunk;
                if(i <= remainder)
                    sum++;
                workIndexes[i] = sum;
            }
        };

        int elementsNum = workIndexes[tid+1] - workIndexes[tid];
        prefixSumsSeq(levelCounts + workIndexes[tid], elementsNum);
        cPrefix[tid] = cTotal[tid] = levelCounts[workIndexes[tid+1] - 1];
        lPrefix[tid] = lTotal[tid] = levelCounts[workIndexes[tid+1] - 1];

        int tidn;
        for(int i = 0; i < num_changes; i++){
            tidn = tid ^ ((int)pow(2, i));
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
                levelCounts[j] += cPrefix[tid - 1];
            }
        }


        // === Placement ===
        //TODO: Remove debug?
#pragma omp barrier
        int levels = height + 1;
        #pragma omp single nowait
        {
            // Allocate and initialize readOffset and writeOffset arrays
            readOffset = malloc((levels+1) * sizeof(int));
            memcpy(readOffset + 1, levelCounts, levels * sizeof(int));
            readOffset[0] = 0;
        };

        #pragma omp single nowait
        {
            // Add root to the permutation
            permutation[0] = root;
            nodeAdded[root] = 1;
        };

        #pragma omp single
        {
            writeOffset = malloc((levels+1) * sizeof(int));
            memcpy(writeOffset + 1, levelCounts, levels * sizeof(int));
            writeOffset[0] = 1;

            // Also add all other level 0 nodes - for multi-component graphs
            // Disabling this might yield better performance on connected graphs
            for(int i = 0; i < vertices; i++){
                if((graph[i].level == 0) && (i!=root)){
                    permutation[writeOffset[0]++] = i;
                    nodeAdded[i] = 1;
                }
            }
        };

        volatile node** children = malloc(width * sizeof(node*));
        int threadsNum = omp_get_team_size(omp_get_level());
        for(int l = tid; l < levels; l += threadsNum){
            while(readOffset[l] < levelCounts[l]){
                while(readOffset[l] >= writeOffset[l]){}
                node N = graph[permutation[readOffset[l]]];
                ++readOffset[l];
                int childCount = neighborsAtLevel(graph, &N, N.level + 1, children);
                quickSortNodes(children, 0, childCount-1);
                for(int c = 0; c < childCount; c++){
                    // Check if child already added
                    if(nodeAdded[children[c]->id] != 1){
                        nodeAdded[children[c]->id] = 1;
                        permutation[writeOffset[l+1]] = children[c]->id;
                        ++writeOffset[l+1];
                    }
                }

            }
        }
        // Free loop memory
        free(children);

    }

    // Free worklist (bfs)
    free(wl);

    // Free memory (height, width, counts)
    free(localcount);
    free(localmax);

    // Free memory (placement)
    free(readOffset);
    free(writeOffset);
    free(nodeAdded);

    return permutation;
}