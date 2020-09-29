#include "rcm.h"

// Revese Cuthill McKee method
int* rcm(int vertices, int* colIDs, int* rowDel){

    // Initialize permutation array (results)
    int* permutation = calloc(vertices, sizeof(int));

    // PermCounter counts how many elements are currently in the permutation
    int permCounter = 0;

    // Initialize array indicating vIDs that have been used already
    int *vIDsUsed = calloc(vertices, sizeof(int));

    // Build graph
    node* graph = buildGraphSequential(vertices, colIDs, rowDel);

    // Select root and add it to the permutation
    int minDegreeNode = 0;
    for(int i = 1; i < vertices; i++){
        if(graph[i].neighborsNum < graph[minDegreeNode].neighborsNum){
            minDegreeNode = i;
        }
    }
    permutation[permCounter++] = minDegreeNode;
    vIDsUsed[minDegreeNode] = 1;

    // External loop - loops more than once in the case of disconnected graphs in the same adjacency matrix
    for(int k = 0; k < vertices; k++) {

        // Get the last added node's neighbors
        int neighborsNum = graph[permutation[k]].neighborsNum;
        node** neighbors = malloc(neighborsNum * sizeof(node*));
        for(int i = 0; i < neighborsNum; i++){
            neighbors[i] = &(graph[graph[permutation[k]].neighbors[i]]);
        }

        // Quicksort the neighbors based on their degree
        if(neighborsNum > 0)
            quickSortNodes(neighbors, 0, neighborsNum - 1);

        // Add the neighbors to the permutation if they are not already in it
        for(int j = 0; ((j < neighborsNum) && (permCounter < vertices)); j++){
            if(vIDsUsed[neighbors[j]->id] == 0){
                permutation[permCounter++] = neighbors[j]->id;
                vIDsUsed[neighbors[j]->id] = 1;
            }
        }

        free(neighbors);

        // In case of disconnected graph, add next root. Look for it in the sorted by degrees array by just picking the
        // next available one
        if(((k+1) == permCounter) && ((k + 1) != vertices)){
            int nextRoot = -1;
            for(int i = 0; i < vertices; i++){
                if(vIDsUsed[i] == 0){
                    nextRoot = i;
                    break;
                }
            }
            permutation[permCounter++] = graph[nextRoot].id;
            vIDsUsed[graph[nextRoot].id] = 1;
        }
    }

    free(vIDsUsed);
    for(int i = 0; i< vertices; i ++)
        free(graph[i].neighbors);
    free(graph);

    return permutation;

}

// Write permutation vector to permutation.txt
void printPermutation(int* permutation, char name[50], int length){
    // Open file for reading
    FILE *fp;
    fp = fopen(name, "a");
//    fprintf(fp, "\n%s Permutation Vector (Copy-Paste into Matlab):\np = [", name);
    for(int i = 0; i < length; i++){
        fprintf(fp, "%d", permutation[i] + 1);
        if(i == length - 1)
            break;
        fprintf(fp, ", ");
    }
//    fprintf(fp, "];\n");
fclose(fp);
}