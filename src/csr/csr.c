#include "csr.h"

// Load CSR Matrix from CSV file
void loadcsr(char filename[50], int vertices, int nnzCount, int** nnz, int** colIDs, int** rowDel){

    // Open file for reading
    FILE *fp;
    fp = fopen(filename, "r");

    // Allocate temporary line buffer and temporary character variable
    char* lineBuff = calloc((nnzCount * 2) * 10, sizeof(char));
    char* token;

    // Position pointer to array
    int arrayPtr = 0;

    // Delimiter for strtok()
    char delimiter[1]= ",";

    // Allocate arrays to write
    (*nnz) = malloc(nnzCount * sizeof(int));
    (*colIDs) = malloc(nnzCount * sizeof(int));
    (*rowDel) = malloc((vertices + 1) * sizeof(int));

    // Read nnz
    fgets(lineBuff, (nnzCount * 2) * 10, fp);
    char* tmp = strdup(lineBuff);
    token = strtok(tmp, delimiter);
    while((token != NULL) && (arrayPtr < nnzCount) && (*token != '\n')){
        (*nnz)[arrayPtr] = atoi(token);
        arrayPtr++;
        token = strtok(NULL, delimiter);
    }
    free(tmp);

    // Read colIDs
    arrayPtr = 0;
    fgets(lineBuff, (nnzCount * 2) * 10, fp);
    tmp = strdup(lineBuff);
    token = strtok(tmp, delimiter);
    while((token != NULL) && (arrayPtr < nnzCount) && (*token != '\n')){
        (*colIDs)[arrayPtr] = atoi(token);
        arrayPtr++;
        token = strtok(NULL, delimiter);
    }
    free(tmp);

    // Read rowDel
    arrayPtr = 0;
    fgets(lineBuff, (vertices * 2) * 10, fp);
    tmp = strdup(lineBuff);
    token = strtok(tmp, delimiter);
    while((token != NULL) && (arrayPtr < vertices + 1) && (*token != '\n')){
        (*rowDel)[arrayPtr] = atoi(token);
        arrayPtr++;
        token = strtok(NULL, delimiter);
    }
    free(tmp);

    free(lineBuff);
}