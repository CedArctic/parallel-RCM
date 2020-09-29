#ifndef RCM_CSR_H
#define RCM_CSR_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
 * Load CSR formatted matrix from csv
 *
 * filename: filename of the csr formatted csv file
 * vertices: dimension of matrix
 * nnzCount: number of non zero elements
 * nnz: pointer to int pointer, at the end it contains the non zero matrix elements
 * colIDs: pointer to int pointer, at the end it contains the column ids of the nnz elements
 * rowDel: row delimiters of the nnz and colIDs arrays. It's length is vertices+1
 *
 * WARNING: If the buffers and line reads are not wide enough (character-wise) elements will not be read and issues will
 * arise (such as uninitialized memory etc) in other functions
 */
void loadcsr(char filename[50], int vertices, int nnzCount, int** nnz, int** colIDs, int** rowDel);

#endif //RCM_CSR_H
