#include "quicksort.h"

void quickSortNodes(node** arr, int low, int high) {
    if (low < high) {
        // idx is the partition point
        int idx = partitionNodes(arr, low, high);

        // Divide and conquer
        quickSortNodes(arr, low, idx - 1);
        quickSortNodes(arr, idx + 1, high);
    }
}

// QuickSort Partition function. low and high are the range of indexes in arr where partition should work
int partitionNodes(node** arr, int low, int high) {

    // Select a pivot and initialize flag to position of smallest element before pivot
    int pivot = arr[high]->neighborsNum;
    int i = (low - 1);

    // Go through the array examining each element
    for (int j = low; j <= high - 1; j++) {
        // If current element is smaller than the pivot, increment i and swap it out with the one currently pointed by i
        if (arr[j]->neighborsNum < pivot) {
            i++;
            // Swap distances and corresponding point ids
            swap_nodes(&arr[i], &arr[j]);
        }
    }

    // Finally place pivot in its correct position in the array and return the position as the middle point
    swap_nodes(&arr[i + 1], &arr[high]);
    return (i + 1);
}

// A utility function to swap two elements
void swap_nodes(node **a, node **b) {
    node* t = *a;
    *a = *b;
    *b = t;
}