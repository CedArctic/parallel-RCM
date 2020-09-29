#ifndef RCM_CONFIG_H
#define RCM_CONFIG_H

// NOTE: CONFIGURE GRAPH_DIM AND NNZ_COUNT OR ELSE THINGS WON'T WORK
#define CHUNK 0.75
#define GRAPH_DIM 221119
#define NNZ_COUNT 7666057
#define PREFIX_SUMS_SEQ_LIMIT 32
#define MAX_THREADS 2
#define PREFIX_SUMS_THREADS 1

#endif //RCM_CONFIG_H