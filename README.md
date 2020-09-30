# parallel-RCM
A parallel Reverse Cuthill-McKee algorithm implementation using OpenMP. Based on [the paper of Rodriguez et al](https://ieeexplore.ieee.org/document/8104594)

## Usage
- Download the latest [release](https://github.com/CedArctic/parallel-RCM/releases)
- Convert the Sparse Matrix to CSR format (currently the CO matrix is included)
- Place the converted matrix as csr.csv in the build directory
- Configure config.h appropriately for your matrix (currently configured for included matrix and test system)
- Use makefile to build using gcc or CMakeLists for CMake+Clang
- Run the executable. Results will be placed in unified.csv and sequential.csv

## Results
Some of the permutations produced using the implementation have been visualized and can be found in the [Graphs](Graphs/) folder

<img height="300" src="https://github.com/CedArctic/parallel-RCM/blob/master/Graphs/offshore.jpg"><img height="300" src="https://github.com/CedArctic/parallel-RCM/blob/master/Graphs/offshore_Reordered.jpg">
