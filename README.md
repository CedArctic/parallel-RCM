# parallel-RCM
A parallel Reverse Cuthill-McKee algorithm implementation using OpenMP. Based on [the paper of Rodriguez et al](https://ieeexplore.ieee.org/document/8104594)

## Usage
- Download one of the [releases](https://github.com/CedArctic/parallel-RCM/releases)
- Convert the Sparse Matrix to CSR format
- Place the converted matrix as csr.csv in the build directory
- Configure config.h appropriately for your matrix
- Use makefile to build using gcc or CMakeLists.txt for CMake+Clang
- Run the executable. Results will be placed in classic.csv
