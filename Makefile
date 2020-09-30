# The ARM cross-compiler: Using abhiTronix's gcc optimized for Raspberry Pi's
# https://github.com/abhiTronix/raspberry-pi-cross-compilers/wiki/Cross-Compiler:-Installation-Instructions
CC = gcc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#  -lm   links math library
#  -pthread   links the pthreads library
CFLAGS  = -g -Wall -O3
MATH = -lm
OPENMP = -fopenmp
RM = rm

# the build target executable:
TARGET = RCM_GCC

all:	$(TARGET)

$(TARGET):	main_unified.c graph.o rcm.o rcm_parallel.o csr.o quicksort.o
	$(CC) $(CFLAGS) -o $(TARGET) $(MATH) $(OPENMP) main_unified.c graph.o rcm.o rcm_parallel.o csr.o quicksort.o

graph.o: src/graph/graph.h src/graph/graph.c
	$(CC) $(CFLAGS) $(MATH) $(OPENMP) -c src/graph/graph.c

rcm.o:  quicksort.o graph.o src/rcm/rcm.h src/rcm/rcm.c
	$(CC) $(CFLAGS) -c src/rcm/rcm.c

quicksort.o:    src/quicksort/quicksort.h src/quicksort/quicksort.c
	$(CC) $(CFLAGS) -c src/quicksort/quicksort.c

rcm_parallel.o:    rcm.o graph.o src/rcm_parallel/rcm_parallel.h src/rcm_parallel/rcm_parallel.c
	$(CC) $(CFLAGS) $(MATH) $(OPENMP) -c src/rcm_parallel/rcm_parallel.c

csr.o:    src/csr/csr.h src/csr/csr.c
	$(CC) $(CFLAGS) -c src/csr/csr.c

clean:
	$(RM) $(TARGET) graph.o rcm.o quicksort.o rcm_parallel.o csr.o

clear:
	$(RM) $(TARGET) graph.o rcm.o quicksort.o rcm_parallel.o csr.o sequential.csv classic.csv
