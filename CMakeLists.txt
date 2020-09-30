cmake_minimum_required(VERSION 3.15)
project(RCM C)

set(CMAKE_C_STANDARD 99)

OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
    FIND_PACKAGE(OpenMP)
    IF(OPENMP_FOUND)
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ENDIF()
ENDIF()

# Set debugging flags to use gprof
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

# OpenMP Environment Variables Optimizations
SET(OMP_WAIT_POLICY "passive")
#SET(OMP_DYNAMIC "false")

add_executable(RCM main_unified.c src/rcm/rcm.h src/rcm/rcm.c src/quicksort/quicksort.h src/quicksort/quicksort.c
        src/rcm_parallel/rcm_parallel.c src/rcm_parallel/rcm_parallel.h src/csr/csr.c src/csr/csr.h src/graph/graph.c
        src/graph/graph.h config.h)

target_link_libraries(RCM m)