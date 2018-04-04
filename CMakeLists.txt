cmake_minimum_required(VERSION 3.6)

#project setting
project(SMG2S)

message([STATUS] "Root File is ${SMG2S_SOURCE_DIR}")

#SET(CMAKE_C_FLAGS /apps/intel/2017.1/bin/icc)
#SET(CMAKE_CXX_FLAGS /apps/intel/2017.1/bin/icpc)

# MPI compiler
find_package(MPI REQUIRED)

if(NOT MPI_CXX_FOUND)
	message([FATAL_ERROR] "MPI CXX compiler was not found and is required ")
endif()

# Set MPICXX as the CXX compiler

SET(CMAKE_CXX_COMPILER mpicxx)

message([STATUS] "MPICXX Compiler is ${CMAKE_CXX_COMPILER}")

option(USE_OPENMP "Do we use OpenMP for Compiler?" OFF)

if(USE_OPENMP)
message([STATUS] "USE OpenMP: ${USE_OPENMP}")
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

else()
message([STATUS] "Don't USE OpenMP: ${USE_OPENMP}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
endif()

# Add parVec subdirectory
add_subdirectory(parVector)

# Add parMat subdirectory
add_subdirectory(parMatrix)

# Add test subdirectory
add_subdirectory(tests)

