cmake_minimum_required(VERSION 2.6)

#project setting
project(SMG2S)

message([STATUS] "Root File is ${SMG2S_SOURCE_DIR}")

# MPI compiler
find_package(MPI REQUIRED)

if(NOT MPI_CXX_FOUND)
	message([FATAL_ERROR] "MPI CXX compiler was not found and is required ")
endif()

# Set MPICXX as the CXX compiler

SET(CMAKE_CXX_COMPILER mpicxx)

# Add parVec subdirectory
add_subdirectory(parVector)

# Add parMat subdirectory
add_subdirectory(parMatrix)

# Add test subdirectory
add_subdirectory(tests)

