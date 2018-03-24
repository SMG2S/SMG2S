cmake_minimum_required(VERSION 2.6)

#project setting
project(smg2s)

# MPI compiler
find_package(MPI REQUIRED)

# Add parVec subdirectory
add_subdirectory(./parVector)

# Add parMat subdirectory
add_subdirectory(./parMatrix)

# Add test subdirectory
add_subdirectory(./tests)

