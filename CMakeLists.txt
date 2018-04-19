cmake_minimum_required(VERSION 3.6)

#project setting
project(SMG2S)

message([STATUS] "Root File is ${SMG2S_SOURCE_DIR}")

list (APPEND CMAKE_MODULE_PATH "${SMG2S_SOURCE_DIR}/cmake-modules")

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

include(CheckCXXCompilerFlag)  
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)  
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)  
if(COMPILER_SUPPORTS_CXX11)  
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        message([STATUS] "MPICXX Compiler found CXX11 Flags -std=c++11") 
elseif(COMPILER_SUPPORTS_CXX0X)  
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
        message([STATUS] "MPICXX Compiler found CXX11 Flags -std=c++0x") 
else()  
     message([FATAL_ERROR] "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")  
endif() 

# Add parVec subdirectory
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/parVector)

# Add parMat subdirectory
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/parMatrix)

# Add utils and config include 
include_directories(${CMAKE_CURRENT_LIST_DIR}/utils ${CMAKE_CURRENT_LIST_DIR}/config)

# Add src subdirectory
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/matgen)

# Generate SMG2S file
message([STATUS] "We use MPICXX Compiler as ${CMAKE_CXX_COMPILER} to generate
the the SMG2S executable")
add_executable(smg2s.exe main.cpp)

link_directories(${SMG2S_SOURCE_FILE}/${DIR_PARVEC_SRC} ${SMG2S_SOURCE_FILE}/${DIR_PARMAT_SRC} ${SMG2S_SOURCE_FILE}/${DIR_MATGEN_SRC})

target_link_libraries(smg2s.exe PRIVATE ${MPI_CXX_LIBRARIES} parvec parmat matgen)

target_include_directories(smg2s.exe PRIVATE ${MPI_CXX_INCLUDE_PATH})

option(USE_SHARED_LIBRARY "Do we use OpenMP for Compiler?" ON)

if(USE_SHARED_LIBRARY)
message([STATUS] "Generate Shared library smg2s: ${USE_SHARED_LIBRARY}")
add_library(smg2s SHARED ${MPI_CXX_LIBRARIES} ${SMG2S_SOURCE_FILE}/${DIR_PARVEC_SRC} ${SMG2S_SOURCE_FILE}/${DIR_PARMAT_SRC} ${SMG2S_SOURCE_FILE}/${DIR_MATGEN_SRC})
target_include_directories(smg2s PRIVATE ${MPI_CXX_INCLUDE_PATH} ${CMAKE_CURRENT_LIST_DIR}/utils ${CMAKE_CURRENT_LIST_DIR}/config)

SET_TARGET_PROPERTIES(smg2s PROPERTIES LINKER_LANGUAGE CXX)
endif()
find_package (PETSc REQUIRED)

if(PETSC_FOUND)
message([STATUS] "PETSC is provided, thus one can use the interface of SMG2S to PETSc.")
else()
message([STATUS] "PETSC is not provided, the interface of SMG2S to PETSc can not be used.")
endif()

