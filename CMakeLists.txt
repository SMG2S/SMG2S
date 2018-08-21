cmake_minimum_required(VERSION 3.6)

#project setting
project(SMG2S)

message([STATUS] "Root File is ${SMG2S_SOURCE_DIR}")

list (APPEND CMAKE_MODULE_PATH "${SMG2S_SOURCE_DIR}/cmake-modules")

# MPI compiler
find_package(MPI REQUIRED)

if(NOT MPI_CXX_FOUND)
	message([FATAL_ERROR] "MPI CXX compiler was not found and is required ")
endif()

# Set MPICXX as the CXX compiler

SET(CMAKE_CXX_COMPILER mpicxx)

message([STATUS] "MPICXX Compiler is ${CMAKE_CXX_COMPILER}")

option(USE_OPENMP "Do we use OpenMP for Compiler?" OFF)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

if(USE_OPENMP AND COMPILER_SUPPORTS_CXX11)
    message([STATUS] "USE OpenMP: ${USE_OPENMP} and CXX 11 Flags -std=c++11")
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

elseif(USE_OPENMP AND COMPILER_SUPPORTS_CXX0X)
    message([STATUS] "USE OpenMP: ${USE_OPENMP} and CXX 11 Flags -std=c++0x")
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

elseif(COMPILER_SUPPORTS_CXX11)
    message([STATUS] "DO NOT use OpenMP, and CXX 11 Flags -std=c++11 is found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

elseif(COMPILER_SUPPORTS_CXX0X)
    message([STATUS] "DO NOT use OpenMP, and CXX 11 Flags -std=c++0x is found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

else()
     message([FATAL_ERROR] "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

configure_file (
  "${PROJECT_SOURCE_DIR}/config/config.h.in"
  "${PROJECT_SOURCE_DIR}/config/config.h"
  )

option (__USE_COMPLEX__ "Use the complex scalar type" ON)

if(__USE_COMPLEX__)
	message([STATUS] "USE the complex scalar to compile the SMG2S")
endif()

option (__USE_DOUBLE__ "Use the double precision" ON)

if(__USE_DOUBLE__)
	message([STATUS] "USE the double precision to compile the SMG2S")
endif()

option (__USE_64BIT__ "Use the long integer" OFF)

if(__USE_64BIT__)
	message([STATUS] "USE the long integer to compile the SMG2S")
endif()


# Add utils and config include
include_directories(${SMG2S_SOURCE_FILE}/utils ${SMG2S_SOURCE_FILE}/config ${SMG2S_SOURCE_FILE}/parVector ${SMG2S_SOURCE_FILE}/parMatrix)

# Generate SMG2S file
message([STATUS] "We use MPICXX Compiler as ${CMAKE_CXX_COMPILER} to generate
the the SMG2S executable")

add_executable(smg2s.exe main.cpp)

target_link_libraries(smg2s.exe PRIVATE ${MPI_CXX_LIBRARIES})

target_include_directories(smg2s.exe PRIVATE ${MPI_CXX_INCLUDE_PATH})

file(GLOB C_WRAPPERS "interface/C/*.cc")
add_library(smg2s2c SHARED ${C_WRAPPERS})

option(INSTALL_TO_USE "Install SMG2S include files?" ON)

if(INSTALL_TO_USE)
message([STATUS] "Copy the SMG2S include files: ${INSTALL_TO_USE}")
install(TARGETS smg2s.exe DESTINATION bin)
INSTALL (
    DIRECTORY ${CMAKE_SOURCE_DIR}/utils ${CMAKE_SOURCE_DIR}/smg2s ${CMAKE_SOURCE_DIR}/parMatrix ${CMAKE_SOURCE_DIR}/parVector
    DESTINATION include)
INSTALL(TARGETS smg2s2c DESTINATION lib)
INSTALL (
    DIRECTORY ${CMAKE_SOURCE_DIR}/interface/C/
    DESTINATION include/interface/C
    FILES_MATCHING PATTERN "*.h*")
INSTALL (
		DIRECTORY ${CMAKE_SOURCE_DIR}/config/
		DESTINATION include/config
		FILES_MATCHING PATTERN "*.h")
INSTALL (
        DIRECTORY ${CMAKE_SOURCE_DIR}/interface/Trilinos/
        DESTINATION include/interface/Trilinos
        FILES_MATCHING PATTERN "trilinos_interface.hpp")
endif()

# Find PETSc and SLEPc for verification
option(USE_PETSC "Find PETSc and SLEPc for verification?" ON)

if(USE_PETSC)
message([STATUS] "Find PETSc and SLEPc for verification: ${USE_PETSC}")
    find_package (PETSc REQUIRED)

    if(PETSC_FOUND)
    message([STATUS] "PETSC is provided, thus one can use the interface of SMG2S to PETSc.")
    if(INSTALL_TO_USE)
    INSTALL (
        DIRECTORY ${CMAKE_SOURCE_DIR}/interface/PETSc/
        DESTINATION include/interface/PETSc
        FILES_MATCHING PATTERN "*.h*")
    message([STATUS] "PETSc interface file is copied to include/interface/PETSc.")
    endif()
    else()
    message([STATUS] "PETSC is not provided, the interface of SMG2S to PETSc can not be used.")
    endif()
endif()
