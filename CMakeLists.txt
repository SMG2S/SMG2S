cmake_minimum_required(VERSION 3.6)

#project setting
project(SMG2S)

message([STATUS] "Root File is ${SMG2S_SOURCE_DIR}")

list (APPEND CMAKE_MODULE_PATH "${SMG2S_SOURCE_DIR}/cmake-modules")

# MPI compiler
find_package(MPI REQUIRED)

message([STATUS] "MPICXX Compiler is ${CMAKE_CXX_COMPILER}")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)

if(COMPILER_SUPPORTS_CXX14)
    message([STATUS] "CXX 14 Flags -std=c++14 is found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else()
     message([FATAL_ERROR] "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()

# Add utils and config include
include_directories(${SMG2S_SOURCE_FILE}/utils ${SMG2S_SOURCE_FILE}/config ${SMG2S_SOURCE_FILE}/parVector ${SMG2S_SOURCE_FILE}/parMatrix)

# Generate SMG2S file
add_executable(smg2s.exe smg2s.cpp)
target_link_libraries(smg2s.exe PRIVATE MPI::MPI_CXX)

file(GLOB C_WRAPPERS "interface/C/*.cc")
add_library(smg2s2c STATIC ${C_WRAPPERS})
target_link_libraries(smg2s2c PRIVATE MPI::MPI_CXX)
add_library(smg2s2cs SHARED ${C_WRAPPERS})
target_link_libraries(smg2s2cs PRIVATE MPI::MPI_CXX)


###example
add_subdirectory(example)

#add_executable(c_interface.exe test_c.c)
#target_link_libraries(c_interface.exe PRIVATE MPI::MPI_CXX smg2s2c)
#target_include_directories(c_interface.exe PRIVATE ${CMAKE_SOURCE_DIR})

message([STATUS] "Copy the SMG2S include files: ${INSTALL_TO_USE}")
install(TARGETS smg2s.exe DESTINATION bin)

INSTALL (
    DIRECTORY ${CMAKE_SOURCE_DIR}/utils ${CMAKE_SOURCE_DIR}/smg2s ${CMAKE_SOURCE_DIR}/parMatrix ${CMAKE_SOURCE_DIR}/parVector
    DESTINATION include)

INSTALL(TARGETS smg2s2c DESTINATION lib)
INSTALL(TARGETS smg2s2cs DESTINATION lib)

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

# Find PETSc and SLEPc for verification
option(USE_PETSC "Find PETSc and SLEPc for verification?" OFF)

if(USE_PETSC)
message([STATUS] "Find PETSc and SLEPc for verification: ${USE_PETSC}")
    find_package (PETSc REQUIRED)

    if(PETSC_FOUND)
    message([STATUS] "PETSC is provided, thus one can use the interface of SMG2S to PETSc.")
    INSTALL (
        DIRECTORY ${CMAKE_SOURCE_DIR}/interface/PETSc/
        DESTINATION include/interface/PETSc
        FILES_MATCHING PATTERN "*.h*")
    message([STATUS] "PETSc interface file is copied to include/interface/PETSc.")
    else()
    message([STATUS] "PETSC is not provided, the interface of SMG2S to PETSc can not be used.")
    endif()
endif()


# SMG2S test
enable_testing()
add_test(Test_Size_10000_w_proc1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_CURRENT_BINARY_DIR}/smg2s.exe -D 10000 -L 5 -C 2)
add_test(Test_Size_20000_w_proc2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/smg2s.exe -D 20000 -L 5 -C 2)

add_test(Test_Size_10000_s_proc1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_CURRENT_BINARY_DIR}/smg2s.exe -D 10000 -L 5 -C 2)
add_test(Test_Size_10000_s_proc2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/smg2s.exe -D 10000 -L 5 -C 2)
