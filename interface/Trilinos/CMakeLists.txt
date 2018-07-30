cmake_minimum_required(VERSION 2.6)

MESSAGE("Starting configure")

#set TRILINOS_PATH of the platform
#set(TRILINOS_PATH /Users/xwu/phd_xwu/trilinos-12.12.1-Source/installation)
set(TRILINOS_PATH /Users/xwu/phd_xwu/trilinos-12.12.1-Source2/installation2)
#set(TRILINOS_PATH /home/xinzhewu/workstations/trilinos/install)

SET(CMAKE_PREFIX_PATH ${TRILINOS_PATH} ${CMAKE_PREFIX_PATH})

#find trilinos
FIND_PACKAGE(Trilinos REQUIRED PATHS ${TRILINOS_PATH}/include ${TRILINOS_PATH})

#Output some info about the Trilinos installation found
IF(Trilinos_FOUND)
   MESSAGE("\nFound Trilinos!  Here are the details: ")
   MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
   MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
   MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
   MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
   MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
   MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
   MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
   MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
   # Will Tramonto work with a Shared version of Trilinos?
   MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
   #These last outputs are temporary
   MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
   MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
   MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
   MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS =
     ${Trilinos_Fortran_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
   MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
   MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
   # Add MPI_EXEC, MPI_EXEC_MAX_NUMPROCS, MPI_EXEC_NUMPROC_FLAGS
   # These are not available in Trilinos 10.2.0, but should be in 10.4.0
   MESSAGE("End of Trilinos details\n")
ELSE()
  # We should never reach this because FIND_PACKAGE(Trilinos now uses REQUIRED
  MESSAGE(FATAL_ERROR "Could not find Trilinos 10.2 or greater!")
ENDIF()

IF(NOT CMAKE_CXX_COMPILER)
  MESSAGE("Setting CXX compiler.")
  SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
ENDIF()

IF(NOT CMAKE_C_COMPILER)
  MESSAGE("Setting C compiler.")
  SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
ENDIF()

IF(NOT CMAKE_LINKER)
  MESSAGE("Setting linker.")
  SET(CMAKE_LINKER ${Trilinos_LINKER})
ENDIF()

IF(NOT CMAKE_EXE_LINKER_FLAGS)
  MESSAGE("Setting linker flags.")
  SET(CMAKE_EXE_LINKER_FLAGS ${Trilinos_EXTRA_LD_FLAGS})
ENDIF()

IF(NOT CMAKE_AR)
  MESSAGE("Setting ar.")
  SET(CMAKE_AR ${CMAKE_AR})
ENDIF()

project(trilinos_tutorial)
SET(CPACK_PACKAGE_VERSION "1.0")
set(TRILINOS_PATH "" CACHE FILEPATH "Path to Trilinos Directory")

IF(NOT CMAKE_CXX_FLAGS)
  MESSAGE("Setting CXX compiler flags.")
  SET(CMAKE_CXX_FLAGS ${Trilinos_CXX_COMPILER_FLAGS})
ENDIF()

IF(NOT CMAKE_C_FLAGS)
  MESSAGE("Setting C compiler flags.")
  SET(CMAKE_C_FLAGS ${Trilinos_C_COMPILER_FLAGS})
ENDIF()

ENABLE_TESTING()

#ADD_SUBDIRECTORY(Solvers)
#ADD_SUBDIRECTORY(Solvers/GMRES)

#Add Trilinos information to the include and link lines
include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
link_directories(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )
# /Library/Frameworks/QtCore.framework /Library/Frameworks/QtGui.framework)

#set trilinos libraries to link (LINK_LIBRARIES)
set(LINK_LIBRARIES ${Teuchos_LIBRARIES} ${Tpetra_LIBRARIES})

#set Anasazi libraries to link (LINK_LIBRARIES_ARNOLDI)
set(LINK_LIBRARIES_ARNOLDI ${Tpetra_LIBRARIES} ${Anasazi_LIBRARIES} ${Amesos_LIBRARIES})

#set Belos libraries to link (LINK_LIBRARIES_GMRES)
set(LINK_LIBRARIES_GMRES ${Tpetra_LIBRARIES} ${Teuchos_LIBRARIES} ${Belos_LIBRARIES})

#add executable
add_executable(test.exe main.cpp)
target_link_libraries(test.exe ${LINK_LIBRARIES})

include(Dart)
include(CPack)
