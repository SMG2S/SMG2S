# - Find PAPI
# Find the PAPI libraries
# PAPI: Performance Application Programming Interface
# Website: https://icl.utk.edu/papi/index.html
#   PAPI_FOUND            : True if PAPI_INCUDE_DIR are found
#   PAPI_INCLUDE_DIR      : where to find papi.h, etc.
#   PAPI_INCLUDE_DIRS     : set when PAPU_INCLUDE_DIR found
#   PAPI_LIBRARIES        : the library to link against.

set(PAPI_ROOT "/usr/local" CACHE PATH "Folder containing PAPI libraries")

if (NOT PAPI_ROOT AND DEFINED ENV{PAPI_ROOT})
  set(PAPI_ROOT $ENV{PAPI_ROOT} CACHE PATH "Folder containing PAPI")
elseif (NOT PAPI_ROOT AND DEFINED ENV{PAPIROOT})
  set(BLIS_ROOT $ENV{PAPIROOT} CACHE PATH "Folder containing PAPI")
endif()

find_path(PAPI_INCLUDE_DIR
    NAMES
        papi.h
    PATHS
        ${BLIS_ROOT}/include
)

find_library(PAPI_LIBRARY libpapi.so libpapi.a papi
   PATHS ${PAPI_ROOT}/lib
)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(BLIS
        DEFAULT_MSG
        PAPI_LIBRARY
        PAPI_INCLUDE_DIR
)

if(PAPI_FOUND)
    set(PAPI_INCLUDE_DIRS ${PAPI_INCLUDE_DIR})
    set(PAPI_LIBRARIES ${PAPI_LIBRARY})
endif()


