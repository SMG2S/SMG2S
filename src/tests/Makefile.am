#AM_CCFLAGS = -I$(top_srcdir)/src @AM_CCFLAGS@
###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

# the library names to build (note we are building static libs only)
lib_LIBRARIES = libtestlib.a

# where to install the headers on the system
libtestlib_adir = $(includedir)/testlib

# the list of header files that belong to the library (to be installed later)
libtestlib_a_HEADERS = test_power.hpp 

# the sources to add to the library and to add to the source distribution
libtestlib_a_SOURCES = $(libtestlib_a_HEADERS) \
			test_power.cpp 						
###############################################################################
