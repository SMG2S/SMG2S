# additional include pathes necessary to compile the C++ library
#AM_CCFLAGS = -I$(top_srcdir)/src @AM_CCFLAGS@

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = testpower

# list of sources for the 'helloWorld' binary
testpower_SOURCES = test.cpp

# the additional libraries needed to link helloWorld
testpower_LDADD = $(top_builddir)/src/tests/libtestlib.a $(AM_LDFLAGS)

###############################################################################
