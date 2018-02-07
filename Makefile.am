AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -I.

bin_PROGRAMS = plotfill2
plotfill2_CXXFLAGS = $(libpng_CFLAGS)
plotfill2_LDADD = $(libpng_LIBS)
plotfill2_SOURCES = src/clipper.cpp \
                    src/main.cpp

noinst_HEADERS = src/clipper.hpp
