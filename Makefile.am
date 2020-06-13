AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -I. -Wall -Wextra -std=cpp11 -O2

bin_PROGRAMS = plotfill2
plotfill2_LDADD = $(LIBPNG)
plotfill2_SOURCES = src/clipper.cpp \
                    src/main.cpp

noinst_HEADERS = src/clipper.hpp
