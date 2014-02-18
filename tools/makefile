###################################################
#
# Makefile for tools
#
###################################################

#
# Macros
#

ZDIR = ../external/zlib-1.2.5
ZLIB = libz.a

BAMTOOLSDIR = ../external/BamTools-2.0/

INCLUDE = -I. -I../external -I$(ZDIR) -I$(BAMTOOLSDIR) -I$(HOME)/Install/include/coin
CC = g++
CXX = g++
CXXFLAGS = -O3 -ggdb -DDEBUG_CHECKS $(INCLUDE)
LDFLAGS = -O3 -ggdb -L$(ZDIR) -lz -lboost_serialization -lboost_iostreams -lipopt -lblas -llapack -lgfortran -lcoinmumps -lcoinmetis

PyBreakCopies : \
		Common.o\
		DebugCheck.o\
		Parsers.o\
		PyBreakCopies.o


