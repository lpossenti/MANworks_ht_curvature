# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2014-2015
#
#                    Copyright D. Notaro 2015
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for test simulations
#   AUTHOR      : Domenico Notaro <domenico.not@gmail.com>
#   DATE        : April 2015
# ==================================================================== 
include ../config.mk

CPPFLAGS=-I. -I$(mkGetfemInc) -I$(mkBoostInc)
CXXFLAGS+=-std=c++14
ifeq ($(WITH_SAMG), 1)
CXXFLAGS+=-DWITH_SAMG 
endif
ifeq ($(DEBUG),yes)
  OPTFLAGS=-g -Wall
else
  OPTFLAGS=-O3 -march=native
  CPPFLAGS+=-DNDEBUG
endif
LIBDIR=../lib
LIBNAME=problem3d1d
LIBFILE=lib$(LIBNAME).a

LIBSRC=$(wildcard *.cpp)
LIBOBJS=$(LIBSRC:.cpp=.o)
LIBHEADERS=$(wildcard *.hpp)

.PHONY: all clean distclean library fastclean

all: library
	@echo
	@echo Library installed!

library: $(LIBOBJS)
	install -d $(LIBDIR)
	ar -r $(LIBDIR)/$(LIBFILE) $(LIBOBJS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -o $@ -c $<
 
problem3d1d.cpp: darcy_preconditioner_mon.hpp darcy_preconditioner.hpp darcy_preconditioner_mon_coup.hpp darcy_preconditioner_tissue_coup.hpp
	@touch $@

clean:
	$(RM) $(LIBOBJS) *~

distclean: clean
	$(RM) $(LIBDIR)/$(LIBFILE)

fastclean:
	$(RM) problemHT.o
