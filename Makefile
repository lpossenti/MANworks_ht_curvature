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
#   DATE        : November 2015
# ====================================================================

.PHONY: all doc clean distclean library

all: library
	$(MAKE) -C src/1_uncoupled
	$(MAKE) -C src/2_singlebranch
	$(MAKE) -C src/3_bifurcation
	$(MAKE) -C src/4_anastomosis
	$(MAKE) -C src/5_rhombus
	$(MAKE) -C src/6_splitted_singlebranch

library: 
	$(MAKE) -C include

doc:
	install -d doc
	doxygen Doxyfile
	
pdf: doc
	$(MAKE) -C doc/latex pdf

clean:
	$(RM) -r *~ *.log
	$(MAKE) -C include clean
	$(MAKE) -C src/1_uncoupled clean
	$(MAKE) -C src/2_singlebranch clean
	$(MAKE) -C src/3_bifurcation clean
	$(MAKE) -C src/4_anastomosis clean
	$(MAKE) -C src/5_rhombus clean
	$(MAKE) -C src/6_splitted_singlebranch clean

distclean: clean
	$(RM) -r doc/*
	$(MAKE) -C include distclean
	$(MAKE) -C src/1_uncoupled distclean
	$(MAKE) -C src/2_singlebranch distclean
	$(MAKE) -C src/3_bifurcation distclean
	$(MAKE) -C src/4_anastomosis distclean
	$(MAKE) -C src/5_rhombus distclean
	$(MAKE) -C src/6_splitted_singlebranch distclean
