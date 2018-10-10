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

.PHONY: all doc clean distclean library fast

all: library
	$(MAKE) -C src/*

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
	$(MAKE) -C src/*

distclean: clean
	$(RM) -r doc/*
	$(MAKE) -C include distclean
	$(MAKE) -C src/Voronoi_Network distclean
	$(MAKE) -C src/Single_branch distclean
	$(MAKE) -C src/Bifurcation distclean

fast:
	$(MAKE) -C include fastclean
	$(MAKE) -C include
	$(MAKE) -C src/Voronoi_Network fastclean
	$(MAKE) -C src/Single_branch fastclean
	$(MAKE) -C src/Bifurcation fastclean

