# ==============================================================================
.PHONY: all clean configure doc clean cleandoc veryclean

SRC_DIR =.

all: 
	cd util; make all
	cd mcMd; make all
	cd main; make all

clean:
	cd util; make clean
	cd mcMd; make clean
	cd main; make clean
 
doc:    
	doxygen Doxyfile

cleandoc:
	-rm -fr ../doc/html/*

veryclean:
	make clean
	make cleandoc
	-rm -f ../lib/*
	-rm -f ../bin/mcSim*
	-rm -f ../bin/mdSim*
	-rm -f mcMd/boundary/Boundary.h
	-rm -f compiler.mk
	-rm -f util/defines.mk
	-rm -f mcMd/defines.mk
	-rm -f main/defines.mk

# ==============================================================================
