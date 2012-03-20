# ==============================================================================
.PHONY: all clean configure doc clean cleandoc veryclean

include compiler.mk

all: 
	cd util;  make all
	cd inter; make all
	cd mcMd;  make all
ifdef UTIL_MPI
	cd ddMd;  make all
endif

clean:
	cd util;  make clean
	cd inter; make clean
	cd mcMd;  make clean
ifdef UTIL_MPI
	cd ddMd;  make clean
endif
 
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
	-rm -f ../bin/ddSim*
	-rm -f util/boundary/Boundary.h
	-rm -f compiler.mk
	-rm -f util/defines.mk
	-rm -f inter/defines.mk
	-rm -f mcMd/defines.mk
	-rm -f ddMd/defines.mk

# ==============================================================================
