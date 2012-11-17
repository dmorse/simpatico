# ==============================================================================
.PHONY: mcMd ddMd mcMd-clean ddMd-clean clean veryclean

include compiler.mk

mcMd: 
	cd util;  make all
	cd inter; make all
	cd mcMd;  make all

ddMd:
	cd util;  make all
	cd inter; make all
	cd ddMd;  make all

clean-mcMd:
	cd util;  make clean
	cd inter; make clean
	cd mcMd;  make clean

clean-ddMd:
	cd util;  make clean
	cd inter; make clean
	cd mcMd;  make clean

clean:
	cd util;  make clean
	cd inter; make clean
	cd mcMd;  make clean
	cd ddMd;  make clean
 
veryclean:
	make clean
	cd ../tests; make clean
	-rm -f ../lib/*
	-rm -f ../bin/mcSim*
	-rm -f ../bin/mdSim*
	-rm -f ../bin/ddSim*
	-rm -f util/boundary/Boundary.h
	-rm -f util/defines.mk
	-rm -f inter/defines.mk
	-rm -f mcMd/defines.mk
	-rm -f ddMd/defines.mk

# ==============================================================================
