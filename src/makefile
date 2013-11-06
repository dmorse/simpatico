# ==============================================================================
.PHONY: mcMd ddMd mcMd-clean ddMd-clean clean veryclean

include config.mk

mcMd: 
	./configure -m0
	cd util;  make all
	cd inter; make all
	cd mcMd;  make all

mcMd-mpi: 
	./configure -m1
	cd util;  make all
	cd inter; make all
	cd mcMd;  make all

clean-mcMd:
	cd util;  make clean
	cd inter; make clean
	cd mcMd;  make clean

ddMd:
	./configure -m1
	cd util;  make all
	cd inter; make all
	cd ddMd;  make all

clean-ddMd:
	cd util;  make clean
	cd inter; make clean
	cd ddMd;  make clean

clean:
	cd util;  make clean
	cd inter; make clean
	cd mcMd;  make clean
	cd ddMd;  make clean

clean-bin:
	-rm -f $(ROOT_DIR)/bin/mcSim*
	-rm -f $(ROOT_DIR)/bin/mdSim*
	-rm -f $(ROOT_DIR)/bin/ddSim*
 
veryclean:
	cd util;  make veryclean
	cd inter; make veryclean
	cd mcMd;  make veryclean
	cd ddMd;  make veryclean
	make clean-bin
	-rm -f util/config.mk
	-rm -f inter/config.mk
	-rm -f mcMd/config.mk
	-rm -f ddMd/config.mk
	-rm -f config.mk

# ==============================================================================
