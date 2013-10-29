# ==============================================================================
.PHONY: mcMd ddMd mcMd-clean ddMd-clean clean veryclean

include compiler.mk

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
	-rm -f util/defines.mk
	-rm -f inter/defines.mk
	-rm -f mcMd/defines.mk
	-rm -f ddMd/defines.mk
	-rm -f compiler.mk

# ==============================================================================
