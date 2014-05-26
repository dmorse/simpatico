# ==============================================================================
.PHONY: mcMd mcMd-mpi ddMd clean veryclean

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

ddMd:
	./configure -m1
	cd util;  make all
	cd inter; make all
	cd ddMd;  make all

clean:
	cd util;  make clean
	cd inter; make clean
	cd mcMd;  make clean
	cd ddMd;  make clean
	cd mdCf;  make clean

veryclean:
	cd util;  make veryclean
	cd inter; make veryclean
	cd mcMd;  make veryclean
	cd ddMd;  make veryclean
	cd mdCf;  make clean
	rm -f util/config.mk
	rm -f inter/config.mk
	rm -f mcMd/config.mk
	rm -f ddMd/config.mk
	rm -f mdCf/config.mk
	rm -f config.mk

# ==============================================================================
