# ==============================================================================
.PHONY: mcMd mcMd-mpi ddMd clean veryclean

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

veryclean:
	make clean
	-rm -f util/config.mk
	-rm -f inter/config.mk
	-rm -f mcMd/config.mk
	-rm -f ddMd/config.mk
	-rm -f config.mk
	-rm -f */*.a

# ==============================================================================
