# ==============================================================================
.PHONY: mcMd ddMd clean veryclean

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
	-rm -f util/defines.mk
	-rm -f inter/defines.mk
	-rm -f mcMd/defines.mk
	-rm -f ddMd/defines.mk
	-rm -f compiler.mk

# ==============================================================================
