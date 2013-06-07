include ../src/compiler.mk
# ==============================================================================
.PHONY: clean veryclean util inter mcMd mcMd_mpi

util:
	cd ..; ./configure -m0
	cd util; make all; ./Test > log; grep failed log

inter:
	cd ..; ./configure -m0
	cd util; make all; ./Test > log; grep failed log

mcMd:
	cd ..; ./configure -m0
	cd mcMd; make all; ./Test > log; grep failed log

mcMd-mpi:
	cd ..; ./configure -m1
	cd util/mpi; make all; $(MPIRUN) 2 ./Test > log; grep failed log
	cd util/param/mpi; make all; $(MPIRUN) 2 ./MpiTest > log; grep failed log

ddMd:
	cd ..; ./configure -m1
	cd ddMd; make all; $(MPIRUN) 6 ./Test > log; grep failed log

clean: 
	cd util;  make clean
	cd inter; make clean
	cd mcMd;  make clean
	cd ddMd;  make clean
	-rm -f log

veryclean: 
	make clean;
	cd test; make clean

# ==============================================================================
