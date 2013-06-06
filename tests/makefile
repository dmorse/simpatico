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

mcMd_mpi:
	cd ..; ./configure -m1
	cd mcMd/param/mpi; make all; mpirun -np 2 ./Test > log; grep failed log
	cd mcMd/param/mpi; make all; mpirun -np 2 ./MpiTest > log; grep failed log

ddMd:
	cd ..; ./configure -m1
	cd ddMd; make all; mpirun -np 6 ./Test > log; grep failed log

clean: 
	cd util;  make clean
	cd inter; make clean
	cd mcMd;  make clean
	cd ddMd;  make clean

veryclean: 
	make clean;
	cd test; make clean

# ==============================================================================
