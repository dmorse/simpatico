# ==============================================================================
.PHONY: mcMd mcMd-mpi ddMd clean veryclean

# Serial versions of mdSim and mcSim MD and MC programs (MPI disabled)
mcMd: 
	./configure -m0
	cd util; make all
	cd inter; make all
	cd mcMd; make all

# Embarassingly parallel mdSim and mcSim MD and MC programs (MPI enabled)
mcMd-mpi: 
	./configure -m1
	cd util; make all
	cd inter; make all
	cd mcMd; make all
 
# Domain-decomposition ddSim molecular dynamics program (MPI enabled)
ddMd:
	./configure -m1
	cd util; make all
	cd inter; make all
	cd ddMd; make all

# Remove object (*.o), dependency (*.d) and library (*.a) files
clean:
	cd util; make clean
	cd inter; make clean
	cd mcMd; make clean
	cd ddMd; make clean

# Remove all automatically generated files, return to as-distributed state
veryclean:
	cd util; make veryclean
	cd inter; make veryclean
	cd mcMd; make veryclean
	cd ddMd; make veryclean
	rm -f util/config.mk
	rm -f inter/config.mk
	rm -f mcMd/config.mk
	rm -f ddMd/config.mk
	rm -f config.mk

# ==============================================================================
