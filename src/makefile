# ==============================================================================
.PHONY: mcMd mcMd-mpi ddMd spAn clean veryclean

# Serial versions of mdSim and mcSim MD and MC programs (MPI disabled)
mcMd: 
	./configure -m0
	cd util; $(MAKE) all
	cd inter; $(MAKE) all
	cd mcMd; $(MAKE) all

# Embarassingly parallel mdSim and mcSim MD and MC programs (MPI enabled)
mcMd-mpi: 
	./configure -m1
	cd util; $(MAKE) all
	cd inter; $(MAKE) all
	cd mcMd; $(MAKE) all
 
# Domain-decomposition ddSim molecular dynamics program (MPI enabled)
ddMd:
	./configure -m1
	cd util; $(MAKE) all
	cd inter; $(MAKE) all
	cd ddMd; $(MAKE) all

# Single-processor analysis program (MPI disabled)
spAn:
	./configure -m0
	cd util; $(MAKE) all
	cd inter; $(MAKE) all
	cd spAn; $(MAKE) all

# Remove object (*.o), dependency (*.d) and library (*.a) files
clean:
	cd util; $(MAKE) clean
	cd inter; $(MAKE) clean
	cd mcMd; $(MAKE) clean
	cd ddMd; $(MAKE) clean
	cd spAn; $(MAKE) clean

# Remove all automatically generated files, return to as-distributed state
veryclean:
	cd util; $(MAKE) veryclean
	cd inter; $(MAKE) veryclean
	cd mcMd; $(MAKE) veryclean
	cd ddMd; $(MAKE) veryclean
	cd spAn; $(MAKE) veryclean
	rm -f util/config.mk
	rm -f inter/config.mk
	rm -f mcMd/config.mk
	rm -f ddMd/config.mk
	rm -f spAn/config.mk
	rm -f config.mk

# ==============================================================================
