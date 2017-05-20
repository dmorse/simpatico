# ==============================================================================
.PHONY: mcMd mcMd-mpi ddMd tools clean veryclean

# Serial versions of mdSim and mcSim MD and MC programs (MPI disabled)
mcMd: 
	./configure -m0
	cd util; $(MAKE) all
	cd simp; $(MAKE) all
	cd mcMd; $(MAKE) all

# Embarassingly parallel mdSim and mcSim MD and MC programs (MPI enabled)
mcMd-mpi: 
	./configure -m1
	cd util; $(MAKE) all
	cd simp; $(MAKE) all
	cd mcMd; $(MAKE) all
 
# Domain-decomposition ddSim molecular dynamics program (MPI enabled)
ddMd:
	./configure -m1
	cd util; $(MAKE) all
	cd simp; $(MAKE) all
	cd ddMd; $(MAKE) all

# Single-processor analysis program (MPI disabled)
tools:
	./configure -m0
	cd util; $(MAKE) all
	cd simp; $(MAKE) all
	cd tools; $(MAKE) all

# Remove object (*.o), dependency (*.d) and library (*.a) files
clean:
	cd util; $(MAKE) clean
	cd simp; $(MAKE) clean
	cd mcMd; $(MAKE) clean
	cd tools; $(MAKE) clean
	cd ddMd; $(MAKE) clean

# Clean unit tests
clean-tests:
	cd util/tests; $(MAKE) clean
	cd simp/tests; $(MAKE) clean
	cd mcMd/tests; $(MAKE) clean
	cd tools/tests; $(MAKE) clean
	cd ddMd/tests; $(MAKE) clean

# Remove all automatically generated files, return to as-distributed state
veryclean:
	cd util; $(MAKE) veryclean
	cd simp; $(MAKE) veryclean
	cd mcMd; $(MAKE) veryclean
	cd ddMd; $(MAKE) veryclean
	cd tools; $(MAKE) veryclean
ifeq ($(BLD_DIR),$(SRC_DIR))
	rm -f simp/user/sources.mk
	rm -f mcMd/user/sources.mk
	rm -f ddMd/user/sources.mk
	rm -f tools/user/sources.mk
	rm -f mcMd/modules/sources.mk
	rm -f mcMd/mcSim.cpp
	rm -f mcMd/mdSim.cpp
endif
	rm -f util/config.mk
	rm -f simp/config.mk
	rm -f mcMd/config.mk
	rm -f ddMd/config.mk
	rm -f tools/config.mk
	rm -f config.mk

# ==============================================================================
