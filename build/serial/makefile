# ==============================================================================
.PHONY: mcMd ddMd mcMd-clean ddMd-clean clean veryclean

include compiler.mk

mcMd: 
	cd $(ROOT_DIR); ./configure -m0
	cd util;  make all
	cd inter; make all
	cd mcMd;  make all

mcMd-mpi: 
	cd $(ROOT_DIR); ./configure -m1
	cd util;  make all
	cd inter; make all
	cd mcMd;  make all
ddMd:
	cd $(ROOT_DIR); ./configure -m1
	cd util;  make all
	cd inter; make all
	cd ddMd;  make all

clean-mcMd:
	cd util;  make clean
	cd inter; make clean
	cd mcMd;  make clean

clean-ddMd:
	cd util;  make clean
	cd inter; make clean
	cd mcMd;  make clean

clean:
	cd util;  make clean
	cd inter; make clean
	cd mcMd;  make clean
	cd ddMd;  make clean

clean-bin:
	-rm -f $(ROOT_DIR)/lib/*
	-rm -f $(ROOT_DIR)/bin/mcSim*
	-rm -f $(ROOT_DIR)/bin/mdSim*
	-rm -f $(ROOT_DIR)/bin/ddSim*
 
veryclean:
	make clean
	make clean-bin
	-rm -f util/boundary/Boundary.h
	-rm -f mcMd/configIos/ConfigIoFactory.cpp
	-rm -f mcMd/diagnostics/mcSystem/McDiagnosticFactory.cpp
	-rm -f mcMd/diagnostics/mdSystem/MdDiagnosticFactory.cpp
	-rm -f mcMd/diagnostics/system/SystemDiagnosticFactory.cpp
	-rm -f mcMd/mcMoves/McMoveFactory.cpp
	-rm -f mcMd/mdIntegrators/MdIntegratorFactory.cpp
	-rm -f mcMd/perturb/mcSystem/McPerturbationFactory.cpp
	-rm -f mcMd/potentials/angle/AngleFactory.cpp
	-rm -f mcMd/potentials/bond/BondFactory.cpp
	-rm -f mcMd/potentials/dihedral/DihedralFactory.cpp
	-rm -f mcMd/potentials/external/ExternalFactory.cpp
	-rm -f mcMd/potentials/link/LinkFactory.cpp
	-rm -f mcMd/potentials/pair/PairFactory.cpp
	-rm -f mcMd/species/SpeciesFactory.cpp
	-rm -f mcMd/trajectoryIos/TrajectoryIoFactory.cpp
	-rm -f ddMd/configIos/ConfigIoFactory.cpp
	-rm -f ddMd/diagnostics/DiagnosticFactory.cpp
	-rm -f ddMd/integrators/IntegratorFactory.cpp
	-rm -f ddMd/potentials/angle/AngleFactory.cpp
	-rm -f ddMd/potentials/bond/BondFactory.cpp
	-rm -f ddMd/potentials/dihedral/DihedralFactory.cpp
	-rm -f ddMd/potentials/external/ExternalFactory.cpp
	-rm -f ddMd/potentials/pair/PairFactory.cpp
	cd ../tests; make clean
	-rm -f util/defines.mk
	-rm -f inter/defines.mk
	-rm -f mcMd/defines.mk
	-rm -f ddMd/defines.mk
	-rm -f compiler.mk

# ==============================================================================
