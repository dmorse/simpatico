#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Include makefile fragments

BLD_DIR_REL =..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

#-----------------------------------------------------------------------
# Makefile variable definitions

# Location of *.o and *.d files for main programs
mdSim=$(BLD_DIR)/mcMd/mdSim
mcSim=$(BLD_DIR)/mcMd/mcSim

#-----------------------------------------------------------------------
# Main targets

# Note: Library and exectable paths are defined in mcMd/config.mk

all: $(mcMd_LIB) $(mcSim_EXE) $(mdSim_EXE)

lib: (mcMd_LIB)

clean:
	rm -f $(mcMd_OBJS) $(mcMd_OBJS:.o=.d) $(mcMd_LIB)
	rm -f $(mdSim).o $(mdSim).d
	rm -f $(mcSim).o $(mcSim).d
	cd tests; $(MAKE) clean
	rm -f *.o */*.o */*/*.o */*/*/*.o */*/*/*/*.o */*/*/*/*/*.o
	rm -f *.d */*.d */*/*.d */*/*/*.d */*/*/*/*.d */*/*/*/*/*.d

veryclean:
	$(MAKE) clean
	rm -f lib*.a
ifeq ($(BLD_DIR),$(SRC_DIR))
	rm -f configIos/ConfigIoFactory.cpp
	rm -f analyzers/mcSystem/McAnalyzerFactory.cpp
	rm -f analyzers/mdSystem/MdAnalyzerFactory.cpp
	rm -f analyzers/system/SystemAnalyzerFactory.cpp
	rm -f commands/McCommandFactory.cpp
	rm -f commands/MdCommandFactory.cpp
	rm -f mcMoves/McMoveFactory.cpp
	rm -f mdIntegrators/MdIntegratorFactory.cpp
	rm -f perturb/mcSystem/McPerturbationFactory.cpp
	rm -f potentials/angle/AngleFactory.cpp
	rm -f potentials/bond/BondFactory.cpp
	rm -f potentials/coulomb/CoulombFactory.cpp
	rm -f potentials/dihedral/DihedralFactory.cpp
	rm -f potentials/external/ExternalFactory.cpp
	rm -f potentials/link/LinkFactory.cpp
	rm -f potentials/pair/PairFactory.cpp
	rm -f potentials/special/SpecialFactory.cpp
	rm -f species/SpeciesFactory.cpp
	rm -f trajectory/TrajectoryReaderFactory.cpp
	rm -f generators/generatorFactory.cpp
	rm -f mdSim.cpp
	rm -f mcSim.cpp
endif

# Executable targets
# Note: mcSim_EXE and mdSim_EXE are defined in mcMd/config.mk

$(mcSim_EXE): $(mcSim).o $(LIBS)
	$(CXX) -o $(mcSim_EXE) $(mcSim).o $(LIBS) $(LDFLAGS)
 
$(mdSim_EXE): $(mdSim).o $(LIBS)
	$(CXX) -o $(mdSim_EXE) $(mdSim).o $(LIBS) $(LDFLAGS)

# Short names for executable targets (for convenience)
 
mcSim:
	$(MAKE) $(mcSim_EXE)

mdSim:
	$(MAKE) $(mdSim_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(util_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(mdSim).d 
-include $(mcSim).d 
