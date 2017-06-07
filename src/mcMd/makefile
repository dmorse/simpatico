#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Makefile includes

BLD_DIR_REL =..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/mcMd/sources.mk

#-----------------------------------------------------------------------
# Variable definitions

MDSIM=$(BLD_DIR)/mcMd/mdSim
MCSIM=$(BLD_DIR)/mcMd/mcSim
#-----------------------------------------------------------------------
# Main targets

all: $(mcMd_OBJS) $(mcMd_LIB) $(MDSIM).o $(MCSIM).o $(MCSIM_EXE) $(MDSIM_EXE)

lib: $(mcMd_OBJS) (mcMd_LIB)

clean:
	rm -f $(mcMd_OBJS) $(mcMd_OBJS:.o=.d) $(mcMd_LIB)
	rm -f $(MDSIM).o $(MDSIM).d
	rm -f $(MCSIM).o $(MCSIM).d
	cd tests; $(MAKE) clean
	rm -f *.o */*.o */*/*.o
	rm -f *.d */*.d */*/*.d

veryclean:
	$(MAKE) clean
	rm -f lib*.a
ifeq ($(BLD_DIR),$(SRC_DIR))
	rm -f configIos/ConfigIoFactory.cpp
	rm -f analyzers/mcSystem/McAnalyzerFactory.cpp
	rm -f analyzers/mdSystem/MdAnalyzerFactory.cpp
	rm -f analyzers/system/SystemAnalyzerFactory.cpp
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
	rm -f species/SpeciesFactory.cpp
	rm -f trajectoryIos/TrajectoryIoFactory.cpp
	rm -f generators/generatorFactory.cpp
	rm -f mdSim.cpp
	rm -f mcSim.cpp
endif

# Executable targets

$(MCSIM_EXE): $(MCSIM).o $(LIBS)
	$(CXX) -o $(MCSIM_EXE) $(MCSIM).o $(LIBS) $(LDFLAGS)
 
$(MDSIM_EXE): $(MDSIM).o $(LIBS)
	$(CXX) -o $(MDSIM_EXE) $(MDSIM).o $(LIBS) $(LDFLAGS)

# Short names for executable targets (for convenience)
 
mcSim:
	$(MAKE) $(MCSIM_EXE)

mdSim:
	$(MAKE) $(MDSIM_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(util_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(MDSIM).d 
-include $(MCSIM).d 
