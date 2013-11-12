#-----------------------------------------------------------------------
# File includes

OBJ_DIR_REL =..

include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/mcMd/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mcMd/sources.mk

#-----------------------------------------------------------------------
# Variable definitions

MDSIM=$(OBJ_DIR)/mcMd/mdSim
MCSIM=$(OBJ_DIR)/mcMd/mcSim
#-----------------------------------------------------------------------
# Main targets

all: $(mcMd_OBJS) $(mcMd_LIB) $(MDSIM).o $(MCSIM).o $(MCSIM_EXE) $(MDSIM_EXE)

lib: $(mcMd_OBJS) $(mcMd_LIB) 

clean:
	-rm -f $(mcMd_OBJS) $(mcMd_OBJS:.o=.d) $(mcMd_LIB)
	-rm -f $(MDSIM).o $(MDSIM).d
	-rm -f $(MCSIM).o $(MCSIM).d
	cd tests; make clean
	-rm -f *.o */*.o */*/*.o
	-rm -f *.d */*.d */*/*.d

veryclean:
	make clean
	-rm -f lib*.a
	-rm -f configIos/ConfigIoFactory.cpp
	-rm -f analyzers/mcSystem/McAnalyzerFactory.cpp
	-rm -f analyzers/mdSystem/MdAnalyzerFactory.cpp
	-rm -f analyzers/system/SystemAnalyzerFactory.cpp
	-rm -f mcMoves/McMoveFactory.cpp
	-rm -f mdIntegrators/MdIntegratorFactory.cpp
	-rm -f perturb/mcSystem/McPerturbationFactory.cpp
	-rm -f potentials/angle/AngleFactory.cpp
	-rm -f potentials/bond/BondFactory.cpp
	-rm -f potentials/dihedral/DihedralFactory.cpp
	-rm -f potentials/external/ExternalFactory.cpp
	-rm -f potentials/link/LinkFactory.cpp
	-rm -f potentials/pair/PairFactory.cpp
	-rm -f species/SpeciesFactory.cpp
	-rm -f trajectoryIos/TrajectoryIoFactory.cpp

# Executable targets

$(MCSIM_EXE): $(MCSIM).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(MCSIM_EXE) $(MCSIM).o $(LIBS)
 
$(MDSIM_EXE): $(MDSIM).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(MDSIM_EXE) $(MDSIM).o $(LIBS)

# Short names for executable targets (for convenience)
 
mcSim:
	make $(MCSIM_EXE)

mdSim:
	make $(MDSIM_EXE)

#-----------------------------------------------------------------------
# Include dependency files

-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(MDSIM).d 
-include $(MCSIM).d 
