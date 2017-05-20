#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Makefile includes 

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/ddMd/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/ddMd/sources.mk

#-----------------------------------------------------------------------
# Variable definition

# Main program *.o and *.d location
DDSIM=$(BLD_DIR)/ddMd/ddSim

#-----------------------------------------------------------------------
# Main targets

all: $(ddMd_OBJS) $(ddMd_LIB) $(DDSIM).o $(DDSIM_EXE)

$(DDSIM_EXE): $(DDSIM).o $(LIBS)
	$(CXX) -o $(DDSIM_EXE) $(DDSIM).o $(LIBS) $(LDFLAGS)

ddSim:
	$(MAKE) $(DDSIM_EXE)

clean:
	rm -f $(ddMd_OBJS) $(ddMd_OBJS:.o=.d) $(ddMd_LIB)
	rm -f $(DDSIM).o $(DDSIM).d
	cd tests; $(MAKE) clean
	rm -f *.o */*.o */*/*.o
	rm -f *.d */*.d */*/*.d

veryclean:
	$(MAKE) clean
	rm -f lib*.a
ifeq ($(BLD_DIR),$(SRC_DIR))
	rm -f configIos/ConfigIoFactory.cpp
	rm -f analyzers/AnalyzerFactory.cpp
	rm -f integrators/IntegratorFactory.cpp
	rm -f potentials/angle/AngleFactory.cpp
	rm -f potentials/bond/BondFactory.cpp
	rm -f potentials/dihedral/DihedralFactory.cpp
	rm -f potentials/external/ExternalFactory.cpp
	rm -f potentials/pair/PairFactory.cpp
endif

#-----------------------------------------------------------------------
# Include dependency files
 
-include $(util_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(DDSIM).d

