#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the obj/serial and obj/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Makefile includes 

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/ddMd/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/ddMd/sources.mk

#-----------------------------------------------------------------------
# Variable definition

# Main program location
DDSIM=$(OBJ_DIR)/ddMd/ddSim

# MD postprocessor program location
MDPP=$(OBJ_DIR)/ddMd/mdPp

#-----------------------------------------------------------------------
# Main targets

all: $(ddMd_OBJS) $(ddMd_LIB) $(DDSIM).o $(DDSIM_EXE) $(MDPP).o $(MDPP_EXE)

$(DDSIM_EXE): $(DDSIM).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(DDSIM_EXE) $(DDSIM).o $(LIBS)

ddSim:
	make $(DDSIM_EXE)

$(MDPP_EXE): $(MDPP).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(MDPP_EXE) $(MDPP).o $(LIBS)

mdPp:
	make $(MDPP_EXE)

clean:
	rm -f $(ddMd_OBJS) $(ddMd_OBJS:.o=.d) $(ddMd_LIB)
	rm -f $(DDSIM).o $(DDSIM).d
	rm -f $(MDPP).o $(MDPP).d
	cd tests; make clean
	rm -f *.o */*.o */*/*.o
	rm -f *.d */*.d */*/*.d

veryclean:
	make clean
	rm -f lib*.a
ifeq ($(OBJ_DIR),$(SRC_DIR))
	rm -f configIos/ConfigIoFactory.cpp
	rm -f analyzers/AnalyzerFactory.cpp
	rm -f integrators/IntegratorFactory.cpp
	rm -f potentials/angle/AngleFactory.cpp
	rm -f potentials/bond/BondFactory.cpp
	rm -f potentials/dihedral/DihedralFactory.cpp
	rm -f potentials/external/ExternalFactory.cpp
	rm -f potentials/pair/PairFactory.cpp
	rm -f sp/configIos/SpConfigIoFactory.cpp
	rm -f sp/analyzers/SpAnalyzerFactory.cpp
endif

#-----------------------------------------------------------------------
# Include dependency files
 
-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(DDSIM).d

