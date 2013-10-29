SRC_DIR_REL =..

include $(SRC_DIR_REL)/compiler.mk
include $(OBJ_DIR)/util/defines.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/ddMd/defines.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/ddMd/sources.mk

#-----------------------------------------------------------------------
# Variable definitions

DDSIM=$(OBJ_DIR)/ddMd/ddSim
DDSIM_LIBS= $(ddMd_LIB) $(inter_LIB) $(util_LIB)

#-----------------------------------------------------------------------
# Main targets

all: $(ddMd_OBJS) $(ddMd_LIB) $(DDSIM).o $(DDSIM_EXE)

clean:
	-rm -f $(ddMd_OBJS) $(ddMd_OBJS:.o=.d) $(ddMd_LIB)
	-rm -f $(DDSIM).o $(DDSIM).d
	cd tests; make clean
	-rm -f *.o */*.o */*/*.o
	-rm -f *.d */*.d */*/*.d

veryclean:
	make clean
	-rm -f lib*.a
	-rm -f configIos/ConfigIoFactory.cpp
	-rm -f diagnostics/DiagnosticFactory.cpp
	-rm -f integrators/IntegratorFactory.cpp
	-rm -f potentials/angle/AngleFactory.cpp
	-rm -f potentials/bond/BondFactory.cpp
	-rm -f potentials/dihedral/DihedralFactory.cpp
	-rm -f potentials/external/ExternalFactory.cpp
	-rm -f potentials/pair/PairFactory.cpp

$(DDSIM_EXE): $(DDSIM).o $(DDSIM_LIBS)
	$(CXX) $(LDFLAGS) -o $(DDSIM_EXE) $(DDSIM).o $(DDSIM_LIBS)

ddSim:
	make $(DDSIM_EXE)

#-----------------------------------------------------------------------
# Include dependency files
 
-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(DDSIM).d

