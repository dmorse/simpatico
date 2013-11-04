OBJ_DIR_REL =..

include $(OBJ_DIR_REL)/compiler.mk
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

#-----------------------------------------------------------------------
# Main targets

all: $(ddMd_OBJS) $(ddMd_LIB) $(DDSIM).o $(DDSIM_EXE)

clean:
	-rm -f $(ddMd_OBJS) $(ddMd_OBJS:.o=.d) $(ddMd_LIB)
	-rm -f $(DDSIM).o $(DDSIM).d
	cd tests; make clean
	-rm -f */*.o */*/*.o
	-rm -f */*.d */*/*.d

uninstall:
	-rm -f $(DDSIM_EXE)

$(DDSIM_EXE): $(DDSIM).o $(LIBS)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $(DDSIM_EXE) $(DDSIM).o $(LIBS)

ddSim:
	make $(DDSIM_EXE)

#-----------------------------------------------------------------------
# Include dependency files
 
-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(DDSIM).d

