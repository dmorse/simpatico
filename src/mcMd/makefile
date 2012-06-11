SRC_DIR_REL =..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk

#-----------------------------------------------------------------------
# Variable definitions

MDSIM=$(SRC_DIR)/mcMd/mdSim
MCSIM=$(SRC_DIR)/mcMd/mcSim
MCMD_EXE_SUFFIX=$(UTIL_MPI_SUFFIX)$(UTIL_SUFFIX)$(INTER_SUFFIX)$(MCMD_SUFFIX)
MDSIM_EXE=$(BIN_DIR)/mdSim$(MCMD_EXE_SUFFIX)
MCSIM_EXE=$(BIN_DIR)/mcSim$(MCMD_EXE_SUFFIX)
MCMD_LIBS=$(mcMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)
#-----------------------------------------------------------------------

all: $(mcMd_OBJS) $(mcMd_LIB) $(MDSIM).o $(MCSIM).o $(MCSIM_EXE) $(MDSIM_EXE)

clean:
	rm -f $(mcMd_OBJS) $(mcMd_LIB)
	rm -f $(mcMd_OBJS:.o=.d)
	rm -f $(MCSIM_EXE) $(MDSIM).o $(MDSIM).d
	rm -f $(MDSIM_EXE) $(MCSIM).o $(MCSIM).d
	rm -f */*.o */*/*.o
	rm -f */*.d */*/*.d

clean-deps:
	rm -f $(mcMd_OBJS:.o=.d)
	rm -f $(MDSIM).d
	rm -f $(MCSIM).d
	rm -f */*.d */*/*.d

#-----------------------------------------------------------------------
# Include dependency files

-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(MDSIM).d 
-include $(MCSIM).d 
#-----------------------------------------------------------------------
# Executables

$(MCSIM_EXE): $(MCSIM).o $(MCMD_LIBS)
	$(CXX) $(LDFLAGS) -o $(MCSIM_EXE) $(MCSIM).o $(MCMD_LIBS) $(LIBS)
 
$(MDSIM_EXE): $(MDSIM).o $(MCMD_LIBS)
	$(CXX) $(LDFLAGS) -o $(MDSIM_EXE) $(MDSIM).o $(MCMD_LIBS)

# Targets without full path (for convenience)
 
mcSim:
	make $(MCSIM_EXE)

mdSim:
	make $(MDSIM_EXE)

