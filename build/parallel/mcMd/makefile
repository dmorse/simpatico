BLD_DIR_REL =..

include $(BLD_DIR_REL)/compiler.mk
include $(BLD_DIR)/util/defines.mk
include $(BLD_DIR)/inter/defines.mk
include $(BLD_DIR)/mcMd/defines.mk

include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mcMd/sources.mk

#-----------------------------------------------------------------------
# Variable definitions

MDSIM=$(BLD_DIR)/mcMd/mdSim
MCSIM=$(BLD_DIR)/mcMd/mcSim
MCMD_LIBS=$(mcMd_LIB) $(inter_LIB) $(util_LIB)

#-----------------------------------------------------------------------
# Main targets

all: $(mcMd_OBJS) $(mcMd_LIB) $(MDSIM).o $(MCSIM).o $(MCSIM_EXE) $(MDSIM_EXE)

lib: $(mcMd_OBJS) $(mcMd_LIB) 

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
# Executables

$(MCSIM_EXE): $(MCSIM).o $(MCMD_LIBS)
	$(CXX) $(LDFLAGS) -o $(MCSIM_EXE) $(MCSIM).o $(MCMD_LIBS) $(LIBS)
 
$(MDSIM_EXE): $(MDSIM).o $(MCMD_LIBS)
	$(CXX) $(LDFLAGS) -o $(MDSIM_EXE) $(MDSIM).o $(MCMD_LIBS)

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
