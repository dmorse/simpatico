OBJ_DIR_REL =..

include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/util/defines.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/mcMd/defines.mk
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
	rm -f $(mcMd_OBJS) $(mcMd_LIB)
	rm -f $(mcMd_OBJS:.o=.d)
	rm -f $(MDSIM).o $(MDSIM).d
	rm -f $(MCSIM).o $(MCSIM).d
	cd tests; make clean
	-rm -f */*.o */*/*.o
	-rm -f */*.d */*/*.d

uninstall:
	rm -f $(MCSIM_EXE)
	rm -f $(MDSIM_EXE) 

#-----------------------------------------------------------------------
# Executables

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
