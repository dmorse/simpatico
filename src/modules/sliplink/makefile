# ------------------------------------------------------------------
# Users may need to change the following Makefile macros:
#  - SRC_DIR_REL, if this directory is moved or copied to another location.
#  - BIN_DIR, to change the directory for executables

# Path to directory containing Simpatico library source files
SRC_DIR_REL=../..

# Include master makefiles
include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk
include $(SRC_DIR_REL)/modules/sliplink/sources.mk

#-----------------------------------------------------------------------
# Executable file name macros

BIN_DIR=$(ROOT_DIR)/bin
BIN_SUFFIX:=$(MCMD_SUFFIX)_sl

MCSIM_EXE=$(BIN_DIR)/mcSim$(BIN_SUFFIX)
MDSIM_EXE=$(BIN_DIR)/mdSim$(BIN_SUFFIX)

#-------------------------------------------------------------------
# Major targets

all:	
	make $(modules_sliplink_OBJS)
	make $(MDSIM_EXE) 
	make $(MCSIM_EXE) 

clean:
	-rm $(MCSIM_EXE) mcSim.o mcSim.d
	-rm $(MDSIM_EXE) mdSim.o mdSim.d
	-rm -f $(modules_sliplink_OBJS) $(modules_sliplink_OBJS:.o=.d)

clean-deps:
	-rm -f $(modules_sliplink_OBJS:.o=.d)

-include $(modules_sliplink_OBJS:.o=.d)
-include mcSim.d
-include mdSim.d

#-------------------------------------------------------------------
# Executables 

$(MCSIM_EXE): mcSim.o $(modules_sliplink_OBJS) \
                 $(mcMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)
	$(CXX) $(LDFLAGS) -o $(MCSIM_EXE) mcSim.o \
                 $(modules_sliplink_OBJS) \
                 $(mcMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)

$(MDSIM_EXE): mdSim.o $(modules_sliplink_OBJS) \
                 $(mcMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)
	$(CXX) $(LDFLAGS) -o $(MDSIM_EXE) mdSim.o \
                 $(modules_sliplink_OBJS) \
                 $(mcMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)
