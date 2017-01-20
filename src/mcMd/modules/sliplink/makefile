# ------------------------------------------------------------------
# Users may need to change the following Makefile macros:
#  - SRC_DIR_REL, if this directory is moved or copied to another location.
#  - BIN_DIR, to change the directory for executables

# Path to directory containing Simpatico library source files
BLD_DIR_REL=../..

# Include master makefiles
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(SRC_DIR)/modules/sliplink/sources.mk

#-----------------------------------------------------------------------
# Executable file name macros

BIN_DIR=$(ROOT_DIR)/bin
BIN_SUFFIX:=$(MCMD_SUFFIX)_sl

MCSIM=$(BLD_DIR)/modules/sliplink/mcSim
MDSIM=$(BLD_DIR)/modules/sliplink/mdSim

MCSIM_EXE=$(BIN_DIR)/mcSim$(BIN_SUFFIX)
MDSIM_EXE=$(BIN_DIR)/mdSim$(BIN_SUFFIX)

#-------------------------------------------------------------------
# Major targets

all:	
	$(MAKE) $(modules_sliplink_OBJS)
	$(MAKE) $(MDSIM).o
	$(MAKE) $(MCSIM).o
	$(MAKE) $(MDSIM_EXE) 
	$(MAKE) $(MCSIM_EXE) 

clean:
	-rm -f $(MCSIM_EXE) mcSim.o mcSim.d
	-rm -f $(MDSIM_EXE) mdSim.o mdSim.d
	-rm -f $(modules_sliplink_OBJS) $(modules_sliplink_OBJS:.o=.d)

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
