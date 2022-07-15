# ----------------------------------------------------------------------
# This makefile compiles main programs for MC and MD simulations
# ----------------------------------------------------------------------

# Relative path to source file root directory
# Must be defined before including makefile.inc files
SRC_DIR=../../../src

include $(SRC_DIR)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/patterns.mk
include makefile.inc

#-----------------------------------------------------------------------
# Executable file name macros

EXTRACT_TOPOLOGY_EXE=$(BIN_DIR)/extract_topology$(BIN_SUFFIX)

#-----------------------------------------------------------------------

all:
	$(MAKE) $(EXTRACT_TOPOLOGY_EXE) 

clean:	
	rm -f $(EXTRACT_TOPOLOGY_EXE) 

# Executables

$(EXTRACT_TOPOLOGY_EXE): $(EXTRACT_TOPOLOGY_SRC)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o $(EXTRACT_TOPOLOGY_EXE) $(EXTRACT_TOPOLOGY_CPP)
 
