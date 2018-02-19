# ---------------------------------------------------------------------
# File: src/tools/patterns.mk
#
# This makefile fragment contains the pattern rule used to compile all 
# C++ sources files in the src/tools directory tree. The src/tools 
# directory contains all the source code for the Tools C++ namespace. 
#
# This file is included by all makefile files in the tools/ directory.
# This pattern file should be included in other makefiles after 
# inclusion of the main config.mk file, and the namespace level config 
# files named util/config.mk, simp/config.mk, and tools/config.mk,
# because the patterns defined in this file use makefile variables
# defined in these configuration files.
#-----------------------------------------------------------------------
# Makefile includes

include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/tools/config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/tools/sources.mk
#-----------------------------------------------------------------------

# All libraries needed by files in src/tools
LIBS=$(tools_LIB) $(simp_LIB) $(util_LIB)

# All C preprocessor macro definitions needed in src/tools
DEFINES=$(UTIL_DEFS) $(SIMP_DEFS) $(TOOLS_DEFS)

# Dependencies of source files in src/tools on makefile fragments
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/simp/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/tools/config.mk

# Pattern rule to compile all *.cpp class source files in src/tools
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(CXX_STD) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile all *.cc test programs in src/tools
$(BLD_DIR)/tools/tests/%.o: $(SRC_DIR)/tools/tests/%.cc
	$(CXX) $(INCLUDES) $(DEFINES) $(TESTFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(CXX_STD) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile link *.cc test programs in src/tools
$(BLD_DIR)/tools/tests/%: $(BLD_DIR)/tools/tests/%.o $(LIBS)
	$(CXX) $(INCLUDES) $(DEFINES) -o $@ $< $(LIBS) $(LDFLAGS)

