# ---------------------------------------------------------------------
# File: src/ddMd/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the src/ddMd directory tree, which contains the source code 
# for the DdMd C++ namespace. This file is included by all makefile 
# files in the src/ddMd directory tree. This file should be included 
# after inclusion of the main config.mk file.
#-----------------------------------------------------------------------
# Makefile includes 

include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/ddMd/config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/ddMd/sources.mk
#-----------------------------------------------------------------------

# All libraries needed by files in src/ddMd
LIBS=$(ddMd_LIB) $(simp_LIB) $(util_LIB)

# All C preprocessor macro definitions needed in src/ddMd
DEFINES=$(UTIL_DEFS) $(SIMP_DEFS) $(DDMD_DEFS)

# Dependencies of source files in src/ddMd on makefile fragments
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/simp/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/ddMd/config.mk

# Pattern rule to compile all *.cpp class source files in src/ddMd
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(INCLUDES) $(DEFINES) $(CXXFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(CXX_STD) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile all *.cc test programs in src/ddMd/tests
$(BLD_DIR)/ddMd/tests/%.o: $(SRC_DIR)/ddMd/tests/%.cc 
	$(CXX) $(INCLUDES) $(DEFINES) $(TESTFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(CXX_STD) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to link all *.cc test programs in src/ddMd/tests
$(BLD_DIR)/ddMd/tests/%: $(BLD_DIR)/ddMd/tests/%.o $(LIBS)
	$(CXX) $(INCLUDES) $(DEFINES) $(TESTFLAGS) -o $@ $@.o $(LIBS) $(LDFLAGS)

