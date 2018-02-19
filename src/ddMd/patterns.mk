# ---------------------------------------------------------------------
# File: src/ddMd/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at directory src/ddMd. This
# directory contains the source code for the DdMd namespace. This file
# is included by all makefile files in the src/ddMd directory tree. 
#
# This file should be included in other makefiles after inclusion of
# the files src/config.mk, src/util/config.mk, src/simp/config.mk, 
# and src/ddMd/config.mk, because pattern defined in this file uses 
# makefile variables defined in those other makefile fragments.
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

