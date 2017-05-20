# ---------------------------------------------------------------------
# File: src/mcMd/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at directory src/mcMd, which
# contains the source code for the McMd namespace. It is included by
# all makefile files in this directory tree. 
#
# This file should be included in other makefiles after inclusion of
# the files src/config.mk, src/util/config.mk, src/simp/config.mk, 
# and src/mcMd/config.mk, because this file uses makefile variables
# defined in those files. 
#-----------------------------------------------------------------------

# All libraries needed for files in src/mcMd
LIBS=$(mcMd_LIB) $(simp_LIB) $(util_LIB)

# C preprocessor macro definitions needed by files in src/mcMd
DEFINES=$(UTIL_DEFS) $(SIMP_DEFS) $(MCMD_DEFS)

# Dependencies on config.mk build configuration files
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/simp/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/mcMd/config.mk

# Pattern rule to compile all *.cpp class source files in src/mcMd
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile all *.cc test programs in src/mcMd/tests
$(BLD_DIR)/% $(BLD_DIR)/%.o:$(SRC_DIR)/%.cc $(LIBS)
	$(CXX) $(CPPFLAGS) $(TESTFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
	$(CXX) $(INCLUDES) $(DEFINES) -o $(@:.o=) $@ $(LIBS) $(LDFLAGS)
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

