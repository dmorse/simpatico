# ---------------------------------------------------------------------
# File: src/inter/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at directory src/inter. This file
# is included by all makefile files in this directory tree. 
#
# This file should be included in other makefiles after inclusion of the 
# files src/config.mk, src/util/config.mk and src/inter/config.mk, 
# because this file uses makefile variables defined in those files.
#-----------------------------------------------------------------------

# All libraries needed in for files in src/inter
LIBS= $(inter_LIB) $(util_LIB)

# C preprocessor macro definitions needed in src/inter
DEFINES=$(UTIL_DEFS) $(INTER_DEFS)

# Dependencies on build configuration files
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/inter/config.mk

# Pattern rule to compile all *.cpp class source files in src/inter
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile all *.cpp test programs in src/inter/tests
$(BLD_DIR)/% $(BLD_DIR)/%.o:$(SRC_DIR)/%.cc $(LIBS)
	$(CXX) $(CPPFLAGS) $(TESTFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $(@:.o=) $@ $(LIBS)
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

