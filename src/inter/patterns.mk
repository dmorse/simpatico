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

# All libraries needed in this namespace
LIBS= $(inter_LIB) $(util_LIB)

# C preprocessor macro definitions
DEFINES=$(UTIL_DEFS) $(INTER_DEFS)

# Dependencies of source files in src/inter on makefile fragments
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/inter/config.mk

# Pattern rule to compile all class source (*.cpp) files in src/inter
$(BLD_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile all class test (*.cpp) files in src/inter
$(BLD_DIR)/% $(BLD_DIR)/%.o:$(SRC_DIR)/%.cc $(LIBS)
	$(CXX) $(CPPFLAGS) $(TESTFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $(@:.o=) $@ $(LIBS)
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

