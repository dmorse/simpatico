# ---------------------------------------------------------------------
# File: src/ddMd/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at directory src/ddMd, which
# contains the source code for the DdMd namespace. It is included by
# all makefile files in this directory tree. 
#
# This file should be included in other makefiles after inclusion of
# the files src/config.mk, src/util/config.mk, src/simp/config.mk, 
# and src/ddMd/config.mk, because this file uses makefile variables
# defined in those files.
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
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(CXX_STD) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile all *.cc test programs in src/ddMd
$(BLD_DIR)/% $(BLD_DIR)/%.o:$(SRC_DIR)/%.cc $(LIBS)
	$(CXX) $(CPPFLAGS) $(TESTFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $(@:.o=) $@ $(LIBS)
ifdef MAKEDEP
	$(MAKEDEP) $(CXX_STD) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

