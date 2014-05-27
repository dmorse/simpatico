# ---------------------------------------------------------------------
# File: src/mdCf/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at directory src/mdCf, which
# contains the source code for the DdPp namespace. It is included by
# all makefile files in this directory tree. 
#
# This file should be included in other makefiles after inclusion of
# the files src/config.mk, src/util/config.mk, src/inter/config.mk, 
# and src/ddMd/config.mk, because this file uses makefile variables
# defined in those files.
#-----------------------------------------------------------------------

# All libraries needed in this namespace
LIBS=$(mdCf_LIB) $(inter_LIB) $(util_LIB)

# All C preprocessor macro definitions for this namespace
DEFINES=$(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS)

# Dependencies of source files in src/ddMd on makefile fragments
MAKE_DEPS= -A$(OBJ_DIR)/config.mk
MAKE_DEPS+= -A$(OBJ_DIR)/util/config.mk
MAKE_DEPS+= -A$(OBJ_DIR)/inter/config.mk
MAKE_DEPS+= -A$(OBJ_DIR)/mdCf/config.mk

# Pattern rule to compile all class source (*.cpp) files in src/mdCf
$(OBJ_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(CXX_STD) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(OBJ_DIR) $<
endif

# Pattern rule to compile all test source (*.cc) files in src/mdCf
$(OBJ_DIR)/% $(OBJ_DIR)/%.o:$(SRC_DIR)/%.cc $(LIBS)
	$(CXX) $(CPPFLAGS) $(TESTFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $(@:.o=) $@ $(LIBS)
ifdef MAKEDEP
	$(MAKEDEP) $(CXX_STD) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(OBJ_DIR) $<
endif

