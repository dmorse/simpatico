# ---------------------------------------------------------------------
# File: src/ddMd/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at directory src/ddMd, which
# contains the source code for the DdMd namespace. It is included by
# all makefile files in this directory tree. 
#
# This file should be included in other makefiles after inclusion of
# the files src/compiler.mk, src/util/defines.mk, src/inter/defines.mk, 
# and src/ddMd/defines.mk, because this file uses makefile variables
# defined in those files.
#-----------------------------------------------------------------------

# All libraries needed in this namespace
LIBS=$(ddMd_LIB) $(inter_LIB) $(util_LIB)

# All C preprocessor macro definitions for this namespace
DEFINES=$(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS)

# Dependencies of source files in src/ddMd on makefile fragments
MAKE_DEPS= -A$(OBJ_DIR)/compiler.mk
MAKE_DEPS+= -A$(OBJ_DIR)/util/defines.mk
MAKE_DEPS+= -A$(OBJ_DIR)/inter/defines.mk
MAKE_DEPS+= -A$(OBJ_DIR)/ddMd/defines.mk

# Pattern rule to compile all class source (*.cpp) files in src/ddMd
$(OBJ_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(OBJ_DIR) $<
endif

# Pattern rule to compile all test source (*.cc) files in src/ddMd
$(OBJ_DIR)/% $(OBJ_DIR)/%.o::$(SRC_DIR)/%.cc $(LIBS)
	$(CXX) $(CPPFLAGS) $(TESTFLAGS) $(INCLUDES) $(DEFINES) -c -o $@ $<
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o $(@:.o=) $< $(LIBS)
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(MAKE_DEPS) -S$(SRC_DIR) -B$(OBJ_DIR) $<
endif

