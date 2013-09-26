# ---------------------------------------------------------------------
# File: src/inter/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at directory src/inter. This file
# is included by all makefile files in this directory tree. 
#
# This file should be included in other makefiles after inclusion of the 
# files src/compiler.mk, src/util/defines.mk and src/inter/defines.mk, 
# because this file uses makefile variables defined in those files.
#-----------------------------------------------------------------------

# C preprocessor macro definitions
CPPDEFS=$(UTIL_DEFS) $(INTER_DEFS)

# Dependencies of source files in src/inter on makefile fragments
MAKE_DEPS= -A$(SRC_DIR)/compiler.mk
MAKE_DEPS+= -A$(SRC_DIR)/util/defines.mk
MAKE_DEPS+= -A$(SRC_DIR)/inter/defines.mk

# Pattern rule to compile all class source (*.cpp) files in src/inter
%.o:%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(CPPDEFS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(CPPDEFS) $(MAKE_DEPS) $<
endif

