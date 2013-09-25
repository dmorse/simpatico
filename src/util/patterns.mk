# ---------------------------------------------------------------------
# File: src/util/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/util directory, which
# contains all source code for the Util namespace. It is included by
# all "makefile" files in this directory tree. 
#
# This file should be included in other makefiles after inclusion of
# the files src/compiler.mk and src/util/defines.mk because this file
# uses makefile variables defined in those files.
#-----------------------------------------------------------------------

# Preprocessor macro definitions
CPPDEFS=$(UTIL_DEFS)

# Dependencies of source files in src/util on makefile fragments
MAKE_DEPS= -A$(SRC_DIR)/compiler.mk
MAKE_DEPS+= -A$(SRC_DIR)/util/defines.mk

# Pattern rule to compile *.cpp class source (*.cpp) files in src/util
%.o:%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(CPPDEFS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(CPPDEFS) $(MAKE_DEPS) $<
endif

#-----------------------------------------------------------------------
