# ---------------------------------------------------------------------
# File: src/util/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at the src/util directory, which
# contains all source code for the Util namespace. It is included by
# all "makefile" files in this directory tree. 
#
# This pattern use makefile variables defined in src/compiler.mk, and
# uses a variable $(UTIL_DEFS) that is defined in src/util/defines.mk.
# This file should thus be be included after src/compiler.mk and 
# src/util/defines.mk in other makefiles.
#-----------------------------------------------------------------------

# Dependencies of source files on makefile fragments
UTIL_ALLDEPS= -A$(SRC_DIR)/compiler.mk
UTIL_ALLDEPS+= -A$(SRC_DIR)/util/defines.mk

# Rule to compile *.cpp class source (*.cpp) files in src/util
%.o:%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(UTIL_DEFS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(UTIL_DEFS) $(UTIL_ALLDEPS) $<
endif

#-----------------------------------------------------------------------
