# ---------------------------------------------------------------------
# File: src/mcMd/patterns.mk
#
# This makefile contains the pattern rules used to compile all sources
# files in the directory tree rooted at directory src/mcMd, which
# contains the source code for the McMd namespace. It is included by
# all makefile files in this directory tree.
#
# These patterns use the string $(MCMD_DEFS) of preprocessor macro
# definitions that is constructed in src/mcMc/defines.mk and various
# other variables defined in src/compiler.mk. It should thus be
# included in other makefiles after these files.
#-----------------------------------------------------------------------
# Compilation pattern rules

# Path(s) to search for header files.
INCLUDES= -I$(SRC_DIR)

# Extra dependencies for all source files
MCMD_ALLDEPS= -A$(SRC_DIR)/compiler.mk
MCMD_ALLDEPS+= -A$(SRC_DIR)/util/defines.mk
MCMD_ALLDEPS+= -A$(SRC_DIR)/inter/defines.mk
MCMD_ALLDEPS+= -A$(SRC_DIR)/mcMd/defines.mk

# Rule to compile all class source (*.cpp) files.
%.o:%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) \
               $(UTIL_DEFS) $(INTER_DEFS) $(MCMD_DEFS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) \
               $(UTIL_DEFS) $(INTER_DEFS) $(MCMD_DEFS) $(MCMD_ALLDEPS) $<
endif

# Rule to compile *.cc main programs for unit tests.
%.o:%.cc
	$(CXX) $(CPPFLAGS) $(TESTFLAGS) $(INCLUDES) \
               $(UTIL_DEFS) $(INTER_DEFS) $(MCMD_DEFS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) \
               $(UTIL_DEFS) $(INTER_DEFS) $(MCMD_DEFS) $(MCMD_ALLDEPS) $<
endif

# Note: The main program files for unit tests must use a file suffix *.cc,
# while all other source files must use *.cpp.

