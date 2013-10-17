# ---------------------------------------------------------------------
# File: src/mcMd/patterns.mk
#
# This makefile contains the pattern rules used to compile all sources
# files in the directory tree rooted at directory src/mcMd, which
# contains the source code for the McMd namespace. It is included by
# all makefile files in this directory tree. It should be included 
# after the file src/compiler.mk and after the defines.mk files from
# the src/util, src/inter, and src/mcMd subdirectories. 
#-----------------------------------------------------------------------
# Compilation pattern rules

# Path(s) to search for header files.
INCLUDES= -I$(SRC_DIR) -I$(TESTS_DIR)

# Dependencies of source files on makefile fragments
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
