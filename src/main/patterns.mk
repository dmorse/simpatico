# ---------------------------------------------------------------------
# File: src/mcMd/patterns.mk
#
# This makefile contains the pattern rules used to compile all sources
# files in the directory tree rooted at directory src/mcMd, which
# contains the source code for the McMd namespace. It is included by
# all makefile files in this directory tree.
#-----------------------------------------------------------------------
# Compilation pattern rules

# Path(s) to search for header files.
INCLUDES= -I$(SRC_DIR)

# Extra dependencies for all source files
MAIN_ALLDEPS= -A$(SRC_DIR)/compiler.mk
MAIN_ALLDEPS+= -A$(SRC_DIR)/util/defines.mk
MAIN_ALLDEPS+= -A$(SRC_DIR)/inter/defines.mk
MAIN_ALLDEPS+= -A$(SRC_DIR)/mcMd/defines.mk

# Rule to compile all class source (*.cpp) files.
%.o:%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(UTIL_DEFS) \
               $(INTER_DEFS) $(MCMD_DEFS) $(DDMD_DEFS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) \
               $(MCMD_DEFS) $(DDMD_DEFS) $(MAIN_ALLDEPS) $<
endif

