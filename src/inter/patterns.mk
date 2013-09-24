# ---------------------------------------------------------------------
# File: src/inter/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at directory src/inter. This file
# is included by all makefile files in this directory tree. 
#
# These patterns use variables defined in src/compiler.mk, and the 
# strings $(INTER_UTIL) and $(INTER_DEFS) that are defined in the
# defines.mk files in src/util and src/inter directories.  This file 
# should thus be 
# included in other makefiles after these files. 
#-----------------------------------------------------------------------

# Dependencies of source files on makefile fragments
INTER_ALLDEPS= -A$(SRC_DIR)/compiler.mk
INTER_ALLDEPS+= -A$(SRC_DIR)/util/defines.mk
INTER_ALLDEPS+= -A$(SRC_DIR)/inter/defines.mk

# Rule to compile all class source (*.cpp) files in src/inter
%.o:%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(INTER_ALLDEPS) $<
endif

