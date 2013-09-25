# ---------------------------------------------------------------------
# File: src/mcMd/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at directory src/mcMd, which
# contains the source code for the McMd namespace. It is included by
# all makefile files in this directory tree. 
#
# This file should be included in other makefiles after inclusion of
# the file src/compiler.mk and after inclusion of the defines.mk files 
# from the src/util, src/inter, and src/mcMd subdirectories. 
#-----------------------------------------------------------------------

# Dependencies of source files on makefile fragments
MCMD_ALLDEPS= -A$(SRC_DIR)/compiler.mk
MCMD_ALLDEPS+= -A$(SRC_DIR)/util/defines.mk
MCMD_ALLDEPS+= -A$(SRC_DIR)/inter/defines.mk
MCMD_ALLDEPS+= -A$(SRC_DIR)/mcMd/defines.mk

# Rule to compile all class source (*.cpp) files in src/mcMd
%.o:%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) \
               $(UTIL_DEFS) $(INTER_DEFS) $(MCMD_DEFS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) \
               $(UTIL_DEFS) $(INTER_DEFS) $(MCMD_DEFS) $(MCMD_ALLDEPS) $<
endif

