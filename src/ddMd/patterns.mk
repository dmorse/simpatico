# ---------------------------------------------------------------------
# File: src/ddMd/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at directory src/ddMd, which
# contains the source code for the DdMd namespace. It is included by
# all makefile files in this directory tree. 
#
# This pattern uses makefile variables defined in src/compiler.mk and 
# the variables $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) constructed in
# the defines.mk files in the src/util, src/mcMd, and src/ddMd 
# directories. This file should thus be included in other makefiles 
# after src/compiler.mk and after these three defines.mk files. 
#-----------------------------------------------------------------------

# Dependencies of source files on makefile fragments
DDMD_ALLDEPS= -A$(SRC_DIR)/compiler.mk
DDMD_ALLDEPS+= -A$(SRC_DIR)/util/defines.mk
DDMD_ALLDEPS+= -A$(SRC_DIR)/inter/defines.mk
DDMD_ALLDEPS+= -A$(SRC_DIR)/ddMd/defines.mk

# Rule to compile all class source (*.cpp) files in src/ddMd
%.o:%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) \
               $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) \
               $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) $(DDMD_ALLDEPS) $<
endif

