# ---------------------------------------------------------------------
# File: src/msDd/patterns.mk
#
# This makefile contains the pattern rules used to compile all sources
# files in the directory tree rooted at directory src/ddMd, which
# contains the source code for the DdMd namespace. It is included by
# all makefile files in this directory tree. 
#
# These patterns pass the compiler strings $(UTIL_DEFS) $(SIMP_DEFS)
# and $(DDMD_DEFS) that are defined in the config.mk files in the 
# util/, simp/, and mcMd/ directories, and in the src/config.mk
# file, and should thus be included after these four files. 
#-----------------------------------------------------------------------

# Path(s) to search for header files. 
INCLUDES= -I$(SRC_DIR)

# Extra dependencies for all source files
MSDD_ALLDEPS= -A$(SRC_DIR)/config.mk
MSDD_ALLDEPS+= -A$(SRC_DIR)/util/config.mk
MSDD_ALLDEPS+= -A$(SRC_DIR)/simp/config.mk
MSDD_ALLDEPS+= -A$(SRC_DIR)/mcMd/config.mk
MSDD_ALLDEPS+= -A$(SRC_DIR)/ddMd/config.mk

# Rule to compile all class source (*.cpp) files.
%.o:%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) \
            $(UTIL_DEFS) $(SIMP_DEFS) $(MCMD_DEFS) $(DDMD_DEFS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) \
            $(UTIL_DEFS) $(SIMP_DEFS) $(MCMD_DEFS) $(DDMD_DEFS) \
            $(MSDD_ALLDEPS) $<
endif

# Rule to compile *.cc main programs for unit tests. 
%.o:%.cc
	$(CXX) $(CPPFLAGS) $(TESTFLAGS) $(INCLUDES) \
            $(UTIL_DEFS) $(SIMP_DEFS) $(MCMD_DEFS) $(DDMD_DEFS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) \
            $(UTIL_DEFS) $(SIMP_DEFS) $(MCMD_DEFS) $(DDMD_DEFS) \
            $(MSDD_ALLDEPS) $<
endif

# Note: The main program files for unit tests must use a file suffix *.cc,
# while all other source files must use *.cpp. 
