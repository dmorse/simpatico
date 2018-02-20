# ---------------------------------------------------------------------
# File: src/ddMd/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the src/ddMd directory tree, which contains the source code 
# for the DdMd C++ namespace. This file is included by all makefile 
# files in the src/ddMd directory tree. This file should be included 
# after inclusion of the main config.mk file.
#-----------------------------------------------------------------------
# Makefile includes 

include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/ddMd/config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/ddMd/sources.mk
#-----------------------------------------------------------------------


# Dependencies of source files in src/ddMd on makefile fragments
UTIL_CFGS= -A$(BLD_DIR)/config.mk
UTIL_CFGS+= -A$(BLD_DIR)/util/config.mk
SIMP_CFGS= $(UTIL_CFGS)
SIMP_CFGS+= -A$(BLD_DIR)/simp/config.mk
DDMD_CFGS= $(SIMP_CFGS)
DDMD_CFGS+= -A$(BLD_DIR)/ddMd/config.mk

# All libraries needed by files in src/ddMd
LIBS= $(ddMd_LIB) $(simp_LIB) $(util_LIB) 

# All C preprocessor macro definitions needed in src/ddMd
DEFINES=$(UTIL_DEFS) $(SIMP_DEFS) $(DDMD_DEFS)

# Pattern rule to compile all *.cpp class source files in src/ddMd
$(BLD_DIR)/ddMd/%.o: $(SRC_DIR)/ddMd/%.cpp
	$(CXX) $(INCLUDES) $(DEFINES) $(CXXFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(CXX_STD) $(DDMD_CFGS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile all *.cpp class source files in src/simp
$(BLD_DIR)/simp/%.o: $(SRC_DIR)/simp/%.cpp
	$(CXX) $(INCLUDES) $(UTIL_DEFS) $(SIMP_DEFS) $(CXXFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(UTIL_DEFS) $(SIMP_DEFS) $(CXX_STD) $(SIMP_CFGS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile all *.cpp class source files in src/util
$(BLD_DIR)/util/%.o: $(SRC_DIR)/util/%.cpp
	$(CXX) $(INCLUDES) $(UTIL_DEFS) $(CXXFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(UTIL_DEFS) $(CXX_STD) $(UTIL_CFGS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to compile all *.cc test programs in src/ddMd/tests
$(BLD_DIR)/ddMd/tests/%.o: $(SRC_DIR)/ddMd/tests/%.cc 
	$(CXX) $(INCLUDES) $(DEFINES) $(TESTFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(CXX_STD) $(DDMD_CFGS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to link all *.cc test programs in src/ddMd/tests
$(BLD_DIR)/ddMd/tests/%: $(BLD_DIR)/ddMd/tests/%.o $(LIBS)
	$(CXX) $(INCLUDES) $(DEFINES) $(TESTFLAGS) -o $@ $@.o $(LIBS) $(LDFLAGS)

