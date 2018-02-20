# ---------------------------------------------------------------------
# File: src/mcMd/patterns.mk
#
# This makefile contains the pattern rule used to compile all sources
# files in the directory tree rooted at directory src/mcMd, which
# contains the source code for the McMd namespace. It is included by
# all makefile files in this directory tree. It should be included
# after the main config.mk in the build directory, because it uses
# makefile variables that are defined in the main configuration file.
#-----------------------------------------------------------------------

include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/mcMd/sources.mk

# All libraries needed for files in src/mcMd
LIBS=$(mcMd_LIB) $(simp_LIB) $(util_LIB)

# C preprocessor macro definitions needed by files in src/mcMd
DEFINES=$(UTIL_DEFS) $(SIMP_DEFS) $(MCMD_DEFS)

# Dependencies on config.mk build configuration files
UTIL_CFGS= -A$(BLD_DIR)/config.mk
UTIL_CFGS+= -A$(BLD_DIR)/util/config.mk
SIMP_CFGS=$(UTIL_CFGS)
SIMP_CFGS+= -A$(BLD_DIR)/simp/config.mk
MCMD_CFGS=$(SIMP_CFGS)
MCMD_CFGS+= -A$(BLD_DIR)/mcMd/config.mk

# Pattern rule to compile all *.cpp class source files in src/mcMd
$(BLD_DIR)/mcMd/%.o: $(SRC_DIR)/mcMd/%.cpp
	$(CXX) $(INCLUDES) $(DEFINES) $(CXXFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(CXX_STD) $(MCMD_CFGS) -S$(SRC_DIR) -B$(BLD_DIR) $<
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

# Pattern rule to compile all *.cc test programs in src/mcMd/tests
$(BLD_DIR)/mcMd/tests/%.o: $(SRC_DIR)/mcMd/tests/%.cc 
	$(CXX) $(INCLUDES) $(DEFINES) $(TESTFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(DEFINES) $(CXX_STD) $(MCMD_CFGS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Pattern rule to link all *.cc test programs in src/mcMd/tests
$(BLD_DIR)/mcMd/tests/%: $(BLD_DIR)/mcMd/tests/%.o $(LIBS)
	$(CXX) $(INCLUDES) $(DEFINES) $(TESTFLAGS) -o $@ $@.o $(LIBS) $(LDFLAGS)

