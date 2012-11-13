# ------------------------------------------------------------------
# Users may need to change the following Makefile macros:
#  - SRC_DIR_REL, if this directory is moved or copied to another location.
#  - BIN_DIR, to change the directory for executables

# Relative path to to Simpatico src directory.
SRC_DIR_REL=../..

# Include master makefiles

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk
include $(SRC_DIR_REL)/modules/hoomd/sources.mk
include $(SRC_DIR_REL)/modules/hoomd/defines.mk

#-----------------------------------------------------------------------
# Executable file name macros

BIN_DIR=$(ROOT_DIR)/bin
BIN_SUFFIX:=$(UTIL_MPI_SUFFIX)$(UTIL_SUFFIX)$(INTER_SUFFIX)$(MCMD_SUFFIX)_aspect_h

MCSIM_EXE=$(BIN_DIR)/mcSim$(BIN_SUFFIX)

# Rule to compile all NVCC source (*.cu) files.
%.cu.o: %.cu
	nvcc -arch sm_20 -c -o $@ $<

#-------------------------------------------------------------------
# Major targets

all:	
	make $(modules_hoomd_NVCC_OBJS)
	make $(modules_hoomd_OBJS)
	make mcSim.o $(MCSIM_EXE)

clean:
	-rm $(MCSIM_EXE) mcSim.o mcSim.d
	-rm -f $(modules_hoomd_OBJS) $(modules_hoomd_OBJS:.o=.d) $(modules_hoomd_NVCC_OBJS)

clean-deps:
	-rm -f $(modules_hoomd_OBJS:.o=.d) $(modules_hoomd_NVCC_OBJS:.o=.d) mcSim.d

-include $(modules_hoomd_OBJS:.o=.d) $(modules_hoomd_NVCC_OBJS:.o=.d)
-include mcSim.d

#-------------------------------------------------------------------
# Executables 

$(MCSIM_EXE): mcSim.o $(modules_hoomd_NVCC_OBJS) $(modules_hoomd_OBJS) $(modules_hoomd_NVCC_OBJS) $(util_LIB) $(inter_LIB) $(mcMd_LIB) $(LIBS)
	$(CXX) $(LDFLAGS) -o $(MCSIM_EXE) mcSim.o \
		 $(modules_hoomd_OBJS) $(modules_hoomd_NVCC_OBJS) $(mcMd_LIB) $(util_LIB) $(inter_LIB) $(LIBPATH) $(LIBS)

