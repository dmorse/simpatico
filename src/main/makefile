# ----------------------------------------------------------------------
# This makefile compiles and links programs for MC and MD simulations
# ----------------------------------------------------------------------

# Relative path to source file root directory
SRC_DIR_REL=..

include $(SRC_DIR_REL)/compiler.mk

include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/main/defines.mk

include $(SRC_DIR_REL)/mcMd/patterns.mk

include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk
include $(SRC_DIR_REL)/main/sources.mk

#-----------------------------------------------------------------------


all: $(main_EXES)

clean:
	-rm -f $(main_OBJS) $(main_OBJS:.o=.d)
	-rm -f $(main_EXES)

clean-deps:
	-rm -f $(main_OBJS:.o=.d)

# Executables

$(MCSIM_EXE): $(MCSIM).o $(util_LIB) $(inter_LIB) $(mcMd_LIB)
	$(CXX) $(LDFLAGS) -o $(MCSIM_EXE) $(MCSIM).o \
               $(mcMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)
 
$(MDSIM_EXE): $(MDSIM).o $(util_LIB) $(inter_LIB) $(mcMd_LIB)
	$(CXX) $(LDFLAGS) -o $(MDSIM_EXE) $(MDSIM).o \
        $(mcMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)
 
$(DDSIM_EXE): $(DDSIM).o $(util_LIB) $(ddMd_LIB) $(inter_LIB)
	$(CXX) $(LDFLAGS) -o $(DDSIM_EXE) $(DDSIM).o \
        $(ddMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)

mcSim:
	make $(MCSIM_EXE)

mdSim:
	make $(MDSIM_EXE)

ddSim:
	make $(DDSIM_EXE)

# Object file targets and dependencies
-include $(main_OBJS:.o=.d)

