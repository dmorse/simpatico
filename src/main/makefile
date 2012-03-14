# ----------------------------------------------------------------------
# This makefile compiles and links programs for MC and MD simulations
# ----------------------------------------------------------------------

# Relative path to source file root directory
SRC_DIR_REL=..

include $(SRC_DIR_REL)/compiler.mk

include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/main/defines.mk

include $(SRC_DIR_REL)/mcMd/patterns.mk

include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk
include $(SRC_DIR_REL)/main/sources.mk

#-----------------------------------------------------------------------


all: $(main_OBJS) $(MCSIM_EXE) $(MDSIM_EXE)

clean:
	-rm -f $(main_OBJS) $(main_OBJS:.o=.d)
	-rm -f $(MCSIM_EXE) $(MDSIM_EXE)

clean-deps:
	-rm -f $(main_OBJS:.o=.d)

# Executables

-include $(main_OBJS:.o=.d)

$(MCSIM_EXE): mcSim.o $(util_LIB) $(inter_LIB) $(mcMd_LIB)
	$(CXX) $(LDFLAGS) -o $(MCSIM_EXE) mcSim.o \
               $(mcMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)
 
$(MDSIM_EXE): mdSim.o $(util_LIB) $(mcMd_LIB)
	$(CXX) $(LDFLAGS) -o $(MDSIM_EXE) mdSim.o \
        $(mcMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)
 
