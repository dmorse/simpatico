SRC_DIR_REL =..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk

#-----------------------------------------------------------------------
MDSIM=$(SRC_DIR)/mcMd/mdSim
MCSIM=$(SRC_DIR)/mcMd/mcSim
BIN_DIR=$(ROOT_DIR)/bin
MCSIM_EXE=$(BIN_DIR)/mcSim$(UTIL_SUFFIX)$(INTER_SUFFIX)$(MCMD_SUFFIX)
MDSIM_EXE=$(BIN_DIR)/mdSim$(UTIL_SUFFIX)$(INTER_SUFFIX)$(MCMD_SUFFIX)
#-----------------------------------------------------------------------
all: $(mcMd_OBJS) $(mcMd_LIB) $(MDSIM).o $(MCSIM).o $(MCSIM_EXE) $(MDSIM_EXE)

clean:
	rm -f $(mcMd_OBJS) $(mcMd_LIB)
	rm -f $(mcMd_OBJS:.o=.d)
	rm -f $(MCSIM_EXE) $(MDSIM).o $(MDSIM).d
	rm -f $(MDSIM_EXE) $(MCSIM).o $(MCSIM).d
	rm -f Test Test.o Test.d

clean-deps:
	rm -f $(mcMd_OBJS:.o=.d)
	rm -f $(MDSIM).d
	rm -f $(MCSIM).d

clean-test:
	rm -f Test Test.o Test.d

#-----------------------------------------------------------------------
# Include dependency files

-include $(MDSIM).d 
-include $(MCSIM).d 
-include $(mcMd_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
#-----------------------------------------------------------------------
# Executables

main_EXES=$(MCSIM_EXE) $(MDSIM_EXE) 

$(MCSIM_EXE): $(MCSIM).o $(util_LIB) $(inter_LIB) $(mcMd_LIB)
	$(CXX) $(LDFLAGS) -o $(MCSIM_EXE) $(MCSIM).o \
               $(mcMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)
 
$(MDSIM_EXE): $(MDSIM).o $(util_LIB) $(inter_LIB) $(mcMd_LIB)
	$(CXX) $(LDFLAGS) -o $(MDSIM_EXE) $(MDSIM).o \
        $(mcMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)
 
mcSim:
	make $(MCSIM_EXE)

mdSim:
	make $(MDSIM_EXE)

#-----------------------------------------------------------------------
# Tests

Test: $(mcMd_LIB) $(SRC_DIR)/mcMd/Test.o
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test $(SRC_DIR)/mcMd/Test.o \
	       $(mcMd_LIB) $(util_LIB)

-include Test.d
