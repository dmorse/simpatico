SRC_DIR_REL =..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk

ddSim=$(SRC_DIR)/ddMd/ddSim
ddSim_EXE=$(ROOT_DIR)/bin/ddSim

all: $(ddMd_OBJS) $(ddMd_LIB) $(ddSim).o $(ddSim_EXE)

clean:
	-rm -f $(ddMd_OBJS) $(ddMd_OBJS:.o=.d) $(ddMd_LIB)
	-rm -f $(ddSim).o $(ddSim).d
	-rm -f $(ddSim_EXE)

clean-deps:
	-rm -f $(ddMd_OBJS:.o=.d)
	-rm -f $(ddSim).d

-include $(ddMd_OBJS:.o=.d)
-include $(ddSim).d

$(ddSim_EXE): $(ddSim).o $(util_LIB) $(inter_LIB) $(ddMd_LIB)
	$(CXX) $(LDFLAGS) -o $(ddSim_EXE) $(ddSim).o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
 
