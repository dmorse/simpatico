SRC_DIR_REL =..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk

DDSIM=$(SRC_DIR)/ddMd/ddSim
#-----------------------------------------------------------------------

all: $(ddMd_OBJS) $(ddMd_LIB) $(DDSIM).o $(DDSIM_EXE)

clean:
	-rm -f $(ddMd_OBJS) $(ddMd_OBJS:.o=.d) $(ddMd_LIB)
	-rm -f $(DDSIM).o $(DDSIM).d
	-rm -f $(DDSIM_EXE)

clean-deps:
	-rm -f $(ddMd_OBJS:.o=.d)
	-rm -f $(DDSIM).d

-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(DDSIM).d

$(DDSIM_EXE): $(DDSIM).o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) -o $(DDSIM_EXE) $(DDSIM).o \
               $(ddMd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)

ddSim:
	make $(DDSIM_EXE)

