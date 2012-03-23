SRC_DIR_REL =..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/msDd/defines.mk

include $(SRC_DIR_REL)/msDd/patterns.mk

include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk
include $(SRC_DIR_REL)/msDd/sources.mk

#-----------------------------------------------------------------------

all: $(msDd_OBJS) $(msDd_LIB)

clean:
	-rm -f $(msDd_OBJS) $(msDd_OBJS:.o=.d)
	-rm -f  $(mdDd_LIB)

clean-deps:
	-rm -f $(msDd_OBJS:.o=.d)

-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(msDd_OBJS:.o=.d)

#$(DDSIM_EXE): $(DDSIM).o $(msDd_LIB) $(inter_LIB) $(util_LIB)
#	$(CXX) $(LDFLAGS) -o $(DDSIM_EXE) $(DDSIM).o \
#               $(msDd_LIB) $(inter_LIB) $(util_LIB) $(LIBS)

