SRC_DIR_REL =..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/simp/config.mk
include $(SRC_DIR_REL)/mcMd/config.mk
include $(SRC_DIR_REL)/ddMd/config.mk
include $(SRC_DIR_REL)/msDd/config.mk

include $(SRC_DIR_REL)/msDd/patterns.mk

include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/simp/sources.mk
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
-include $(simp_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(msDd_OBJS:.o=.d)

#$(DDSIM_EXE): $(DDSIM).o $(msDd_LIB) $(simp_LIB) $(util_LIB)
#	$(CXX) $(LDFLAGS) -o $(DDSIM_EXE) $(DDSIM).o \
#               $(msDd_LIB) $(simp_LIB) $(util_LIB) $(LIBS)

