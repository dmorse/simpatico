SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/ddMd/config.mk

include $(SRC_DIR_REL)/msDd/patterns.mk

include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk
include $(SRC_DIR_REL)/msDd/sources.mk

#-----------------------------------------------------------------------

all: $(msDd_simulation_OBJS) 

clean:
	-rm -f $(msDd_simulation_OBJS) $(msDd_simulation_OBJS:.o=.d)

clean-deps:
	-rm -f $(msDd_simulation_OBJS:.o=.d)

-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(msDd_simulation_OBJS:.o=.d)


