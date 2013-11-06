SRC_DIR_REL =../../../

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/mcMd/config.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk
include sources.mk

all: $(mcMd_eeSimulation_OBJS)

clean:
	rm -f $(mcMd_eeSimulation_OBJS) $(mcMd_eeSimulation_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_eeSimulation_OBJS:.o=.d)

-include $(mcMd_eeSimulation_OBJS:.o=.d)

