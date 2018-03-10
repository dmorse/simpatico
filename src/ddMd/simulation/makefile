SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk

all: $(ddMd_simulation_OBJS)

clean:
	rm -f $(ddMd_simulation_OBJS) $(ddMd_simulation_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_simulation_OBJS:.o=.d)

-include $(ddMd_simulation_OBJS:.o=.d)

