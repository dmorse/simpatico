BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

include sources.mk

all: $(mcMd_simulation_OBJS)

clean:
	rm -f $(mcMd_simulation_OBJS) $(mcMd_simulation_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_simulation_OBJS:.o=.d)

-include $(mcMd_simulation_OBJS:.o=.d)

