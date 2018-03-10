BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/tools/patterns.mk

all: $(tools_trajectory_OBJS)

clean:
	rm -f $(tools_trajectory_OBJS) $(tools_trajectory_OBJS:.o=.d)

clean-deps:
	rm -f $(tools_trajectory_OBJS:.o=.d)

-include $(tools_trajectory_OBJS:.o=.d)

