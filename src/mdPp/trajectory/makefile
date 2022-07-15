BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mdPp/patterns.mk

all: $(mdPp_trajectory_OBJS)

clean:
	rm -f $(mdPp_trajectory_OBJS) $(mdPp_trajectory_OBJS:.o=.d)

clean-deps:
	rm -f $(mdPp_trajectory_OBJS:.o=.d)

-include $(mdPp_trajectory_OBJS:.o=.d)

