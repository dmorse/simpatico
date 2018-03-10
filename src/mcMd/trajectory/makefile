BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_trajectory_OBJS)

clean:
	rm -f $(mcMd_trajectory_OBJS) $(mcMd_trajectory_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_trajectory_OBJS:.o=.d)

-include $(mcMd_trajectory_OBJS:.o=.d)

