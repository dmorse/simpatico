BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/tools/patterns.mk

all: $(tools_neighbor_OBJS)

clean:
	rm -f $(tools_neighbor_OBJS) $(tools_neighbor_OBJS:.o=.d)

clean-deps:
	rm -f $(tools_neighbor_OBJS:.o=.d)

-include $(tools_neighbor_OBJS:.o=.d)

