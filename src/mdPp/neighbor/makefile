BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mdPp/patterns.mk

all: $(mdPp_neighbor_OBJS)

clean:
	rm -f $(mdPp_neighbor_OBJS) $(mdPp_neighbor_OBJS:.o=.d)

clean-deps:
	rm -f $(mdPp_neighbor_OBJS:.o=.d)

-include $(mdPp_neighbor_OBJS:.o=.d)

