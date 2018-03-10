BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_neighbor_OBJS)

clean:
	rm -f $(mcMd_neighbor_OBJS) $(mcMd_neighbor_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_neighbor_OBJS:.o=.d)

-include $(mcMd_neighbor_OBJS:.o=.d)

