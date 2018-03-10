BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk


all: $(ddMd_neighbor_OBJS)

clean:
	rm -f $(ddMd_neighbor_OBJS) $(ddMd_neighbor_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_neighbor_OBJS:.o=.d)

-include $(ddMd_neighbor_OBJS:.o=.d)

