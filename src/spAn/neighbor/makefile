SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/spAn/config.mk
include $(SRC_DIR)/spAn/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/spAn/sources.mk

all: $(spAn_neighbor_OBJS)

clean:
	rm -f $(spAn_neighbor_OBJS) $(spAn_neighbor_OBJS:.o=.d)

clean-deps:
	rm -f $(spAn_neighbor_OBJS:.o=.d)

-include $(spAn_neighbor_OBJS:.o=.d)

