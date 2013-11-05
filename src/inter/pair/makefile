SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/inter/patterns.mk
include $(SRC_DIR_REL)/inter/pair/sources.mk

all: $(inter_pair_OBJS)

clean:
	rm -f $(inter_pair_OBJS) $(inter_pair_OBJS:.o=.d)

clean-deps:
	rm -f $(inter_pair_OBJS:.o=.d)

-include $(inter_pair_OBJS:.o=.d)

