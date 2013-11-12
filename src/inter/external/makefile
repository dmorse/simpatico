SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/inter/patterns.mk
include $(SRC_DIR_REL)/inter/external/sources.mk

all: $(inter_external_OBJS)

clean:
	rm -f $(inter_external_OBJS) $(inter_external_OBJS:.o=.d)

clean-deps:
	rm -f $(inter_external_OBJS:.o=.d)

-include $(inter_external_OBJS:.o=.d)

