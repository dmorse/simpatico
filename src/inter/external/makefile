SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/inter/patterns.mk
include $(SRC_DIR_REL)/inter/external/sources.mk

all: $(inter_external_OBJS)

clean:
	rm -f $(inter_external_OBJS) $(inter_external_OBJS:.o=.d)

clean-deps:
	rm -f $(inter_external_OBJS:.o=.d)

-include $(inter_external_OBJS:.o=.d)

