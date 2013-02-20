SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/space/sources.mk

all: $(util_space_OBJS)

clean:
	rm -f $(util_space_OBJS) $(util_space_OBJS:.o=.d)

clean-deps:
	rm -f $(util_space_OBJS:.o=.d)

-include $(util_space_OBJS:.o=.d)

