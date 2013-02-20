SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/misc/sources.mk

all: $(util_misc_OBJS)

clean:
	rm -f $(util_misc_OBJS) $(util_misc_OBJS:.o=.d)

clean-deps:
	rm -f $(util_misc_OBJS:.o=.d)

-include $(util_misc_OBJS:.o=.d)

