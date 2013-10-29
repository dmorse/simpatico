OBJ_DIR_REL =../..

include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/misc/sources.mk

all: $(util_misc_OBJS)

clean:
	rm -f $(util_misc_OBJS) $(util_misc_OBJS:.o=.d)

clean-deps:
	rm -f $(util_misc_OBJS:.o=.d)

-include $(util_misc_OBJS:.o=.d)

