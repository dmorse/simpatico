SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/simp/config.mk
include $(SRC_DIR_REL)/tools/config.mk
include $(SRC_DIR_REL)/tools/patterns.mk
include $(SRC_DIR_REL)/tools/sources.mk

all: $(tools_user_OBJS)

clean:
	rm -f $(tools_user_OBJS) $(tools_user_OBJS:.o=.d)

clean-deps:
	rm -f $(tools_user_OBJS:.o=.d)

-include $(tools_user_OBJS:.o=.d)

