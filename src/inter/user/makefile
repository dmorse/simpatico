SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/inter/patterns.mk
include $(SRC_DIR_REL)/inter/user/sources.mk

all: $(inter_user_OBJS)

clean:
	rm -f $(inter_user_OBJS) $(inter_user_OBJS:.o=.d)

clean-deps:
	rm -f $(inter_user_OBJS:.o=.d)

-include $(inter_user_OBJS:.o=.d)

