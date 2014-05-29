SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/mdCf/config.mk
include $(SRC_DIR_REL)/mdCf/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/mdCf/user/sources.mk

all: $(mdCf_user_OBJS)

clean:
	rm -f $(mdCf_user_OBJS) $(mdCf_user_OBJS:.o=.d)

clean-deps:
	rm -f $(mdCf_user_OBJS:.o=.d)

-include $(mdCf_user_OBJS:.o=.d)

