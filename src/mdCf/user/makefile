SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/mdCf/config.mk
include $(SRC_DIR_REL)/mdCf/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/mdCf/user/sources.mk

all: $(ddMd_sp_user_OBJS)

clean:
	rm -f $(ddMd_sp_user_OBJS) $(ddMd_sp_user_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_sp_user_OBJS:.o=.d)

-include $(ddMd_sp_user_OBJS:.o=.d)

