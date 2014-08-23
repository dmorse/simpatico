SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/spAn/config.mk
include $(SRC_DIR_REL)/spAn/patterns.mk
include $(SRC_DIR_REL)/spAn/sources.mk

all: $(spAn_user_OBJS)

clean:
	rm -f $(spAn_user_OBJS) $(spAn_user_OBJS:.o=.d)

clean-deps:
	rm -f $(spAn_user_OBJS:.o=.d)

-include $(spAn_user_OBJS:.o=.d)

