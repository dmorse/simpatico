SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/simp/config.mk
include $(SRC_DIR_REL)/mcMd/config.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/mcMd/sources.mk

all: $(mcMd_modules_OBJS)

clean:
	rm -f $(mcMd_modules_OBJS) $(mcMd_modules_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_modules_OBJS:.o=.d)

-include $(mcMd_modules_OBJS:.o=.d)

