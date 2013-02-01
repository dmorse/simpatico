SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include sources.mk

all: $(mcMd_misc_OBJS)

clean:
	rm -f $(mcMd_misc_OBJS) $(mcMd_misc_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_misc_OBJS:.o=.d)

-include $(mcMd_misc_OBJS:.o=.d)

