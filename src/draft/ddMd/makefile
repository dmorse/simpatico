SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk
include sources.mk

all: $(draft_ddMd_OBJS)

clean:
	rm -f $(draft_ddMd_OBJS) $(draft_ddMd_OBJS:.o=.d)

clean-deps:
	rm -f $(draft_ddMd_OBJS:.o=.d)

-include $(draft_ddMd_OBJS:.o=.d)

