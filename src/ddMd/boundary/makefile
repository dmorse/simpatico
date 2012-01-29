SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk

all: $(ddMd_boundary_OBJS)

clean:
	rm -f $(ddMd_boundary_OBJS) $(ddMd_boundary_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_boundary_OBJS:.o=.d)

-include $(ddMd_boundary_OBJS:.o=.d)

