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

all: $(ddMd_actors_OBJS)

clean:
	rm -f $(ddMd_actors_OBJS) $(ddMd_actors_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_actors_OBJS:.o=.d)

-include $(ddMd_actors_OBJS:.o=.d)

