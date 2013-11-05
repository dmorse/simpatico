SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/ddMd/config.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk

all: $(ddMd_potentials_pair_OBJS)

clean:
	rm -f $(ddMd_potentials_pair_OBJS) $(ddMd_potentials_pair_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_potentials_pair_OBJS:.o=.d)

-include $(ddMd_potentials_pair_OBJS:.o=.d)

