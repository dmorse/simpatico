SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/inter/patterns.mk
include $(SRC_DIR_REL)/inter/dihedral/sources.mk

all: $(inter_dihedral_OBJS)

clean:
	rm -f $(inter_dihedral_OBJS) $(inter_dihedral_OBJS:.o=.d)

clean-deps:
	rm -f $(inter_dihedral_OBJS:.o=.d)

-include $(inter_dihedral_OBJS:.o=.d)

