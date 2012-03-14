SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk

all: $(inter_dihedral_OBJS)

clean:
	rm -f $(inter_dihedral_OBJS) $(inter_dihedral_OBJS:.o=.d)

clean-deps:
	rm -f $(inter_dihedral_OBJS:.o=.d)

-include $(inter_dihedral_OBJS:.o=.d)

