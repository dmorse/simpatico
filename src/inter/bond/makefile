SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/inter/patterns.mk
include $(SRC_DIR_REL)/inter/bond/sources.mk

all: $(inter_bond_OBJS)

clean:
	rm -f $(inter_bond_OBJS) $(inter_bond_OBJS:.o=.d)

clean-deps:
	rm -f $(inter_bond_OBJS:.o=.d)

-include $(inter_bond_OBJS:.o=.d)

