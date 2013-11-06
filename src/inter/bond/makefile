SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/inter/patterns.mk
include $(SRC_DIR_REL)/inter/bond/sources.mk

all: $(inter_bond_OBJS)

clean:
	rm -f $(inter_bond_OBJS) $(inter_bond_OBJS:.o=.d)

clean-deps:
	rm -f $(inter_bond_OBJS:.o=.d)

-include $(inter_bond_OBJS:.o=.d)

