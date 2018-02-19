BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk

all: $(simp_interaction_dihedral_OBJS)

clean:
	rm -f $(simp_interaction_dihedral_OBJS) $(simp_interaction_dihedral_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_interaction_dihedral_OBJS:.o=.d)

-include $(simp_interaction_dihedral_OBJS:.o=.d)

