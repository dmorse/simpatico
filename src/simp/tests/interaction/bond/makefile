BLD_DIR_REL =../../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/tests/interaction/bond/sources.mk

all: $(simp_tests_interaction_bond_EXES)

clean:
	rm -f $(simp_tests_interaction_bond_EXES)
	rm -f $(simp_tests_interaction_bond_OBJS)
	rm -f $(simp_tests_interaction_bond_OBJS:.o=.d)
	rm -f *.rst
	rm -f out/*

-include $(simp_tests_interaction_bond_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

