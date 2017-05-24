BLD_DIR_REL =../../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/tests/interaction/bond/sources.mk

all: $(simp_tests_interaction_bond_OBJS)

clean:
	rm -f $(simp_tests_interaction_bond_OBJS)
	rm -f $(simp_tests_interaction_bond_OBJS:.o=.d)
	rm -f $(simp_tests_interaction_bond_OBJS:.o=)
	rm -f *.rst
	rm -f out/*

-include $(simp_tests_interaction_bond_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

