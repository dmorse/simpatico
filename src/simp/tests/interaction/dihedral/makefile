BLD_DIR_REL =../../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/tests/interaction/dihedral/sources.mk

TEST= $(SRC_DIR)/simp/tests/interaction/dihedral/Test

all: $(simp_tests_interaction_dihedral_EXES)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(simp_tests_interaction_dihedral_EXES)
	rm -f $(simp_tests_interaction_dihedral_OBJS)
	rm -f $(simp_tests_interaction_dihedral_OBJS:.o=.d)

-include $(simp_tests_interaction_dihedral_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

