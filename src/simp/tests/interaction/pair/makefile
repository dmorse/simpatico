BLD_DIR_REL =../../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/tests/interaction/pair/sources.mk

TEST=$(SRC_DIR)/simp/tests/interaction/pair/Test


all: $(simp_tests_interaction_pair_EXES) 

run: $(TEST)
	$(TEST)

clean:
	rm -f $(simp_tests_interaction_pair_EXES)
	rm -f $(simp_tests_interaction_pair_OBJS)
	rm -f $(simp_tests_interaction_pair_OBJS:.o=.d)
	$(MAKE) clean-outputs

clean-outputs:
	rm -f out/*

-include $(simp_tests_interaction_pair_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

