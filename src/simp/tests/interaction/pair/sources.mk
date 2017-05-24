simp_tests_interaction_pair_=\
    simp/tests/interaction/pair/Test.cc

simp_tests_interaction_pair_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_interaction_pair_))
simp_tests_interaction_pair_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_interaction_pair_:.cc=.o))

