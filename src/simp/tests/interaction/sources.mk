simp_tests_interaction_=\
     simp/tests/interaction/Test.cc

simp_tests_interaction_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_interaction_))
simp_tests_interaction_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_interaction_:.cc=.o))

