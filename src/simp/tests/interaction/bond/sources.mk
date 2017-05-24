simp_tests_interaction_bond_=\
     simp/tests/interaction/bond/Test.cc

simp_tests_interaction_bond_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_interaction_bond_))
simp_tests_interaction_bond_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_interaction_bond_:.cc=.o))

