simp_tests_interaction_dihedral_=\
     simp/tests/interaction/dihedral/Test.cc

simp_tests_interaction_dihedral_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_interaction_dihedral_))
simp_tests_interaction_dihedral_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_interaction_dihedral_:.cc=.o))

