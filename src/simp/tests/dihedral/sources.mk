simp_tests_dihedral_=simp/tests/dihedral/Test.cc

simp_tests_dihedral_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_dihedral_))
simp_tests_dihedral_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_dihedral_:.cc=.o))

