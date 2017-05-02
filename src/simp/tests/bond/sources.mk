simp_tests_bond_=simp/tests/bond/Test.cc

simp_tests_bond_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_bond_))
simp_tests_bond_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_bond_:.cc=.o))

