inter_tests_dihedral_=inter/tests/dihedral/Test.cc

inter_tests_dihedral_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_tests_dihedral_))
inter_tests_dihedral_OBJS=\
     $(addprefix $(BLD_DIR)/, $(inter_tests_dihedral_:.cc=.o))

