inter_tests_dihedral_=inter/tests/dihedral/Test.cpp

inter_tests_dihedral_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_tests_dihedral_))
inter_tests_dihedral_OBJS=\
     $(addprefix $(BLD_DIR)/, $(inter_tests_dihedral_:.cpp=.o))

