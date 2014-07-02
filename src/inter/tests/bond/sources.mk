inter_tests_bond_=inter/tests/bond/Test.cc

inter_tests_bond_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_tests_bond_))
inter_tests_bond_OBJS=\
     $(addprefix $(BLD_DIR)/, $(inter_tests_bond_:.cc=.o))

