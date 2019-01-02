mdPp_tests_neighbor_= mdPp/tests/neighbor/Test.cc

mdPp_tests_neighbor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_tests_neighbor_))
mdPp_tests_neighbor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mdPp_tests_neighbor_:.cc=.o))

