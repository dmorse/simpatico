tools_tests_neighbor_= tools/tests/neighbor/Test.cc

tools_tests_neighbor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_tests_neighbor_))
tools_tests_neighbor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_tests_neighbor_:.cc=.o))

