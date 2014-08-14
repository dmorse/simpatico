spAn_tests_neighbor_= spAn/tests/neighbor/Test.cc

spAn_tests_neighbor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_tests_neighbor_))
spAn_tests_neighbor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_tests_neighbor_:.cc=.o))

