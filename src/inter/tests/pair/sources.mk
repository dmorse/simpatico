inter_tests_pair_=inter/tests/pair/Test.cc

inter_tests_pair_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_tests_pair_))
inter_tests_pair_OBJS=\
     $(addprefix $(BLD_DIR)/, $(inter_tests_pair_:.cc=.o))

