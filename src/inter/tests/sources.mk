inter_tests_=inter/tests/Test.cc

inter_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_tests_))
inter_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(inter_tests_:.cc=.o))

