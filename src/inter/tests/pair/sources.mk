inter_tests_pair_=inter/tests/pair/Test.cpp

inter_tests_pair_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_tests_pair_))
inter_tests_pair_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(inter_tests_pair_:.cpp=.o))

