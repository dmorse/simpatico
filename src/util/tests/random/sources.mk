util_tests_random_=util/tests/random/Test.cpp

util_tests_random_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_random_))
util_tests_random_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(util_tests_random_:.cpp=.o))

