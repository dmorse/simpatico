util_tests_=util/tests/Test.cpp

util_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_))
util_tests_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(util_tests_:.cpp=.o))

