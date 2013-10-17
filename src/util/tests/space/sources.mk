util_tests_space_=util/tests/space/Test.cpp

util_tests_space_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_space_))
util_tests_space_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_space_:.cpp=.o))

