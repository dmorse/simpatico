util_tests_format_=util/tests/format/Test.cpp

util_tests_format_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_format_))
util_tests_format_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(util_tests_format_:.cpp=.o))

