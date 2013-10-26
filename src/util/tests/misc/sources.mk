util_tests_misc_=util/tests/misc/Test.cpp

util_tests_misc_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_misc_))
util_tests_misc_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(util_tests_misc_:.cpp=.o))

