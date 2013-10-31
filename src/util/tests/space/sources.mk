util_tests_space_=util/tests/space/Test.cc

util_tests_space_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_space_))
util_tests_space_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(util_tests_space_:.cc=.o))

