util_tests_containers_=util/tests/containers/Test.cc

util_tests_containers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_containers_))
util_tests_containers_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(util_tests_containers_:.cc=.o))

