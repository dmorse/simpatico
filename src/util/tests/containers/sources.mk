util_tests_containers_=util/tests/containers/Test.cpp

util_tests_containers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_containers_))
util_tests_containers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_containers_:.cpp=.o))

