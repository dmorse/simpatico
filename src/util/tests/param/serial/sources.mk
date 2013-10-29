util_tests_param_serial_=util/tests/param/serial/Test.cpp

util_tests_param_serial_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_param_serial_))
util_tests_param_serial_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(util_tests_param_serial_:.cpp=.o))

