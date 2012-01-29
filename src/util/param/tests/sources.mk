include $(SRC_DIR)/util/param/tests/serial/sources.mk
include $(SRC_DIR)/util/param/tests/mpi/sources.mk

util_param_tests_SRCS=$(util_param_tests_serial_SRCS) \
    $(util_param_tests_mpi_SRCS) 

util_param_tests_OBJS=$(util_param_tests_SRCS:.cc=.o)

