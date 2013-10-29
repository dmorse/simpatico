util_tests_accumulators_unit_= \
    util/tests/accumulators/unit/Test.cpp

util_tests_accumulators_unit_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_accumulators_unit_))
util_tests_accumulators_unit_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(util_tests_accumulators_unit_:.cpp=.o))

